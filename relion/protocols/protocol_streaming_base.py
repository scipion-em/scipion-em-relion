# **************************************************************************
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import json

import pyworkflow.protocol.constants as cons
import os


class _StepArgScan:
    """Incremental scan state for one ``_collectStepArgKeys`` query.

    Steps are only ever appended to the step lists, and a FINISHED step
    never goes back, so a scan can resume where the previous one stopped
    instead of re-parsing the whole graph on every poll.
    """

    def __init__(self):
        self.keys = set()
        self.cursors = {}
        # Positions, not step objects: the step lists can be rebuilt, and a
        # held-on-to object would keep reporting the status it had back then.
        self.pending = set()


class RelionStreamingBase:
    """Shared streaming helpers for Relion SPA protocols.

    The persisted protocol step graph is the completion authority for every
    streaming protocol here (no DONE marker files), so the logic that reads
    it lives in one place instead of being repeated per protocol.

    Everything here is built to cost O(new items) per poll rather than
    O(items seen so far): a streaming protocol can accumulate millions of
    items, and a full rescan on every poll would make it slower the longer
    it runs.
    """

    # --------------------------- termination ---------------------------

    def _streamingMustStop(self):
        """True when a streaming loop has to give up polling.

        A failed step makes pyworkflow mark the protocol as FAILED and
        the executor break out of its own loop - and then join every
        running thread. A loop that keeps polling is never joined, so the
        run hangs with nothing left to do. The same applies once it has
        been aborted.
        """
        status = getattr(self, 'status', None)
        value = status.get() if hasattr(status, 'get') else status

        return value in (cons.STATUS_FAILED, cons.STATUS_ABORTED)

    # ------------------------ per-item artefacts -----------------------

    def _itemScopedName(self, item, baseName):
        """A per-item artefact name two items can never share.

        A Set can hold two movies whose files differ only in their
        directory. Named after the basename alone they collide: one
        output overwrites the other, and both end up pointing at the
        survivor. The id keeps them apart and the original basename stays
        in the name so logs remain readable.
        """
        return '%06d__%s' % (item.getObjId(), baseName)

    def _itemScopedPath(self, item, baseName, pathFunc=None):
        """Path for a per-item artefact, scoped by the item's own id.

        Projects written before this was scoped hold the unscoped name,
        and those files are the user's results: when only the old one is
        on disk it is still the one returned, so Continue keeps working.
        """
        pathFunc = pathFunc or self._getExtraPath
        scoped = pathFunc(self._itemScopedName(item, baseName))

        if not os.path.exists(scoped):
            legacy = pathFunc(baseName)

            if os.path.exists(legacy):
                return legacy

        return scoped

    # ------------------------- input discovery -------------------------
    def _discoverIdsAfter(self, inputSet, lastId):
        """Discover logical ids above the current streaming watermark.

        The ``id > N`` filter is pushed down to the backend so an input Set
        that already holds millions of rows is not walked again on every
        poll. Backends that cannot filter fall back to a full id listing,
        which still beats hydrating every item.
        """
        try:
            ids = list(inputSet.getUniqueValues('id', where='id > %d' % lastId))
        except (NotImplementedError, TypeError):
            ids = [itemId for itemId in inputSet.getUniqueValues('id')
                   if itemId > lastId]

        ids = sorted(ids)

        if ids:
            lastId = max(ids)

        return ids, lastId

    # How many polls a closed producer may keep showing exactly the same
    # incomplete view before the protocol gives up on it.
    TERMINAL_STALL_POLLS = 10

    def _hasActiveStreamingWork(self):
        """Whether something is still in flight for this protocol.

        Work in flight is progress, however long it takes, so it must
        never be counted towards a terminal stall. A protocol that does
        not know how to answer this says so by returning True, which
        simply means the stall detector stays out of its way.
        """
        return True

    def _recordTerminalProgress(self, inputSet, knownIds, watermarkAttr,
                                terminalConsistent):
        """Refuse to poll forever for rows that are never coming.

        A producer can close declaring more items than the consumer can
        see, and usually the rest turn up a moment later. When they do
        not - the declared size, what is known, and the watermark all
        stay exactly as they were, poll after poll, with nothing in
        flight - the protocol would otherwise sit there RUNNING for the
        rest of time. Say what is missing and fail instead.
        """
        if terminalConsistent:
            self._terminalStallSignature = None
            self._terminalStallCount = 0

            return

        if self._hasActiveStreamingWork():
            self._terminalStallCount = 0

            return

        signature = (inputSet.getSize(), len(knownIds),
                     getattr(self, watermarkAttr, 0))

        if signature == getattr(self, '_terminalStallSignature', None):
            self._terminalStallCount = getattr(
                self, '_terminalStallCount', 0) + 1
        else:
            self._terminalStallSignature = signature
            self._terminalStallCount = 1

        if self._terminalStallCount >= self.TERMINAL_STALL_POLLS:
            raise RuntimeError(
                "The input stream closed declaring %d items but only %d "
                "are visible, and that has not changed in %d polls with "
                "nothing left to process. Refusing to wait for rows that "
                "are not coming."
                % (inputSet.getSize(), len(knownIds),
                   self._terminalStallCount))

    def _reconcileClosedStreamIds(self, inputSet, discoveredIds, knownIds,
                                  producerClosed, watermarkAttr):
        """Recover late-visible ids, but only once the producer is closed.

        A producer can declare itself closed while some of its rows are
        not visible yet. Rescanning every poll to cover that would defeat
        the watermark, so the full listing happens only in this terminal
        reconciliation, and only while the declared size still exceeds
        what has actually been seen.
        """
        discoveredIds = list(discoveredIds)

        if not producerClosed:
            return discoveredIds, True

        expectedSize = inputSet.getSize()
        knownIds = set(knownIds)
        visibleKnownIds = knownIds.union(discoveredIds)

        if len(visibleKnownIds) >= expectedSize:
            self._recordTerminalProgress(inputSet, visibleKnownIds,
                                         watermarkAttr, True)

            return discoveredIds, True

        reconciledIds = list(inputSet.getUniqueValues('id'))

        if reconciledIds:
            setattr(self, watermarkAttr,
                    max(getattr(self, watermarkAttr, 0), max(reconciledIds)))

        visibleIds = set(discoveredIds)
        visibleIds.update(reconciledIds)

        newIds = [itemId for itemId in sorted(visibleIds)
                  if itemId not in knownIds]

        reconciledKnownIds = knownIds.union(visibleIds)
        terminalConsistent = len(reconciledKnownIds) >= expectedSize

        self._recordTerminalProgress(inputSet, reconciledKnownIds,
                                     watermarkAttr, terminalConsistent)

        return newIds, terminalConsistent

    @staticmethod
    def _refreshLogicalSet(inputSet):
        """Reload the Set's own properties, once per poll."""
        loadAllProperties = getattr(inputSet, 'loadAllProperties', None)

        if callable(loadAllProperties):
            loadAllProperties()

    def _loadLogicalSetItemsByIds(self, inputSet, itemIds, batchSize=500,
                                  refresh=True):
        """Load only the items behind ``itemIds``, never the whole Set."""
        itemIds = sorted(set(itemIds))

        if not itemIds:
            return []

        if refresh:
            self._refreshLogicalSet(inputSet)

        def cloneItem(item):
            clone = getattr(item, 'clone', None)
            return clone() if callable(clone) else item

        iterItems = getattr(inputSet, 'iterItems', None)

        if callable(iterItems):
            items = []

            try:
                for offset in range(0, len(itemIds), batchSize):
                    batch = itemIds[offset:offset + batchSize]
                    where = 'id IN (%s)' % ','.join(str(i) for i in batch)

                    for item in iterItems(orderBy='id', direction='ASC',
                                          where=where):
                        items.append(cloneItem(item))

                return items
            except (NotImplementedError, TypeError):
                pass

        getItem = getattr(inputSet, 'getItem', None)

        if callable(getItem):
            items = []

            for itemId in itemIds:
                try:
                    item = getItem('id', itemId)
                except (UnboundLocalError, KeyError):
                    item = None

                if item is not None:
                    items.append(cloneItem(item))

            return items

        wantedIds = set(itemIds)

        return [cloneItem(item) for item in inputSet
                if item.getObjId() in wantedIds]

    def _resumeWatermarkWithGaps(self, inputSet, processedIds):
        """Where to restart discovery after a Continue, plus what it skips.

        Starting a resumed run at the highest already-processed id keeps the
        first poll from hydrating everything a previous run dealt with. On
        its own that would silently drop any lower id that was never
        processed - one that arrived late, or whose step failed - so the
        ids below the watermark are listed once (ids only, nothing
        hydrated) and the unprocessed ones are handed back as gaps.
        """
        processedIds = set(processedIds)

        if not processedIds:
            return 0, set()

        watermark = max(processedIds)

        try:
            seenIds = set(inputSet.getUniqueValues(
                'id', where='id <= %d' % watermark))
        except (NotImplementedError, TypeError):
            seenIds = {itemId for itemId in inputSet.getUniqueValues('id')
                       if itemId <= watermark}

        return watermark, seenIds - processedIds

    def _getKnownStreamIds(self, watermarkAttr):
        """Ids already seen on the stream tracked by ``watermarkAttr``.

        A protocol can watch several inputs at once - coordinates,
        micrographs, CTFs - and each advances on its own, so each gets its
        own watermark and its own set of seen ids.
        """
        allKnown = getattr(self, '_knownStreamIds', None)

        if allKnown is None:
            allKnown = {}
            self._knownStreamIds = allKnown

        return allKnown.setdefault(watermarkAttr, set())

    def _discoverNewInputItems(self, inputSet, watermarkAttr, knownIds):
        """Discover and load the items added since the last poll.

        Returns ``(items, producerClosed, terminalConsistent)``.
        """
        # Refresh the Set's own properties first: isStreamClosed() reads
        # one of them, and a stale value would keep the protocol polling
        # forever after the producer has actually closed. Once here is
        # enough for the whole poll.
        self._refreshLogicalSet(inputSet)

        lastId = getattr(self, watermarkAttr, 0)
        newIds, lastId = self._discoverIdsAfter(inputSet, lastId)
        setattr(self, watermarkAttr, lastId)

        producerClosed = inputSet.isStreamClosed()

        newIds, terminalConsistent = self._reconcileClosedStreamIds(
            inputSet, newIds, knownIds, producerClosed, watermarkAttr)

        items = self._loadLogicalSetItemsByIds(inputSet, newIds,
                                               refresh=False)

        return items, producerClosed, terminalConsistent

    # --------------------------- output state --------------------------
    @staticmethod
    def _getOutputIdSet(outputSet):
        """Published ids, as an id query instead of hydrating every item."""
        if outputSet is None:
            return set()

        getIdSet = getattr(outputSet, 'getIdSet', None)

        if callable(getIdSet):
            return set(getIdSet())

        return {item.getObjId() for item in outputSet}

    @staticmethod
    def _getOutputUniqueValues(outputSet, attribute):
        """One column of the output Set, without hydrating its items.

        Falls back to walking the Set only when the backend cannot answer
        the query, which keeps the fakes used by the tests working.
        """
        if outputSet is None:
            return None

        getUniqueValues = getattr(outputSet, 'getUniqueValues', None)

        if not callable(getUniqueValues):
            return None

        try:
            values = getUniqueValues(attribute)
        except (NotImplementedError, TypeError, KeyError, AttributeError):
            return None

        if values is None:
            return None

        if isinstance(values, dict):
            values = values.get(attribute)

            if values is None:
                return None

        return {value for value in values if value is not None}

    # ------------------------- step graph scans ------------------------
    def _iterKnownStreamingSteps(self):
        """Iterate the steps of this run and of the ones restored on resume.

        Steps restored from a previous execution live in ``_prevSteps``, so
        reading only ``_steps`` would make work that already finished before
        a Continue invisible and schedule it again.
        """
        seen = set()

        for attrName in ('_steps', '_prevSteps'):
            for step in getattr(self, attrName, None) or []:
                stepKey = id(step)

                if stepKey in seen:
                    continue

                seen.add(stepKey)
                yield step

    def _iterUnscannedSteps(self, scan):
        """Yield the steps a previous scan has not consumed yet.

        Steps are appended, never reordered, so each list is walked from
        where the last scan left off. Steps that matched but had not
        finished are carried over explicitly, since those are the only
        already-seen ones whose answer can still change.
        """
        for attrName in ('_steps', '_prevSteps'):
            steps = getattr(self, attrName, None) or []
            start = scan.cursors.get(attrName, 0)

            if start > len(steps):
                # The list was replaced by a shorter one: rescan it.
                start = 0

            for index in range(start, len(steps)):
                yield (attrName, index), steps[index]

            scan.cursors[attrName] = len(steps)

        # Re-read the pending ones from the current lists, so a step that
        # has finished since is seen as finished now.
        for stepKey in list(scan.pending):
            attrName, index = stepKey
            steps = getattr(self, attrName, None) or []

            if index >= len(steps):
                scan.pending.discard(stepKey)
                continue

            yield stepKey, steps[index]

    @staticmethod
    def _getStepFuncName(step):
        funcName = getattr(step, 'funcName', None)

        if hasattr(funcName, 'get'):
            funcName = funcName.get()

        return funcName

    @staticmethod
    def _parseStepArgKeys(step, dictField, keyType):
        """Read the per-item keys a single persisted step carries.

        ``args[0]`` is either the key itself, a list of keys (batch step
        variant), or a serialized object - all three shapes are handled.
        """
        argsStr = getattr(step, 'argsStr', None)

        if hasattr(argsStr, 'get'):
            argsStr = argsStr.get('[]')

        try:
            args = json.loads(argsStr or '[]')
        except (TypeError, ValueError):
            return set()

        if not args:
            return set()

        value = args[0]

        if dictField is not None:
            if isinstance(value, dict):
                itemKey = value.get(dictField)

                if itemKey is not None:
                    return {itemKey}

            return set()

        candidates = value if isinstance(value, list) else [value]

        return {candidate for candidate in candidates
                if keyType is None or isinstance(candidate, keyType)}

    def _getStepArgScan(self, funcNames, onlyFinished, dictField, keyType):
        scans = getattr(self, '_streamingStepScans', None)

        if scans is None:
            scans = {}
            self._streamingStepScans = scans

        scanKey = (tuple(sorted(funcNames)), bool(onlyFinished), dictField,
                   keyType)

        scan = scans.get(scanKey)

        if scan is None:
            scan = _StepArgScan()
            scans[scanKey] = scan

        return scan

    def _collectStepArgKeys(self, funcNames, onlyFinished=True,
                            dictField=None, keyType=None):
        """Collect the per-item keys carried by matching persisted steps.

        :param funcNames: step function names to match.
        :param onlyFinished: when True, only steps already FINISHED count;
            when False, every persisted step counts (i.e. "scheduled").
        :param dictField: set it when the step serializes the item as a dict
            (``args[0][dictField]`` holds the key) instead of carrying the
            key directly.
        :param keyType: optional type filter for the collected keys.

        The scan is incremental: the step graph grows with the stream, so
        re-parsing every step on every poll would cost O(steps processed
        so far) and get slower the longer the protocol runs.
        """
        scan = self._getStepArgScan(funcNames, onlyFinished, dictField,
                                    keyType)

        for stepKey, step in self._iterUnscannedSteps(scan):
            if self._getStepFuncName(step) not in funcNames:
                scan.pending.discard(stepKey)
                continue

            if onlyFinished and not step.isFinished():
                # Keep looking at it: it can still finish later.
                scan.pending.add(stepKey)
                continue

            scan.pending.discard(stepKey)
            scan.keys.update(self._parseStepArgKeys(step, dictField, keyType))

        return set(scan.keys)
