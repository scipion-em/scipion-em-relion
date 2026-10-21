# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import threading

import pyworkflow.object as pwobj
from pyworkflow.protocol.constants import STATUS_NEW
import pyworkflow.utils as pwutils
from pwem.protocols import ProtParticlePickingAuto

import relion.convert as convert
from .protocol_base import ProtRelionBase
from .protocol_streaming_base import RelionStreamingBase

# Module level: these helpers are called unbound on light test harnesses.
PICKING_STEP_NAMES = ('pickMicrographStep', 'pickMicrographListStep')
# The micrographs whose picking failed, as a durable list in extra/: a
# micrograph with no coordinates is otherwise indistinguishable from one
# without particles, and Continue picks the listed ones again.
FAILED_PICKING_LIST = 'FAILED_all.TXT'
_failedListLock = threading.Lock()


class ProtRelionAutopickBase(RelionStreamingBase, ProtParticlePickingAuto,
                             ProtRelionBase):
    """ Base class for auto-picking protocols in Relion.
    """
    _label = None

    def _loadSet(self, inputSet, SetClass, getKeyFunc, watermarkAttr=None,
                 knownIds=None):
        """Discover the items added to one input stream since the last poll.

        pwem rebuilds the Set from its storage filename and walks it whole.
        This queries the ids above that stream's watermark and loads only
        those, so a poll costs what just arrived rather than everything
        the stream has produced so far.
        """
        watermarkAttr = watermarkAttr or '_lastInputId'

        if knownIds is None:
            knownIds = self._getKnownStreamIds(watermarkAttr)

        newItems, producerClosed, terminalConsistent = (
            self._discoverNewInputItems(inputSet, watermarkAttr, knownIds))

        newItemDict = {}

        for item in newItems:
            itemId = item.getObjId()

            if itemId in knownIds:
                continue

            knownIds.add(itemId)
            itemKey = getKeyFunc(item)

            if itemKey not in self.micDict:
                newItemDict[itemKey] = item

        return newItemDict, producerClosed and terminalConsistent

    def _loadMics(self, micSet):
        """Micrographs discovered on their own watermark."""
        return self._loadSet(micSet, None, lambda mic: mic.getMicName(),
                             watermarkAttr='_lastMicId')

    def _loadCTFs(self, ctfSet):
        """CTFs discovered on a separate watermark: they are their own
        stream and advance independently of the micrographs."""
        return self._loadSet(ctfSet, None,
                             lambda ctf: ctf.getMicrograph().getMicName(),
                             watermarkAttr='_lastCtfId')

    def _insertNewMicsSteps(self, inputMics):
        return self._insertKeepingBatchRemainder(
            inputMics,
            lambda mics: ProtParticlePickingAuto._insertNewMicsSteps(self, mics))

    def _hasActiveStreamingWork(self):
        """Micrographs waiting for a batch, or picked and not published."""
        return self._hasBatchRemainder() or not getattr(
            self, '_allPicksPublished', False)

    def _checkNewInput(self):
        # Refresh logical input state directly. Do not gate discovery on
        # storage filenames or filesystem modification times.
        micDict, self.streamClosed = self._loadInputList()
        outputStep = self._getFirstJoinStep()

        if micDict or self._hasBatchRemainder():
            deps = self._insertNewMicsSteps(micDict.values())
            if outputStep is not None:
                outputStep.addPrerequisites(*deps)
            self.updateSteps()

    def _getFinishedPickingMicNames(self):
        """Micrograph names carried by finished picking steps."""
        return self._collectStepArgKeys(PICKING_STEP_NAMES, keyType=str)

    def _getPublishedPickingMicIds(self):
        """Micrograph ids already represented in the output coordinates.

        An id query on the output, not a walk over it.
        """
        micIds = self._getOutputUniqueValues(
            getattr(self, 'outputCoordinates', None), '_micId')

        return set() if micIds is None else micIds

    def _checkNewOutput(self):
        """Publish finished picking without DONE sidecars.

        pwem decides this with extra/DONE/mic_*.TXT plus DONE/all.TXT.
        The persisted step graph already records what finished and the
        output Set what was published, so neither file is consulted and a
        cleaned extra/ cannot make finished work look pending.
        """
        if getattr(self, 'finished', False):
            return

        publishedMicIds = self._getPublishedPickingMicIds()
        finishedMicNames = self._getFinishedPickingMicNames()
        # Read in this run already: a micrograph without coordinates (no
        # particles, or a failed picking) is never in the output, and would
        # otherwise be read again every poll.
        handledMicIds = self.__dict__.setdefault('_handledPickingMicIds',
                                                 set())

        listOfMics = list(self.micDict.values())
        newDone = [mic for mic in listOfMics
                   if mic.getMicName() in finishedMicNames
                   and mic.getObjId() not in publishedMicIds
                   and mic.getObjId() not in handledMicIds]

        doneCount = len([mic for mic in listOfMics
                         if mic.getObjId() in publishedMicIds
                         or mic.getObjId() in handledMicIds])
        allDone = doneCount + len(newDone)

        self.debug('_checkNewOutput: mics=%s published=%s newDone=%s'
                   % (len(listOfMics), doneCount, len(newDone)))

        self._allPicksPublished = allDone == len(listOfMics)
        self.finished = (self.streamClosed and self._allPicksPublished
                         and not self._hasBatchRemainder())
        streamMode = (pwobj.Set.STREAM_CLOSED
                      if self.finished else pwobj.Set.STREAM_OPEN)

        if newDone:
            self._updateOutputCoordSet(newDone, streamMode)
            handledMicIds.update(mic.getObjId() for mic in newDone)
        elif not self.finished:
            self._streamingSleepOnWait()
            return

        if self.finished:
            self._updateStreamState(streamMode)
            self._failIfEveryMicrographFailed()

            # The lifecycle here is still the classic one: createOutputStep
            # was scheduled with wait=True and only runs once something
            # releases it. Dropping this leaves it WAITING for good and the
            # protocol never finishes, however complete the output is.
            outputStep = self._getFirstJoinStep()

            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(STATUS_NEW)

    def _failIfEveryMicrographFailed(self):
        """A failed micrograph is left out; the run fails only when that
        happened to every micrograph."""
        handled = getattr(self, '_handledPickingMicIds', set())
        output = getattr(self, 'outputCoordinates', None)

        if (handled and handled <= self._failedPickingMicIds()
                and (output is None or output.getSize() == 0)):
            raise Exception("Picking failed for every one of the %d "
                            "micrographs; see the errors above. Continue "
                            "picks them again." % len(handled))

    # ------------------------- failed micrographs ------------------------
    def _failedPickingMicIds(self):
        """Ids of the micrographs whose picking failed, from extra/."""
        cached = getattr(self, '_failedPickingMicIdsCache', None)

        if cached is None:
            cached = set()
            path = self._getExtraPath(FAILED_PICKING_LIST)

            if os.path.exists(path):
                with open(path) as handle:
                    cached = {int(line) for line in handle if line.strip()}

            self._failedPickingMicIdsCache = cached

        return set(cached)

    def _writeFailedPicking(self, micIds):
        with open(self._getExtraPath(FAILED_PICKING_LIST), 'w') as handle:
            handle.writelines('%d\n' % micId for micId in sorted(micIds))

        self._failedPickingMicIdsCache = set(micIds)

    def _recordFailedPicking(self, mic):
        with _failedListLock:
            failed = self._failedPickingMicIds()

            if mic.getObjId() not in failed:
                self._writeFailedPicking(failed | {mic.getObjId()})

    def _forgetFailedPicking(self, micIds):
        """Picked well now: no longer listed as failed."""
        with _failedListLock:
            failed = self._failedPickingMicIds()

            if failed & set(micIds):
                self._writeFailedPicking(failed - set(micIds))

    def _forgetDoneMarkersOfFailed(self, micNames):
        """pwem marks a micrograph DONE whatever its picking gave, and a
        Continue skips the ones marked: a failed one must be picked again."""
        failed = self._failedPickingMicIds()

        for micName in micNames:
            mic = self.micDict.get(micName)

            if mic is not None and mic.getObjId() in failed:
                pwutils.cleanPath(self._getMicDone(mic))

    def pickMicrographStep(self, micName, *args):
        ProtParticlePickingAuto.pickMicrographStep(self, micName, *args)
        self._forgetDoneMarkersOfFailed([micName])

    def pickMicrographListStep(self, micNameList, *args):
        ProtParticlePickingAuto.pickMicrographListStep(self, micNameList,
                                                       *args)
        self._forgetDoneMarkersOfFailed(micNameList)

    def _pickMicrograph(self, mic, *args):
        """ This method should be invoked only when working in streaming mode.
        """
        self._pickMicrographList([mic], *args)

    def _pickMicrographList(self, micList, *args):
        if not micList:
            return

        micsDir = self._createTmpMicsDir(micList)
        micStar = os.path.join(micsDir, 'input_micrographs.star')
        writer = convert.createWriter(rootDir=micsDir, outputDir=micsDir)
        writer.writeSetOfMicrographs(micList, micStar)
        try:
            # pickMicrographListStep (pwem) has no exception boundary of
            # its own around this call - a single relion_autopick crash
            # for the whole batch would otherwise propagate uncaught and
            # abort the entire streaming run. Its micrographs are recorded
            # as failed instead: left out of the output, picked again on
            # Continue, and the run fails only if every one failed.
            self._pickMicrographsFromStar(micStar, micsDir, *args)
        except Exception as e:
            self.error(
                "ERROR: Autopick failed for micrograph batch starting at "
                "%s with the exception %s" % (micList[0].getObjId(), e)
            )

        # Move each expected coordinates file explicitly. Do not let a
        # successful wildcard move hide a missing output for one mic in a batch.
        pickedIds = []

        for mic in micList:
            fileName = 'mic_%06d_autopick.star' % mic.getObjId()
            srcFile = os.path.join(micsDir, fileName)
            dstFile = self._getTmpPath(fileName)

            if os.path.exists(srcFile):
                pwutils.moveFile(srcFile, dstFile)
                pickedIds.append(mic.getObjId())
            elif os.path.exists(dstFile):
                pickedIds.append(mic.getObjId())
            else:
                self.error(
                    "ERROR: Missing autopick output for micrograph %s: %s / %s"
                    % (mic.getObjId(), srcFile, dstFile)
                )
                self._recordFailedPicking(mic)

        self._forgetFailedPicking(pickedIds)

    def _createSetOfCoordinates(self, micSet, suffix=''):
        """ Override this method to set the box size. """
        coordSet = ProtParticlePickingAuto._createSetOfCoordinates(
            self, micSet, suffix=suffix)
        coordSet.setBoxSize(self.getBoxSize())
        return coordSet

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        """ Parse back the output star files and populate the SetOfCoordinates.
        """
        template = self._getTmpPath("mic_%06d_autopick.star")
        failedIds = self._failedPickingMicIds()

        for mic in micList:
            if mic.getObjId() in failedIds:
                continue

            # A micrograph is published whole or not at all: its rows are
            # kept apart until its whole file has been read.
            micCoords = _Collected()
            starFile = template % mic.getObjId()

            try:
                if not os.path.exists(starFile):
                    raise FileNotFoundError(starFile)

                convert.readCoordinates(mic, starFile, micCoords)
            except Exception as e:
                self.error("ERROR: Could not read the coordinates of "
                           "micrograph %s: %s" % (mic.getObjId(), e))
                self._recordFailedPicking(mic)
                continue

            for coord in micCoords.items:
                coordSet.append(coord)

    def _pickMicrographsFromStar(self, micStar, micsDir, *args):
        """ Should be defined in subclasses. """
        pass

    def getBoxSize(self):
        """ Return a reasonable box-size in pixels. """
        return None

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    def __getMicListPrefix(self, micList):
        n = len(micList)
        if n == 0:
            raise ValueError("Empty micrographs list!")
        micsPrefix = 'mic_%06d' % micList[0].getObjId()
        if n > 1:
            micsPrefix += "-%06d" % micList[-1].getObjId()
        return micsPrefix

    def _createTmpMicsDir(self, micList):
        """ Create a temporary path to work with a list of micrographs. """
        micsDir = self._getTmpPath(self.__getMicListPrefix(micList))
        pwutils.makePath(micsDir)
        return micsDir


class _Collected:
    """Holds what a reader appends, to be published only if it all reads."""

    def __init__(self):
        self.items = []

    def append(self, item):
        self.items.append(item)
