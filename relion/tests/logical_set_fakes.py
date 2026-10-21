# **************************************************************************
# *
# * Test doubles for the logical Set API used by Relion streaming.
# *
# **************************************************************************

"""Fakes that answer id queries the way a real backend does.

Streaming discovery must cost what just arrived, not everything the stream
has produced, so these doubles implement ``getUniqueValues``/``getIdSet``
and the ``where`` clauses that go with them. They also count how many items
were hydrated, which is what lets a test assert that a poll did not walk
the whole Set.
"""

import re


_ID_GT = re.compile(r'^\s*id\s*>\s*(-?\d+)\s*$')
_ID_LE = re.compile(r'^\s*id\s*<=\s*(-?\d+)\s*$')
_ID_IN = re.compile(r'^\s*id\s+IN\s*\(([\d,\s]*)\)\s*$', re.IGNORECASE)


def _matchesWhere(itemId, where):
    """Evaluate the small subset of SQL the streaming helpers emit."""
    if where is None:
        return True

    match = _ID_GT.match(where)

    if match:
        return itemId > int(match.group(1))

    match = _ID_LE.match(where)

    if match:
        return itemId <= int(match.group(1))

    match = _ID_IN.match(where)

    if match:
        wanted = {int(part) for part in match.group(1).split(',')
                  if part.strip()}
        return itemId in wanted

    raise AssertionError('Unsupported where clause: %r' % where)


class LogicalSetFake:
    """A logical Set that answers by id without being walked."""

    # A polling loop that never sees the stream close would hang the
    # suite instead of failing it, so these doubles refuse to be polled
    # forever. See isStreamClosed below.
    MAX_POLLS = 50

    def __init__(self, items, streamClosed=False):
        self._items = list(items)
        self._streamClosed = streamClosed
        self.loadCalls = 0
        self.reloads = 0
        self.hydratedItems = 0
        self.fullScans = 0
        self.closedChecks = 0

    # The storage filename is never part of the streaming contract.
    def getFileName(self):
        raise AssertionError(
            'Streaming discovery must not depend on a storage filename.'
        )

    def loadAllProperties(self):
        self.loadCalls += 1
        self.reloads += 1

    def isStreamClosed(self):
        self.closedChecks += 1

        if self.closedChecks > self.MAX_POLLS:
            raise AssertionError(
                "Polled %d times without ever seeing the stream close. "
                "A protocol that misses the producer closing must fail "
                "this suite, not hang it."
                % self.closedChecks
            )

        return self._streamClosed

    def isStreamOpen(self):
        return not self._streamClosed

    def setStreamClosed(self, streamClosed):
        self._streamClosed = streamClosed

    def getSize(self):
        return len(self._items)

    def addItems(self, items):
        self._items.extend(items)

    def getUniqueValues(self, attributes, where=None):
        if attributes != 'id':
            raise NotImplementedError(attributes)

        return [item.getObjId() for item in self._items
                if _matchesWhere(item.getObjId(), where)]

    def getIdSet(self):
        return set(self.getUniqueValues('id'))

    def iterItems(self, orderBy=None, direction=None, where=None):
        if where is None:
            self.fullScans += 1

        items = [item for item in self._items
                 if _matchesWhere(item.getObjId(), where)]

        if orderBy == 'id':
            items.sort(key=lambda item: item.getObjId(),
                       reverse=direction == 'DESC')

        self.hydratedItems += len(items)

        return iter(items)

    def __iter__(self):
        return self.iterItems()
