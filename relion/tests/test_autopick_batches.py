# **************************************************************************
# *
# * With a streaming batch size above one, the micrographs that do not fill
# * a batch yet wait for the next poll instead of being dropped.
# *
# **************************************************************************

import unittest
from collections import OrderedDict

from relion.protocols.protocol_autopick import ProtRelionAutopickBase


class _Mic:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return "mic_%03d" % self._objId


class _Batches(ProtRelionAutopickBase):
    def __init__(self, batchSize):
        self.micDict = OrderedDict()
        self.initialIds = []
        self.streamClosed = False
        self._batchSize = batchSize
        self.batches = []

    def _getStreamingBatchSize(self):
        return self._batchSize

    def _insertPickMicrographListStep(self, micList, prerequisites, *args):
        self.batches.append([mic.getObjId() for mic in micList])
        return len(self.batches)

    def _insertPickMicrographStep(self, mic, prerequisites, *args):
        return self._insertPickMicrographListStep([mic], prerequisites)

    def _getPickArgs(self):
        return []


class TestBatchRemainders(unittest.TestCase):
    def test_TheRemainderOfABatchWaitsForTheNextPoll(self):
        protocol = _Batches(batchSize=3)

        protocol._insertNewMicsSteps([_Mic(i) for i in (1, 2, 3, 4, 5)])
        self.assertEqual([[1, 2, 3]], protocol.batches)

        # Discovery offers only what is new: 4 and 5 must not be lost.
        protocol._insertNewMicsSteps([_Mic(6)])
        self.assertEqual([[1, 2, 3], [4, 5, 6]], protocol.batches)

    def test_ClosingTheStreamTakesTheLastIncompleteBatch(self):
        protocol = _Batches(batchSize=3)
        protocol._insertNewMicsSteps([_Mic(i) for i in (1, 2, 3, 4)])

        protocol.streamClosed = True
        protocol._insertNewMicsSteps([])

        self.assertEqual([[1, 2, 3], [4]], protocol.batches)
        self.assertFalse(protocol._hasBatchRemainder())

    def test_AMicrographIsNeverScheduledTwice(self):
        protocol = _Batches(batchSize=2)
        mics = [_Mic(1), _Mic(2)]

        protocol._insertNewMicsSteps(mics)
        protocol._insertNewMicsSteps(mics)

        self.assertEqual([[1, 2]], protocol.batches)


if __name__ == '__main__':
    unittest.main()
