# **************************************************************************
# *
# * Cost of one streaming poll in the Relion SPA protocols.
# *
# **************************************************************************

"""A poll must cost what just arrived, not everything seen so far.

These protocols run for as long as their producer does and their input
Sets keep growing, so anything done per poll over every item seen makes
them slower the longer they run - exactly when there is most data.
"""

import unittest

from relion.protocols.protocol_autopick import ProtRelionAutopickBase
from relion.protocols.protocol_streaming_base import RelionStreamingBase

from .logical_set_fakes import LogicalSetFake


class _Movie:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return "mic_%03d" % self._objId

    def clone(self):
        return _Movie(self._objId)


class _Harness(RelionStreamingBase):
    def __init__(self):
        self._lastInputId = 0

    def debug(self, *args, **kwargs):
        pass


class _StaleCloseSet(LogicalSetFake):
    """Its closed flag only becomes visible after a reload."""

    def __init__(self, items, closedAfterReloads=2):
        super().__init__(items, streamClosed=False)
        self.closedAfterReloads = closedAfterReloads

    def loadAllProperties(self):
        super().loadAllProperties()

        if self.reloads >= self.closedAfterReloads:
            self._streamClosed = True


def _items(firstId, count):
    return [_Movie(objId) for objId in range(firstId, firstId + count)]


class TestRelionStreamingPollCost(unittest.TestCase):
    def testPollOnlyHydratesWhatJustArrived(self):
        inputSet = LogicalSetFake(_items(1, 500))
        protocol = _Harness()
        known = set()

        items, _, _ = protocol._discoverNewInputItems(
            inputSet, '_lastInputId', known)
        known.update(item.getObjId() for item in items)

        self.assertEqual(500, len(items))
        self.assertEqual(500, inputSet.hydratedItems)

        inputSet.addItems(_items(501, 3))
        hydratedBefore = inputSet.hydratedItems

        items, _, _ = protocol._discoverNewInputItems(
            inputSet, '_lastInputId', known)

        self.assertEqual(3, len(items))
        self.assertEqual(3, inputSet.hydratedItems - hydratedBefore)
        self.assertEqual(0, inputSet.fullScans)

    def testIdlePollHydratesNothingAtAll(self):
        inputSet = LogicalSetFake(_items(1, 500))
        protocol = _Harness()
        known = set()

        items, _, _ = protocol._discoverNewInputItems(
            inputSet, '_lastInputId', known)
        known.update(item.getObjId() for item in items)
        hydratedBefore = inputSet.hydratedItems

        protocol._discoverNewInputItems(inputSet, '_lastInputId', known)

        self.assertEqual(hydratedBefore, inputSet.hydratedItems)

    def testStreamClosingIsSeenEvenWhenNothingNewArrives(self):
        # isStreamClosed() reads a Set property, so a poll that skipped the
        # reload would never notice the producer closing and would spin
        # forever.
        inputSet = _StaleCloseSet(_items(1, 2))
        protocol = _Harness()
        known = set()

        _, producerClosed, _ = protocol._discoverNewInputItems(
            inputSet, '_lastInputId', known)
        self.assertFalse(producerClosed)

        known = {1, 2}
        _, producerClosed, _ = protocol._discoverNewInputItems(
            inputSet, '_lastInputId', known)
        self.assertTrue(producerClosed)

    def testEachInputStreamKeepsItsOwnWatermark(self):
        # A protocol watching coordinates, micrographs and CTFs must not
        # let one stream's progress hide another's new items.
        protocol = _Harness()

        self.assertIsNot(
            protocol._getKnownStreamIds('_lastMicId'),
            protocol._getKnownStreamIds('_lastCtfId'),
        )

        mics = LogicalSetFake(_items(1, 5), streamClosed=True)
        ctfs = LogicalSetFake(_items(1, 5), streamClosed=True)

        micItems, _, _ = protocol._discoverNewInputItems(
            mics, '_lastMicId', protocol._getKnownStreamIds('_lastMicId'))
        ctfItems, _, _ = protocol._discoverNewInputItems(
            ctfs, '_lastCtfId', protocol._getKnownStreamIds('_lastCtfId'))

        self.assertEqual(5, len(micItems))
        self.assertEqual(5, len(ctfItems))

    def testResumeStartsAboveWhatIsDoneButKeepsTheGaps(self):
        inputSet = LogicalSetFake(_items(1, 100))
        protocol = _Harness()

        # Everything was processed except item 42.
        processedIds = set(range(1, 101)) - {42}

        watermark, gapIds = protocol._resumeWatermarkWithGaps(
            inputSet, processedIds)

        self.assertEqual(100, watermark)
        self.assertEqual({42}, gapIds)
        self.assertEqual(0, inputSet.hydratedItems)

    def testStepGraphScanDoesNotReparseStepsItAlreadyRead(self):
        protocol = _Harness()

        parsed = []
        original = RelionStreamingBase._parseStepArgKeys

        def countingParse(step, dictField, keyType):
            parsed.append(step)
            return original(step, dictField, keyType)

        # A plain function, not staticmethod(): an instance attribute is
        # never bound, and a staticmethod object is only callable itself
        # from Python 3.10 on.
        protocol._parseStepArgKeys = countingParse
        protocol._steps = [_FinishedStep('["mic_%03d"]' % i)
                           for i in range(1, 201)]

        self.assertEqual(200, len(protocol._collectStepArgKeys(('pickStep',),
                                                               keyType=str)))
        self.assertEqual(200, len(parsed))

        protocol._steps.append(_FinishedStep('["mic_201"]'))

        self.assertEqual(201, len(protocol._collectStepArgKeys(('pickStep',),
                                                               keyType=str)))
        self.assertEqual(201, len(parsed))

    def testStepGraphScanReadsStepsRestoredOnResume(self):
        # Work finished before a Continue lives in _prevSteps; reading only
        # _steps would make it invisible and schedule it all over again.
        protocol = _Harness()
        protocol._steps = [_FinishedStep('["mic_002"]')]
        protocol._prevSteps = [_FinishedStep('["mic_001"]')]

        self.assertEqual(
            {"mic_001", "mic_002"},
            protocol._collectStepArgKeys(('pickStep',), keyType=str),
        )


class _FinishedStep:
    funcName = 'pickStep'

    def __init__(self, argsStr):
        self.argsStr = argsStr

    def isFinished(self):
        return True


if __name__ == "__main__":
    unittest.main()


class _Mic:
    def __init__(self, objId, micName):
        self._objId = objId
        self._micName = micName

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return self._micName


class _CoordOutput(LogicalSetFake):
    def __init__(self, micIds):
        super().__init__([], streamClosed=False)
        self._micIds = list(micIds)

    def getUniqueValues(self, attributes, where=None):
        if attributes == '_micId':
            return list(self._micIds)

        return super().getUniqueValues(attributes, where=where)


class _NoSidecarPickingHarness(ProtRelionAutopickBase):
    """Fails loudly if completion is read from or written to a DONE file."""

    def __init__(self, steps, publishedMicIds):
        self.micDict = {mic.getMicName(): mic
                        for mic in (_Mic(1, "mic_001"), _Mic(2, "mic_002"))}
        self._steps = steps
        self.streamClosed = False
        self.published = []
        self.outputCoordinates = _CoordOutput(publishedMicIds)

    def _readDoneList(self):
        raise AssertionError("Completion must not be read from DONE/all.TXT.")

    def _writeDoneList(self, micList):
        raise AssertionError("Completion must not be written to DONE/all.TXT.")

    def _isMicDone(self, mic):
        raise AssertionError(
            "Completion must not depend on a per-micrograph DONE file.")

    def _updateOutputCoordSet(self, micList, streamMode):
        self.published.append([mic.getMicName() for mic in micList])
        return list(micList)

    def _streamingSleepOnWait(self):
        pass

    def debug(self, *args, **kwargs):
        pass


class TestRelionCompletionWithoutSidecars(unittest.TestCase):
    def testPickingCompletionComesFromTheStepGraph(self):
        protocol = _NoSidecarPickingHarness(
            steps=[_FinishedPickStep('["mic_001", {}]')],
            publishedMicIds=[],
        )

        protocol._checkNewOutput()

        self.assertEqual([["mic_001"]], protocol.published)

    def testAlreadyPublishedMicrographIsNotPublishedTwice(self):
        protocol = _NoSidecarPickingHarness(
            steps=[_FinishedPickStep('["mic_001", {}]')],
            publishedMicIds=[1],
        )

        protocol._checkNewOutput()

        self.assertEqual([], protocol.published)


class _FinishedPickStep:
    funcName = 'pickMicrographStep'

    def __init__(self, argsStr):
        self.argsStr = argsStr

    def isFinished(self):
        return True


class _JoinStep:
    """Stands in for the createOutputStep scheduled with wait=True."""

    def __init__(self):
        self.status = 'waiting'

    def isWaiting(self):
        return self.status == 'waiting'

    def setStatus(self, value):
        self.status = value


class _FinishingPickingHarness(_NoSidecarPickingHarness):
    """Everything picked and published, with the stream closed."""

    def __init__(self):
        super().__init__(
            steps=[_FinishedPickStep('["mic_001", {}]'),
                   _FinishedPickStep('["mic_002", {}]')],
            publishedMicIds=[],
        )
        self.streamClosed = True
        self.joinStep = _JoinStep()
        self.streamStates = []

    def _getFirstJoinStep(self):
        return self.joinStep

    def _updateStreamState(self, streamMode):
        self.streamStates.append(streamMode)


class TestRelionReleasesTheWaitingOutputStep(unittest.TestCase):
    """These protocols keep the classic lifecycle, so the output step is
    scheduled with wait=True and stays WAITING until _checkNewOutput puts
    it back to NEW. Forget that and the protocol hangs with its output
    complete - which is exactly what happened to a real autopick run."""

    def testFinishingReleasesTheWaitingOutputStep(self):
        protocol = _FinishingPickingHarness()

        protocol._checkNewOutput()

        self.assertTrue(protocol.finished)
        self.assertEqual(
            'new',
            protocol.joinStep.status,
            "A finished stream must release createOutputStep, or the "
            "protocol waits for ever with nothing left to do.",
        )

    def testOutputStepIsNotReleasedWhileTheStreamIsStillOpen(self):
        protocol = _FinishingPickingHarness()
        protocol.streamClosed = False

        protocol._checkNewOutput()

        self.assertFalse(protocol.finished)
        self.assertEqual('waiting', protocol.joinStep.status)
