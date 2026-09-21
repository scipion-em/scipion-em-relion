import os
from unittest import TestCase
from unittest.mock import patch

from pwem.protocols import ProtParticlePickingAuto
from relion.protocols.protocol_autopick_ref import ProtRelion2Autopick


class _ClosedMicrographs:
    def isStreamOpen(self):
        return False

    def strId(self):
        return "mics-1"

    def __iter__(self):
        return iter(())


class _OpenCtfSet:
    def isStreamOpen(self):
        return True


class _Pointer:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _References:
    def strId(self):
        return "refs-1"


class _AutopickHarness(ProtRelion2Autopick):
    def __init__(self):
        self.streamingBatchSize = 0
        self.ctfRelations = _Pointer(_OpenCtfSet())
        self._mics = _ClosedMicrographs()

    def getInputMicrographs(self):
        return self._mics

    def getInputReferences(self):
        return _References()

    def _loadInputList(self):
        return {}, False

    def _insertFunctionStep(self, *args, **kwargs):
        return 1

    def _getPickArgs(self):
        return []

    def usesGpu(self):
        return False

    def createOutputStep(self):
        pass


class TestRelionAutopickStreamingFailures(TestCase):
    def test_BatchZeroKeepsStreamingWhenCtfStreamIsOpen(self):
        protocol = _AutopickHarness()

        with patch.object(
            ProtParticlePickingAuto,
            "_insertAllSteps",
        ) as baseInsert:
            protocol._insertAllSteps()

        baseInsert.assert_called_once_with(protocol)
class _ClosedCtfSet:
    def isStreamOpen(self):
        return False


class _PersistedStreamingStep:
    funcName = "_doNothing"

    def isWaiting(self):
        return True


class _AutopickResumeHarness(_AutopickHarness):
    def __init__(self):
        super().__init__()
        self.ctfRelations = _Pointer(_ClosedCtfSet())

    def isContinued(self):
        return True

    def loadSteps(self):
        return [_PersistedStreamingStep()]


class TestRelionAutopickStreamingResume(TestCase):
    def test_ContinueKeepsOriginalStreamingModeAfterInputsClose(self):
        protocol = _AutopickResumeHarness()

        with patch.object(
            ProtParticlePickingAuto,
            "_insertAllSteps",
        ) as baseInsert:
            protocol._insertAllSteps()

        baseInsert.assert_called_once_with(protocol)
class _BatchMic:
    def __init__(self, mic_id):
        self._mic_id = mic_id

    def getObjId(self):
        return self._mic_id


class _BatchWriter:
    def writeSetOfMicrographs(self, *args, **kwargs):
        pass


class _AutopickBatchOutputHarness(ProtRelion2Autopick):
    def _createTmpMicsDir(self, micList):
        return "/work/relion-autopick-batch"

    def _getTmpPath(self, *paths):
        base = "/tmp/relion-autopick-test"
        return os.path.join(base, *paths) if paths else base

    def _pickMicrographsFromStar(self, *args, **kwargs):
        pass


class TestRelionAutopickBatchOutputs(TestCase):
    def test_BatchFailsWhenOneMicOutputIsMissing(self):
        protocol = _AutopickBatchOutputHarness()
        micList = [_BatchMic(1), _BatchMic(2), _BatchMic(3)]

        def outputExists(path):
            return path in {
                "/tmp/relion-autopick-test/mic_000001_autopick.star",
                "/tmp/relion-autopick-test/mic_000002_autopick.star",
            }

        module = "relion.protocols.protocol_autopick"
        with patch(
            module + ".convert.createWriter",
            return_value=_BatchWriter(),
        ), patch(
            module + ".os.system",
            return_value=0,
        ), patch(
            module + ".os.path.exists",
            side_effect=outputExists,
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "Missing autopick output",
            ):
                protocol._pickMicrographList(micList)
class _NamedMic:
    def __init__(self, mic_id):
        self._mic_id = mic_id

    def getMicName(self):
        return "mic_%03d" % self._mic_id


class _ClosedBatchMicrographs:
    def __init__(self, micList):
        self._micList = micList

    def isStreamOpen(self):
        return False

    def strId(self):
        return "mics-batch"

    def __iter__(self):
        return iter(self._micList)


class _AutopickFilteredBatchHarness(ProtRelion2Autopick):
    def __init__(self):
        self.streamingBatchSize = 0
        self.ctfRelations = _Pointer(_ClosedCtfSet())
        self.allMics = [_NamedMic(1), _NamedMic(2), _NamedMic(3)]
        self.readyMics = self.allMics[:2]
        self._mics = _ClosedBatchMicrographs(self.allMics)
        self.pickNames = None

    def isContinued(self):
        return False

    def getInputMicrographs(self):
        return self._mics

    def getInputReferences(self):
        return _References()

    def _loadInputList(self):
        return {
            mic.getMicName(): mic
            for mic in self.readyMics
        }, True

    def _insertFunctionStep(self, func, *args, **kwargs):
        if getattr(func, "__name__", None) == "pickMicrographListStep":
            self.pickNames = args[0]
        return 1

    def _getPickArgs(self):
        return []

    def usesGpu(self):
        return False

    def createOutputStep(self):
        pass


class TestRelionAutopickBatchFiltering(TestCase):
    def test_BatchZeroSchedulesOnlyMicrographsReadyForCtf(self):
        protocol = _AutopickFilteredBatchHarness()

        protocol._insertAllSteps()

        self.assertEqual(
            ["mic_001", "mic_002"],
            protocol.pickNames,
            "The batch-zero path must schedule only micrographs returned by "
            "_loadInputList(), since those are the ones whose required CTF "
            "input is ready.",
        )
