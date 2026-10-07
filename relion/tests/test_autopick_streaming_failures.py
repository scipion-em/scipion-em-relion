import os
from unittest import TestCase
from unittest.mock import patch

from pwem.protocols import ProtParticlePickingAuto
from .logical_set_fakes import LogicalSetFake
from relion.protocols.protocol_autopick_ref import ProtRelion2Autopick
from relion.protocols.protocol_autopick_log import ProtRelionAutopickLoG


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
    def __init__(self):
        self.warnings = []
        self.errors = []

    def _createTmpMicsDir(self, micList):
        return "/work/relion-autopick-batch"

    def _getTmpPath(self, *paths):
        base = "/tmp/relion-autopick-test"
        return os.path.join(base, *paths) if paths else base

    def _pickMicrographsFromStar(self, *args, **kwargs):
        pass

    def warning(self, message):
        self.warnings.append(message)

    def error(self, message):
        self.errors.append(message)


class TestRelionAutopickBatchOutputs(TestCase):
    def test_BatchContinuesWhenOneMicOutputIsMissing(self):
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
            module + ".os.path.exists",
            side_effect=outputExists,
        ):
            protocol._pickMicrographList(micList)

        self.assertTrue(
            any(
                "micrograph 3" in message
                for message in protocol.warnings
            ),
            "A missing autopick output should be reported without aborting "
            "the remaining streaming work.",
        )


class _CrashingAutopickBatchHarness(ProtRelion2Autopick):
    # Regression test harness: relion_autopick (invoked through
    # _pickMicrographsFromStar/runJob) can fail for the whole batch, not
    # just produce a missing output file for one mic. Before this fix,
    # that exception propagated uncaught out of _pickMicrographList and
    # crashed the entire streaming run instead of being reported like
    # any other batch failure.
    def __init__(self):
        self.warnings = []
        self.errors = []

    def _createTmpMicsDir(self, micList):
        return "/work/relion-autopick-batch"

    def _getTmpPath(self, *paths):
        base = "/tmp/relion-autopick-test"
        return os.path.join(base, *paths) if paths else base

    def _pickMicrographsFromStar(self, *args, **kwargs):
        raise RuntimeError("relion_autopick crashed for this batch")

    def warning(self, message):
        self.warnings.append(message)

    def error(self, message):
        self.errors.append(message)


class TestRelionAutopickBatchCrashRegression(TestCase):
    def test_BatchCrashIsReportedInsteadOfAbortingTheRun(self):
        protocol = _CrashingAutopickBatchHarness()
        micList = [_BatchMic(1), _BatchMic(2)]

        module = "relion.protocols.protocol_autopick"
        with patch(
            module + ".convert.createWriter",
            return_value=_BatchWriter(),
        ), patch(
            module + ".os.path.exists",
            return_value=False,
        ):
            # Must not raise.
            protocol._pickMicrographList(micList)

        self.assertEqual(1, len(protocol.errors))
        self.assertIn("batch", protocol.errors[0])
        # The per-mic file-check loop still runs after the caught crash,
        # reporting each mic as missing output rather than silently
        # losing them.
        self.assertEqual(2, len(protocol.warnings))


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

class _NoStorageMicSet:
    def getFileName(self):
        raise AssertionError(
            "Autopick streaming must not depend on a storage filename."
        )


class _AutopickLogicalInputHarness(ProtRelion2Autopick):
    def __init__(self):
        self._mics = _NoStorageMicSet()
        self.micDict = {}
        self.streamClosed = False
        self.loadCalls = 0

    def getInputMicrographs(self):
        return self._mics

    def _loadInputList(self):
        self.loadCalls += 1
        return {}, False

    def _getFirstJoinStep(self):
        return None

    def debug(self, *args, **kwargs):
        pass

    def updateSteps(self):
        raise AssertionError("No new micrographs should have been scheduled.")


class TestRelionAutopickBackendIndependentInputCheck(TestCase):
    def test_CheckNewInputDoesNotDependOnStorageMtime(self):
        protocol = _AutopickLogicalInputHarness()

        protocol._checkNewInput()

        self.assertEqual(
            1,
            protocol.loadCalls,
            "Autopick streaming must refresh logical input state directly "
            "instead of gating discovery on file modification times.",
        )

class _LogicalAutopickMic:
    def __init__(self, objId, name):
        self._objId = objId
        self._name = name

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return self._name

    def clone(self):
        return _LogicalAutopickMic(self._objId, self._name)


class _LogicalAutopickMicSet(LogicalSetFake):
    def __init__(self, items, closed=False):
        super().__init__(items, streamClosed=closed)


class _AutopickLogicalLoadHarness(ProtRelion2Autopick):
    def __init__(self):
        self.micDict = {}
        self.ctfRelations = _Pointer(None)
        self._mics = _LogicalAutopickMicSet(
            [_LogicalAutopickMic(7, "mic_007")],
            closed=True,
        )

    def getInputMicrographs(self):
        return self._mics

    def debug(self, *args, **kwargs):
        pass


class TestRelionAutopickLogicalSetLoading(TestCase):
    def test_LoadInputListUsesLogicalSets(self):
        protocol = _AutopickLogicalLoadHarness()

        newMics, closed = protocol._loadInputList()

        self.assertEqual(["mic_007"], list(newMics.keys()))
        self.assertEqual(1, protocol._mics.loadCalls)
        self.assertTrue(closed)

class _LoGNoStorageMicSet:
    def getFileName(self):
        raise AssertionError(
            "LoG streaming must not depend on a storage filename."
        )


class _LoGLogicalInputHarness(ProtRelionAutopickLoG):
    def __init__(self):
        self._mics = _LoGNoStorageMicSet()
        self.micDict = {}
        self.streamClosed = False
        self.loadCalls = 0

    def getInputMicrographs(self):
        return self._mics

    def _loadInputList(self):
        self.loadCalls += 1
        return {}, False

    def _getFirstJoinStep(self):
        return None

    def debug(self, *args, **kwargs):
        pass

    def updateSteps(self):
        raise AssertionError("No new micrographs should have been scheduled.")


class TestRelionAutopickLoGBackendIndependentInputCheck(TestCase):
    def test_LoGCheckNewInputDoesNotDependOnStorageMtime(self):
        protocol = _LoGLogicalInputHarness()

        protocol._checkNewInput()

        self.assertEqual(
            1,
            protocol.loadCalls,
            "LoG streaming must refresh logical input state directly instead "
            "of gating discovery on file modification times.",
        )

class _LoGLogicalMic:
    def __init__(self, objId, name):
        self._objId = objId
        self._name = name

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return self._name

    def clone(self):
        return _LoGLogicalMic(self._objId, self._name)


class _LoGLogicalMicSet(LogicalSetFake):
    def __init__(self, items, closed=False):
        super().__init__(items, streamClosed=closed)


class _LoGLogicalLoadHarness(ProtRelionAutopickLoG):
    def __init__(self):
        self.micDict = {}
        self._mics = _LoGLogicalMicSet(
            [_LoGLogicalMic(9, "mic_009")],
            closed=True,
        )

    def getInputMicrographs(self):
        return self._mics

    def debug(self, *args, **kwargs):
        pass


class TestRelionAutopickLoGLogicalSetLoading(TestCase):
    def test_LoGLoadInputListUsesLogicalSets(self):
        protocol = _LoGLogicalLoadHarness()

        newMics, closed = protocol._loadInputList()

        self.assertEqual(["mic_009"], list(newMics.keys()))
        self.assertEqual(1, protocol._mics.loadCalls)
        self.assertTrue(closed)
