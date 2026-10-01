from unittest import TestCase
from unittest.mock import MagicMock, patch

from relion.protocols.protocol_extract_particles import (
    ProtRelionExtractParticles,
)


class _Acquisition:
    def getMagnification(self):
        return 10000

    def getVoltage(self):
        return 300

    def getAmplitudeContrast(self):
        return 0.1

    def getSphericalAberration(self):
        return 2.7


class _InputMicrographs:
    def getAcquisition(self):
        return _Acquisition()


class _Mic:
    def getObjId(self):
        return 7

    def getAttributeValue(self, name, default=None):
        return default

    def getCTF(self):
        return None


class _ProtocolHarness(ProtRelionExtractParticles):
    def __init__(self):
        self.coordDict = {7: []}

    def getInputMicrographs(self):
        return _InputMicrographs()

    def _getTmpPath(self, *paths):
        base = "/tmp/relion-extract-test"
        return os.path.join(base, *paths) if paths else base

    def _getExtraPath(self, *paths):
        base = "/extra/relion-extract-test"
        return os.path.join(base, *paths) if paths else base

    def warning(self, *args, **kwargs):
        pass


class TestRelionExtractParticlesStreamingFailures(TestCase):
    def test_ResumeReusesAlreadyMovedParticleStack(self):
        protocol = _ProtocolHarness()
        mic = _Mic()
        outputParts = MagicMock()

        tmpStack = "/tmp/relion-extract-test/mic_000007.mrcs"
        extraStack = "/extra/relion-extract-test/mic_000007.mrcs"

        def exists(path):
            if path == tmpStack:
                return False
            if path == extraStack:
                return True
            return False

        with patch(
            "relion.protocols.protocol_extract_particles.relion.convert.Table",
            return_value=[],
        ), patch(
            "relion.protocols.protocol_extract_particles.os.path.exists",
            side_effect=exists,
        ), patch(
            "relion.protocols.protocol_extract_particles.pwutils.moveFile"
        ) as moveFile:
            protocol.readPartsFromMics([mic], outputParts)

        moveFile.assert_not_called()


class _MultiMic:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId

    def getAttributeValue(self, name, default=None):
        return default

    def getCTF(self):
        return None


class _PerMicIsolationHarness(ProtRelionExtractParticles):
    def __init__(self):
        self.coordDict = {1: [], 2: []}
        self.errors = []

    def getInputMicrographs(self):
        return _InputMicrographs()

    def _getTmpPath(self, *paths):
        base = "/tmp/relion-extract-isolation"
        return os.path.join(base, *paths) if paths else base

    def _getExtraPath(self, *paths):
        base = "/extra/relion-extract-isolation"
        return os.path.join(base, *paths) if paths else base

    def warning(self, *args, **kwargs):
        pass

    def error(self, message):
        self.errors.append(message)


class TestRelionExtractParticlesPerMicIsolation(TestCase):
    # Regression test: readPartsFromMics used to let one mic's exception
    # (e.g. a missing particle stack) propagate out of the whole batch
    # loop, losing every OTHER mic's work in the same batch too. pwem's
    # own outer try/except (_updateOutputPartSet) only prevents a
    # protocol-wide crash - it still treats the whole pending batch as
    # lost. Per-mic isolation means mic 2 must still be attempted (and
    # have its coordDict entry cleaned up) even though mic 1 fails.
    def test_OneMicFailureDoesNotStopTheRestOfTheBatch(self):
        protocol = _PerMicIsolationHarness()
        mic1 = _MultiMic(1)
        mic2 = _MultiMic(2)
        outputParts = MagicMock()

        def exists(path):
            # mic 1's stack is missing everywhere -> must fail and be
            # reported. mic 2's stack already exists in extra -> must
            # still be reached and succeed.
            return path == "/extra/relion-extract-isolation/mic_000002.mrcs"

        with patch(
            "relion.protocols.protocol_extract_particles.relion.convert.Table",
            return_value=[],
        ), patch(
            "relion.protocols.protocol_extract_particles.os.path.exists",
            side_effect=exists,
        ), patch(
            "relion.protocols.protocol_extract_particles.pwutils.moveFile"
        ):
            # Must not raise.
            protocol.readPartsFromMics([mic1, mic2], outputParts)

        self.assertEqual(
            1,
            len(protocol.errors),
            "Exactly mic 1 should have failed and been reported.",
        )
        self.assertIn("micrograph 1 ", protocol.errors[0])
        self.assertEqual(
            {},
            protocol.coordDict,
            "Both mics' coordinates must be dropped after this one-shot "
            "read attempt, whether they succeeded or failed.",
        )


class _BatchValue:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value

    def set(self, value):
        self.value = value


class _BatchInputMics:
    def getObjId(self):
        return 17


class _BatchResumeHarness(ProtRelionExtractParticles):
    def __init__(self):
        self.streamingBatchSize = _BatchValue(1)

    def _setupBasicProperties(self):
        pass

    def _isStreamOpen(self):
        return False

    def isContinued(self):
        return True

    def _getStreamingBatchSize(self):
        return self.streamingBatchSize.get()

    def getInputMicrographs(self):
        return _BatchInputMics()

    def _insertFunctionStep(self, *args, **kwargs):
        return 1

    def info(self, *args, **kwargs):
        pass


class TestRelionExtractParticlesBatchResume(TestCase):
    def test_ContinuePreservesBatchSizeAfterInputStreamCloses(self):
        protocol = _BatchResumeHarness()

        with patch(
            "relion.protocols.protocol_extract_particles.ImageHandler"
        ):
            protocol._insertInitialSteps()

        self.assertEqual(
            1,
            protocol.streamingBatchSize.get(),
            "Continue must preserve the original streaming batch size even "
            "when the input stream has closed since the previous execution.",
        )
class _PendingBatchMic:
    def __init__(self, mic_id):
        self.mic_id = mic_id

    def getMicName(self):
        return "mic_%03d" % self.mic_id


class _PendingBatchHarness(ProtRelionExtractParticles):
    def __init__(self):
        self.coordsClosed = True
        self.micsClosed = True
        self.ctfsClosed = False
        self.initialIds = []
        self.micDict = {}

    def _getStreamingBatchSize(self):
        return 3


class TestRelionExtractParticlesBatchClosure(TestCase):
    def test_PartialBatchWaitsUntilAllRequiredStreamsClose(self):
        protocol = _PendingBatchHarness()
        protocol.streamClosed = protocol._isStreamClosed()

        inserted_batches = []

        def insert_single(mic, prerequisites, *args):
            raise AssertionError("batchSize=3 must not insert single-mic steps")

        def insert_batch(mics, prerequisites, *args):
            inserted_batches.append([mic.getMicName() for mic in mics])
            return 99

        deps = protocol._insertNewMics(
            [_PendingBatchMic(1), _PendingBatchMic(2)],
            lambda mic: mic.getMicName(),
            insert_single,
            insert_batch,
        )

        self.assertEqual(
            [],
            inserted_batches,
            "A partial batch must remain pending while any required input "
            "stream (for example CTF) is still open.",
        )
        self.assertEqual([], deps)
        self.assertEqual({}, protocol.micDict)

class _NoStorageFilenameSet:
    def getFileName(self):
        raise AssertionError(
            "Streaming input discovery must not depend on a storage filename."
        )


class _LogicalCoordinates(_NoStorageFilenameSet):
    def __init__(self, micrographs):
        self._micrographs = micrographs

    def getMicrographs(self):
        return self._micrographs


class _LogicalPointer:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _BackendIndependentCheckHarness(ProtRelionExtractParticles):
    def __init__(self):
        self._micrographs = _NoStorageFilenameSet()
        self._coordinates = _LogicalCoordinates(self._micrographs)
        self.inputCoordinates = _LogicalPointer(self._coordinates)
        self.micDict = {}
        self.loadCalls = 0

    def _micsOther(self):
        return False

    def _useCTF(self):
        return False

    def _loadInputList(self):
        self.loadCalls += 1
        return {}

    def _getFirstJoinStep(self):
        return None

    def debug(self, *args, **kwargs):
        pass

    def updateSteps(self):
        raise AssertionError("No new micrographs should have been scheduled.")


class TestRelionExtractParticlesBackendIndependentInputCheck(TestCase):
    def test_CheckNewInputDoesNotDependOnStorageMtime(self):
        protocol = _BackendIndependentCheckHarness()

        protocol._checkNewInput()

        self.assertEqual(
            1,
            protocol.loadCalls,
            "Streaming must refresh logical input state directly instead of "
            "gating discovery on SQLite/file modification times.",
        )

class _LogicalMic:
    def __init__(self, objId, name):
        self._objId = objId
        self._name = name

    def getObjId(self):
        return self._objId

    def getMicName(self):
        return self._name

    def clone(self):
        return _LogicalMic(self._objId, self._name)


class _LogicalMicSet:
    def __init__(self, items, closed=False):
        self._items = list(items)
        self._closed = closed
        self.loadCalls = 0

    def getFileName(self):
        raise AssertionError(
            "Logical streaming sets must not be reopened from a storage filename."
        )

    def loadAllProperties(self):
        self.loadCalls += 1

    def iterItems(self, *args, **kwargs):
        return iter(self._items)

    def isStreamClosed(self):
        return self._closed


class _LogicalCoordsForLoad:
    def __init__(self, micrographs):
        self._micrographs = micrographs

    def getMicrographs(self):
        return self._micrographs

    def getFileName(self):
        raise AssertionError(
            "Coordinates must be consumed through the logical Set API."
        )


class _LogicalLoadHarness(ProtRelionExtractParticles):
    def __init__(self):
        self.micDict = {}
        self.coordDict = {}
        self._mics = _LogicalMicSet(
            [_LogicalMic(7, "mic_007")],
            closed=True,
        )
        self._coords = _LogicalCoordsForLoad(self._mics)
        self.inputCoordinates = _LogicalPointer(self._coords)

    def _micsOther(self):
        return False

    def _useCTF(self):
        return False

    def _loadInputCoords(self, micDict):
        self.coordsClosed = True
        return micDict

    def debug(self, *args, **kwargs):
        pass


class TestRelionExtractParticlesLogicalSetLoading(TestCase):
    def test_LoadInputListUsesLogicalSets(self):
        protocol = _LogicalLoadHarness()

        newMics = protocol._loadInputList()

        self.assertEqual(["mic_007"], list(newMics.keys()))
        self.assertEqual(1, protocol._mics.loadCalls)
        self.assertTrue(protocol.micsClosed)
        self.assertTrue(protocol.ctfsClosed)
        self.assertTrue(protocol.streamClosed)
