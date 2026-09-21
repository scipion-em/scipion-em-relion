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
