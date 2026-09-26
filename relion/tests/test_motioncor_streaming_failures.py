from unittest import TestCase
from unittest.mock import patch

from relion.protocols.protocol_motioncor import ProtRelionMotioncor


class _Value:
    def __init__(self, value=None):
        self.value = value

    def get(self):
        return self.value

    def hasValue(self):
        return self.value not in (None, "")


class _InputMoviesSet:
    def getGain(self):
        return None

    def getSamplingRate(self):
        return 1.0


class _Pointer:
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _Movie:
    def __init__(self, index=1):
        self._fileName = "/data/movie_%03d.mrcs" % index

    def getFileName(self):
        return self._fileName

    def setFileName(self, fileName):
        self._fileName = fileName

    def getSamplingRate(self):
        return 1.0


class _Writer:
    def writeSetOfMovies(self, *args, **kwargs):
        pass


class _MotioncorHarness(ProtRelionMotioncor):
    def __init__(self):
        self.inputMovies = _Pointer(_InputMoviesSet())
        self.binFactor = 1.0
        self.bfactor = 150
        self.patchX = 5
        self.patchY = 5
        self.groupFrames = 1
        self.numberOfThreads = 1
        self.defectFile = _Value(None)
        self.doDW = False
        self.isEER = False
        self.saveFloat16 = False
        self.extraParams = _Value(None)
        self.runCalls = 0
        self.errors = []

    def _getOutputMovieFolder(self, movie):
        return "/tmp/relion-motioncor-test"

    def _getMovieRoot(self, movie):
        return "movie_001"

    def _getFrameRange(self, n=None, prefix=None):
        return 1, 2

    def _savePsSum(self):
        return False

    def _runProgram(self, *args, **kwargs):
        self.runCalls += 1
        if self.runCalls == 1:
            raise RuntimeError("relion_run_motioncorr failed")

    def _saveAlignmentPlots(self, *args, **kwargs):
        pass

    def _computeExtra(self, *args, **kwargs):
        pass

    def _moveFiles(self, *args, **kwargs):
        pass

    def error(self, message, *args, **kwargs):
        self.errors.append(message)


class TestRelionMotioncorStreamingFailures(TestCase):
    def test_MotioncorProcessingFailureDoesNotAbortFollowingMovie(self):
        protocol = _MotioncorHarness()
        firstMovie = _Movie(1)
        secondMovie = _Movie(2)

        with patch(
            "relion.protocols.protocol_motioncor.pwutils.makePath"
        ), patch(
            "relion.protocols.protocol_motioncor.OpticsGroups.fromImages",
            return_value=object(),
        ), patch(
            "relion.protocols.protocol_motioncor.convert.createWriter",
            return_value=_Writer(),
        ):
            protocol._processMovie(firstMovie)
            protocol._processMovie(secondMovie)

        self.assertEqual(
            2,
            protocol.runCalls,
            "A failed movie must not prevent the following movie from being processed.",
        )
        self.assertTrue(
            any(
                "ERROR processing movie" in message
                for message in protocol.errors
            ),
            "The failed movie should still be reported in the protocol log.",
        )
