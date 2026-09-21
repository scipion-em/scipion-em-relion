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
    def __init__(self):
        self._fileName = "/data/movie_001.mrcs"

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

    def _getOutputMovieFolder(self, movie):
        return "/tmp/relion-motioncor-test"

    def _getMovieRoot(self, movie):
        return "movie_001"

    def _getFrameRange(self, n=None, prefix=None):
        return 1, 2

    def _savePsSum(self):
        return False

    def _runProgram(self, *args, **kwargs):
        raise RuntimeError("relion_run_motioncorr failed")

    def error(self, *args, **kwargs):
        pass


class TestRelionMotioncorStreamingFailures(TestCase):
    def test_MotioncorProcessingFailurePropagatesBeforeDoneCheckpoint(self):
        protocol = _MotioncorHarness()
        movie = _Movie()

        with patch(
            "relion.protocols.protocol_motioncor.pwutils.makePath"
        ), patch(
            "relion.protocols.protocol_motioncor.OpticsGroups.fromImages",
            return_value=object(),
        ), patch(
            "relion.protocols.protocol_motioncor.convert.createWriter",
            return_value=_Writer(),
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "relion_run_motioncorr failed",
            ):
                protocol._processMovie(movie)
