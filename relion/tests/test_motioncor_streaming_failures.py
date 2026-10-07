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

class _DoneMovie:
    def __init__(self, objId=1):
        self._objId = objId

    def getObjId(self):
        return self._objId


class _PersistenceFailureHarness(ProtRelionMotioncor):
    def __init__(self):
        self.listOfMovies = [_DoneMovie(1)]
        self.streamClosed = False
        self.persistenceCalls = 0
        # The step for movie 1 finished; nothing has been published yet.
        self._steps = [_FinishedMovieStep(1)]

    def _writeDoneList(self, movies):
        raise AssertionError(
            "Completion must not be recorded in a DONE sidecar."
        )

    def _updateOutputSets(self, newDone, streamMode):
        self.persistenceCalls += 1
        raise RuntimeError("simulated output persistence failure")

    def debug(self, *args, **kwargs):
        pass


class _FinishedMovieStep:
    funcName = 'processMovieStep'

    def __init__(self, movieId):
        self.argsStr = '[{"object.id": %d}, false]' % movieId

    def isFinished(self):
        return True


class TestRelionMotioncorPersistenceResume(TestCase):
    def test_PersistenceFailureDoesNotCommitMovieDoneList(self):
        protocol = _PersistenceFailureHarness()

        with self.assertRaisesRegex(
            RuntimeError,
            "persistence",
        ):
            protocol._checkNewOutput()

        self.assertEqual(
            1,
            protocol.persistenceCalls,
            "The protocol must attempt to persist the processed movie.",
        )
        # There is no done list to commit any more: what counts as
        # published is the output Set itself, so a failed persist leaves
        # the movie pending and Continue retries it, with nothing to undo.
        self.assertEqual(
            {1},
            protocol._getFinishedMovieIds(),
        )
        self.assertEqual(set(), protocol._getPublishedMovieIds())


class _DoseAcquisition:
    def __init__(self, doseInitial=None, dosePerFrame=None):
        self._doseInitial = doseInitial
        self._dosePerFrame = dosePerFrame

    def getDoseInitial(self):
        return self._doseInitial

    def getDosePerFrame(self):
        return self._dosePerFrame


class _DoseInputMoviesSet(_InputMoviesSet):
    def __init__(self, acquisition):
        self._acquisition = acquisition

    def getAcquisition(self):
        return self._acquisition

    def getFramesRange(self):
        return [1, 2, 1]


class _MissingDoseMotioncorHarness(_MotioncorHarness):
    # Regression test harness: _getCorrectedDose is inherited from
    # pwem's ProtAlignMovies, which does
    # "preExp += dose * (firstFrame - 1)" unconditionally and crashes
    # with TypeError when the acquisition has no dose per frame.
    # _validate() blocks a direct launch with doDW on and no dose, but
    # that is not necessarily re-enforced on every launch path (e.g. a
    # resumed/chained workflow) - the same scenario already fixed in
    # scipion-em-motioncorr.
    def __init__(self):
        super().__init__()
        self.inputMovies = _Pointer(
            _DoseInputMoviesSet(
                _DoseAcquisition(doseInitial=None, dosePerFrame=None))
        )
        self.doDW = True
        self.eerGroup = _Value(32)
        self.saveNonDW = _Value(False)

    def _runProgram(self, *args, **kwargs):
        self.runCalls += 1


class TestRelionMotioncorMissingDoseRegression(TestCase):
    def test_ProcessMovieDoesNotCrashWhenDoseIsUnknown(self):
        protocol = _MissingDoseMotioncorHarness()
        movie = _Movie(1)

        with patch(
            "relion.protocols.protocol_motioncor.pwutils.makePath"
        ), patch(
            "relion.protocols.protocol_motioncor.OpticsGroups.fromImages",
            return_value=object(),
        ), patch(
            "relion.protocols.protocol_motioncor.convert.createWriter",
            return_value=_Writer(),
        ):
            protocol._processMovie(movie)

        self.assertEqual(
            1,
            protocol.runCalls,
            "The protocol must still attempt to run motioncor even when "
            "the dose is unknown - only the dose arguments should default "
            "safely to 0.0 instead of crashing before runJob is reached.",
        )
        self.assertEqual(
            [],
            protocol.errors,
            "A missing dose must not be treated as an argument-building "
            "failure; it should default safely instead.",
        )


class TestRelionMotioncorCalcPsDoseZeroDoseRegression(TestCase):
    def test_CalcPsDoseDoesNotCrashWhenDoseIsUnknown(self):
        protocol = _MissingDoseMotioncorHarness()
        protocol.dosePSsum = _Value(4.0)

        # Must not raise ZeroDivisionError.
        result = protocol._calcPsDose()

        self.assertEqual(1, result)
