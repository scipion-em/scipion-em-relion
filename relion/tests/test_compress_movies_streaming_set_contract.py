from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

from .logical_set_fakes import LogicalSetFake
from relion.protocols.protocol_compress_movies_tasks import (
    ProtRelionCompressMoviesTasks,
)


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _Pointer:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _LogicalMovies(LogicalSetFake):
    """Empty, already-closed input that must never be asked for a file."""

    def __init__(self):
        super().__init__([], streamClosed=True)
        self.iterCalls = 0

    def iterItems(self, *args, **kwargs):
        self.iterCalls += 1
        return super().iterItems(*args, **kwargs)


class _ProtocolHarness(ProtRelionCompressMoviesTasks):
    def __init__(self, movies):
        self.inputMovies = _Pointer(movies)
        self.streamingSleepOnWait = _Value(0)
        self.streamingBatchSize = _Value(1)
        self.numberOfThreads = _Value(0)

    def info(self, *args, **kwargs):
        pass

    def _runProgram(self, *args, **kwargs):
        pass

    def _getTmpPath(self, *args):
        return "/tmp"

    def _linkGain(self):
        return None

    def _getCmd(self):
        return ""

    def _processBatch(self, batch):
        return batch

    def _outputFromBatch(self, batch):
        pass

    def _updateOutputSet(self, *args, **kwargs):
        pass


class _EmptyBatchManager:
    def __init__(self, batchSize, itemsIter, *args, **kwargs):
        self._itemsIter = itemsIter

    def generate(self):
        list(self._itemsIter)
        return iter(())


class _EmptyPipeline:
    # A real Pipeline's generator node fully drains its generator
    # function inside pipe.run() (TaskGenerator.run() iterates it to
    # completion before notifyGeneratorEnds()). This fake must do the
    # same instead of being a no-op, otherwise nothing ever actually
    # iterates batchMgr.generate()/the underlying movies iterator in
    # these tests.
    def __init__(self):
        self._generatorFunc = None

    def addGenerator(self, generatorFunc, *args, **kwargs):
        self._generatorFunc = generatorFunc
        return SimpleNamespace(outputQueue=None)

    def addProcessor(self, *args, **kwargs):
        return SimpleNamespace(outputQueue=None)

    def run(self):
        if self._generatorFunc is not None:
            for _ in self._generatorFunc():
                pass


class TestRelionCompressMoviesStreamingSetContract(TestCase):
    def testCompressMoviesUsesLogicalSetStateInsteadOfStorageFile(self):
        movies = _LogicalMovies()
        protocol = _ProtocolHarness(movies)

        module = "relion.protocols.protocol_compress_movies_tasks"
        with patch(module + ".BatchManager", _EmptyBatchManager), \
                patch(module + ".Pipeline", _EmptyPipeline):
            protocol._processAllMoviesStep()

        self.assertGreater(
            movies.loadCalls,
            0,
            "Streaming should refresh the logical Set state.",
        )
        # Discovery asks the Set for the ids above the watermark rather
        # than walking it, so an empty input is never iterated at all.
        self.assertEqual(
            0,
            movies.fullScans,
            "Streaming must not walk the whole input Set to discover items.",
        )
        self.assertGreater(
            movies.closedChecks,
            0,
            "Streaming should determine completion through Set.isStreamClosed().",
        )

class TestRelionCompressMoviesFailedBatchCompletion(TestCase):
    def testFailedBatchFailsProtocolAfterPipelineDrains(self):
        class _FailedBatchProtocol(_ProtocolHarness):
            def __init__(self, movies):
                super().__init__(movies)
                self.numberOfThreads = _Value(1)
                self.closedOutput = False

            def _processBatch(self, batch):
                batch["error"] = "simulated compression failure"
                return batch

            def _outputFromBatch(self, batch):
                return ProtRelionCompressMoviesTasks._outputFromBatch(
                    self, batch
                )

            def _updateOutputSet(self, outputName, outputSet, state):
                self.closedOutput = True

        class _FailurePipeline:
            def __init__(self):
                self.processors = []

            def addGenerator(self, *args, **kwargs):
                return SimpleNamespace(outputQueue=object())

            def addProcessor(self, inputQueue, processor, outputQueue=None):
                self.processors.append(processor)
                return SimpleNamespace(outputQueue=object())

            def run(self):
                batch = {
                    "id": "batch-1",
                    "index": 1,
                    "items": [object()],
                    "path": "/tmp/batch-1",
                }
                for processor in self.processors:
                    batch = processor(batch)

        movies = _LogicalMovies()
        protocol = _FailedBatchProtocol(movies)

        module = "relion.protocols.protocol_compress_movies_tasks"
        with patch(module + ".BatchManager", _EmptyBatchManager), \
                patch(module + ".Pipeline", _FailurePipeline):
            with self.assertRaisesRegex(
                RuntimeError,
                "batch",
            ):
                protocol._processAllMoviesStep()

        self.assertFalse(
            protocol.closedOutput,
            "A failed compression batch must prevent normal output closure.",
        )

    def testLaterBatchesContinueAfterEarlierBatchFails(self):
        processed = []

        class _MixedBatchProtocol(_ProtocolHarness):
            def __init__(self, movies):
                super().__init__(movies)
                self.numberOfThreads = _Value(1)

            def _processBatch(self, batch):
                processed.append(batch["id"])
                if batch["id"] == "batch-1":
                    batch["error"] = "simulated compression failure"
                return batch

            def _outputFromBatch(self, batch):
                return None

        class _MixedPipeline:
            def __init__(self):
                self.processors = []

            def addGenerator(self, *args, **kwargs):
                return SimpleNamespace(outputQueue=object())

            def addProcessor(self, inputQueue, processor, outputQueue=None):
                self.processors.append(processor)
                return SimpleNamespace(outputQueue=object())

            def run(self):
                batches = [
                    {
                        "id": "batch-1",
                        "index": 1,
                        "items": [object()],
                        "path": "/tmp/batch-1",
                    },
                    {
                        "id": "batch-2",
                        "index": 2,
                        "items": [object()],
                        "path": "/tmp/batch-2",
                    },
                ]

                for batch in batches:
                    current = batch
                    for processor in self.processors:
                        current = processor(current)

        movies = _LogicalMovies()
        protocol = _MixedBatchProtocol(movies)

        module = "relion.protocols.protocol_compress_movies_tasks"
        with patch(module + ".BatchManager", _EmptyBatchManager), \
                patch(module + ".Pipeline", _MixedPipeline):
            with self.assertRaisesRegex(
                RuntimeError,
                "batch",
            ):
                protocol._processAllMoviesStep()

        self.assertEqual(
            ["batch-1", "batch-2"],
            processed,
            "A failed streaming batch must not stop later batches from "
            "being processed before the protocol reports the failure.",
        )


class _MissingTiffMovie:
    def __init__(self, fileName="/data/movie_001.mrcs"):
        self._fileName = fileName

    def getFileName(self):
        return self._fileName

    def setFileName(self, fileName):
        self._fileName = fileName


class _FakeStarFile:
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, traceback):
        return False

    def writeTable(self, *args, **kwargs):
        pass


class _FakeTable:
    def __init__(self, *args, **kwargs):
        pass

    def addRowValues(self, *args, **kwargs):
        pass


class _MissingTiffHarness(ProtRelionCompressMoviesTasks):
    def __init__(self):
        self.cmd = ""

    def info(self, *args, **kwargs):
        pass

    def error(self, *args, **kwargs):
        pass

    def _runProgram(self, *args, **kwargs):
        pass

    def _getExtraPath(self, *paths):
        return "/extra/" + "/".join(paths)


class TestRelionCompressMoviesMissingTiff(TestCase):
    def testMissingTiffMarksBatchAsFailed(self):
        protocol = _MissingTiffHarness()
        movie = _MissingTiffMovie()
        batch = {
            "id": "batch-missing-tiff",
            "items": [movie],
            "path": "/tmp/relion-missing-tiff",
        }

        module = "relion.protocols.protocol_compress_movies_tasks"
        with patch(module + ".StarFile", _FakeStarFile), \
                patch(module + ".Table", _FakeTable), \
                patch(module + ".pwutils.createLink"), \
                patch(module + ".pwutils.envVarOn", return_value=True), \
                patch(module + ".os.path.exists", return_value=False):
            result = protocol._processBatch(batch)

        self.assertIn(
            "error",
            result,
            "A missing expected TIFF must fail the batch so the protocol "
            "cannot close successfully while silently dropping a movie.",
        )
        self.assertRegex(
            result["error"],
            "TIFF|tiff",
            "The batch error must specifically report the missing TIFF, "
            "not an unrelated harness failure.",
        )

class _PersistedMovie:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId


class _ClosedResumeInput(LogicalSetFake):
    def __init__(self, items):
        super().__init__(items, streamClosed=True)


class _PersistedOutput(LogicalSetFake):
    def __init__(self, items):
        super().__init__(items, streamClosed=False)


class _ResumeCompletionHarness(ProtRelionCompressMoviesTasks):
    def __init__(self, inputMovies, outputMovies):
        self.inputMovies = _Pointer(inputMovies)
        self.outputMovies = outputMovies
        self.streamingSleepOnWait = _Value(0)
        self.streamingBatchSize = _Value(1)
        self.numberOfThreads = _Value(0)
        self.closedWith = None

    def info(self, *args, **kwargs):
        pass

    def _runProgram(self, *args, **kwargs):
        pass

    def _getTmpPath(self, *args):
        return "/tmp"

    def _linkGain(self):
        return None

    def _getCmd(self):
        return ""

    def _updateOutputSet(self, outputName, outputSet, state):
        self.closedWith = (outputName, outputSet, state)


class TestRelionCompressMoviesResumeCompletion(TestCase):
    def testResumeClosesPersistedOutputWhenNoNewMoviesRemain(self):
        movies = [_PersistedMovie(1), _PersistedMovie(2)]
        inputMovies = _ClosedResumeInput(movies)
        outputMovies = _PersistedOutput(
            [_PersistedMovie(1), _PersistedMovie(2)]
        )
        protocol = _ResumeCompletionHarness(inputMovies, outputMovies)

        module = "relion.protocols.protocol_compress_movies_tasks"
        with patch(module + ".BatchManager", _EmptyBatchManager),                 patch(module + ".Pipeline", _EmptyPipeline):
            protocol._processAllMoviesStep()

        self.assertIsNotNone(
            protocol.closedWith,
            "Resume must finalize the already-persisted output even when "
            "there are no new movies to process.",
        )
        self.assertIs(
            outputMovies,
            protocol.closedWith[1],
            "Resume must close the existing persisted output, not None.",
        )


class _DoseAcquisition:
    def __init__(self, dosePerFrame=None):
        self._dosePerFrame = dosePerFrame

    def getDosePerFrame(self):
        return self._dosePerFrame

    def setDosePerFrame(self, value):
        self._dosePerFrame = value


class _DoseOutputMovies:
    def __init__(self):
        self._acquisition = _DoseAcquisition(dosePerFrame=None)
        self.appended = []

    def setStreamState(self, state):
        pass

    def copyInfo(self, other):
        pass

    def setDim(self, dim):
        pass

    def getAcquisition(self):
        return self._acquisition

    def setFramesRange(self, r):
        pass

    def getFramesRange(self):
        return [1, 1, 1]

    def setGain(self, gain):
        pass

    def enableAppend(self):
        pass

    def append(self, movie):
        self.appended.append(movie)


class _DoseMovie:
    def getDim(self):
        return (10, 10, 2)

    def getFileName(self):
        return "/data/movie_0001.tif"

    def setAcquisition(self, acq):
        pass

    def setFramesRange(self, r):
        pass


class _DoseBatchHarness(ProtRelionCompressMoviesTasks):
    # Regression harness: the input movies' acquisition may not carry a
    # dose per frame at all (e.g. an import that didn't set it).
    # _outputFromBatch used to do "acq.getDosePerFrame() * eerGroup"
    # unconditionally, crashing with TypeError on None - and because this
    # runs inside an emtools Pipeline worker thread with no exception
    # boundary of its own, that crash would silently hang the whole
    # pipeline instead of failing cleanly.
    def __init__(self):
        self.eerGroup = _Value(32)
        self._outputMovies = None
        self.inputMovies = _Pointer(None)

    def _createSetOfMovies(self):
        return _DoseOutputMovies()

    def _getExtraPath(self, *paths):
        return "/extra/" + "/".join(paths)

    def _updateOutputSet(self, *args, **kwargs):
        pass

    def _defineSourceRelation(self, *args, **kwargs):
        pass


class TestRelionCompressMoviesMissingDoseRegression(TestCase):
    def testOutputFromBatchDoesNotCrashWhenDoseIsUnknown(self):
        protocol = _DoseBatchHarness()
        movie = _DoseMovie()
        batch = {"id": "batch-1", "items": [movie]}

        # Must not raise.
        protocol._outputFromBatch(batch)

        self.assertIsNone(
            protocol._outputMovies.getAcquisition().getDosePerFrame(),
            "An unknown dose per frame must stay None (honestly unknown), "
            "not crash and not be fabricated into a fake 0.0 dose.",
        )
        self.assertEqual(
            [movie],
            protocol._outputMovies.appended,
            "The movie must still be appended to the output even though "
            "its dose could not be regrouped.",
        )


class _CapturingPipeline:
    """ Fake Pipeline that records the last processor added (the
    _updateOutput stage) and actually invokes it during run(), unlike
    _EmptyPipeline above which only records wiring. """
    def __init__(self):
        self.lastProcessor = None

    def addGenerator(self, *args, **kwargs):
        return SimpleNamespace(outputQueue=None)

    def addProcessor(self, inputQueue, processor, outputQueue=None):
        self.lastProcessor = processor
        return SimpleNamespace(outputQueue=None)

    def run(self):
        if self.lastProcessor is not None:
            self.lastProcessor({"id": "batch-1", "items": []})


class _CrashingOutputHarness(_ProtocolHarness):
    def __init__(self, movies):
        super().__init__(movies)
        self.errors = []

    def _outputFromBatch(self, batch):
        raise RuntimeError("boom from _outputFromBatch")

    def error(self, message, *args, **kwargs):
        self.errors.append(message)


class TestRelionCompressMoviesUpdateOutputPipelineSafety(TestCase):
    def testOutputFromBatchExceptionFailsTheBatchInsteadOfHangingThePipeline(self):
        movies = _LogicalMovies()
        protocol = _CrashingOutputHarness(movies)

        module = "relion.protocols.protocol_compress_movies_tasks"
        with patch(module + ".BatchManager", _EmptyBatchManager), \
                patch(module + ".Pipeline", _CapturingPipeline):
            with self.assertRaises(RuntimeError) as ctx:
                protocol._processAllMoviesStep()

        # The RuntimeError must be the controlled, post-pipe.run() failure
        # (proving _updateOutput caught the exception and recorded it in
        # failedBatches), not the raw "boom" escaping uncaught.
        self.assertIn(
            "streaming batches",
            str(ctx.exception),
            "_outputFromBatch raising must be converted into a failed "
            "batch by _updateOutput, not propagate out of the pipeline "
            "processor thread (which would hang it instead of failing).",
        )
