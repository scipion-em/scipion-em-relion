from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

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


class _LogicalMovies:
    def __init__(self):
        self.loadCalls = 0
        self.iterCalls = 0
        self.closedChecks = 0

    def getFileName(self):
        raise AssertionError(
            "Streaming must not inspect the Set storage filename to detect "
            "logical input changes."
        )

    def loadAllProperties(self):
        self.loadCalls += 1

    def iterItems(self):
        self.iterCalls += 1
        return iter(())

    def isStreamClosed(self):
        self.closedChecks += 1
        return True


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
    def addGenerator(self, *args, **kwargs):
        return SimpleNamespace(outputQueue=None)

    def addProcessor(self, *args, **kwargs):
        return SimpleNamespace(outputQueue=None)

    def run(self):
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
        self.assertGreater(
            movies.iterCalls,
            0,
            "Streaming should discover input items through Set.iterItems().",
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
