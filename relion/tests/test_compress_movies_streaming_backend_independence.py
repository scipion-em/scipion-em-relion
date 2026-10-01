# **************************************************************************
# * Regression tests for backend-independent streaming behavior.
# **************************************************************************

from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

from relion.protocols.protocol_compress_movies_tasks import ProtRelionCompressMoviesTasks


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
    def getFileName(self):
        raise AssertionError(
            "Streaming must not inspect the compatibility SQLite filename."
        )

    def isPostgresqlRuntimeOutput(self):
        raise AssertionError(
            "The plugin must not branch on the persistence backend."
        )

    def loadAllProperties(self):
        pass

    def iterItems(self):
        return iter(())

    def isStreamClosed(self):
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
        return '/tmp'

    def _linkGain(self):
        return None

    def _getCmd(self):
        return ''

    def _processBatch(self, batch):
        return batch

    def _outputFromBatch(self, batch):
        pass

    def _updateOutputSet(self, *args, **kwargs):
        pass




class _EmptyBatchManager:
    def __init__(self, *args, **kwargs):
        pass

    def generate(self):
        return iter(())


class _EmptyPipeline:
    def addGenerator(self, *args, **kwargs):
        return SimpleNamespace(outputQueue=None)

    def addProcessor(self, *args, **kwargs):
        return SimpleNamespace(outputQueue=None)

    def run(self):
        pass


class TestRelionCompressMoviesBackendIndependence(TestCase):
    def testCompressMoviesDoesNotInspectPersistenceBackend(self):
        movies = _LogicalMovies()
        protocol = _ProtocolHarness(movies)

        module = 'relion.protocols.protocol_compress_movies_tasks'
        with patch(module + '.BatchManager', _EmptyBatchManager), \
                patch(module + '.Pipeline', _EmptyPipeline):
            protocol._processAllMoviesStep()
