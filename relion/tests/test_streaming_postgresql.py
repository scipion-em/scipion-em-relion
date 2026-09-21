# **************************************************************************
# * Regression tests for PostgreSQL streaming compatibility.
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


class _PostgresqlMovies:
    def getFileName(self):
        return '/tmp/compatibility.sqlite'

    def isPostgresqlRuntimeOutput(self):
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


class _RecordingSetMonitor:
    sources = []

    def __init__(self, setClass, source, *args, **kwargs):
        self.sources.append(source)

    def iterProtocolInput(self, *args, **kwargs):
        return iter(())


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


class TestRelionPostgresqlStreaming(TestCase):
    def test_PostgresqlCompressMoviesDoesNotMonitorCompatibilitySqlite(self):
        movies = _PostgresqlMovies()
        protocol = _ProtocolHarness(movies)
        _RecordingSetMonitor.sources = []

        module = 'relion.protocols.protocol_compress_movies_tasks'
        with patch(module + '.SetMonitor', _RecordingSetMonitor), \
                patch(module + '.BatchManager', _EmptyBatchManager), \
                patch(module + '.Pipeline', _EmptyPipeline):
            protocol._processAllMoviesStep()

        self.assertEqual(
            [],
            _RecordingSetMonitor.sources,
            'PostgreSQL streaming must not be discovered through the '
            'compatibility SQLite filename.'
        )
