from unittest import TestCase

from relion.protocols.protocol_compress_movies_tasks import (
    ProtRelionCompressMoviesTasks,
)


class _Movie:
    def __init__(self, fileName="/data/movie_001.eer"):
        self._fileName = fileName
        self.acquisition = None
        self.framesRange = None

    def getFileName(self):
        return self._fileName

    def setAcquisition(self, acquisition):
        self.acquisition = acquisition

    def setFramesRange(self, framesRange):
        self.framesRange = framesRange


class _OutputMovies:
    def __init__(self):
        self.appended = []
        self.acquisition = object()
        self.framesRange = (1, 40, 1)

    def enableAppend(self):
        pass

    def getAcquisition(self):
        return self.acquisition

    def getFramesRange(self):
        return self.framesRange

    def append(self, movie):
        self.appended.append(movie)


class _ProtocolHarness(ProtRelionCompressMoviesTasks):
    def __init__(self):
        self._outputMovies = _OutputMovies()
        self.outputMovies = self._outputMovies
        self.outputUpdates = []

    def _updateOutputSet(self, outputName, outputSet, state):
        self.outputUpdates.append((outputName, outputSet, state))


class TestRelionPostgresqlStreamingFailures(TestCase):
    def test_FailedCompressBatchIsNotPersistedAsSuccessfulOutput(self):
        protocol = _ProtocolHarness()
        movie = _Movie()

        protocol._outputFromBatch({
            "id": "batch-1",
            "items": [movie],
            "error": "relion_convert_to_tiff failed",
        })

        self.assertEqual(
            [],
            protocol._outputMovies.appended,
            "Movies from a failed compression batch must not be persisted "
            "as successful output because Resume uses persisted output as "
            "the processed-item blacklist.",
        )
