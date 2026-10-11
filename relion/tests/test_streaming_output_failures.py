from unittest import TestCase

from pyworkflow.protocol.constants import STATUS_ABORTED, STATUS_RUNNING

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


class TestRelionStreamingOutputFailures(TestCase):
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


class _Status:
    """Mimics the String attribute pyworkflow keeps the status in."""

    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


class _Mov:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId


# A real hang cannot be asserted on, so the harness caps its own polling:
# reaching the cap is the failure signal.
MAX_POLLS = 50


class _PollHarness(ProtRelionCompressMoviesTasks):
    """Drives the real input-polling loop over a fake stream."""

    def __init__(self, status=STATUS_RUNNING, closeAfter=None):
        self.status = _Status(status)
        self.polls = 0
        self._failedBatches = []
        self._closeAfter = closeAfter
        self._lastInputId = 0

    def info(self, *args):
        pass

    def _getOutputIdSet(self, outputSet):
        return set()

    def _resumeWatermarkWithGaps(self, moviesSet, seenIds):
        return 0, set()

    def _discoverNewInputItems(self, moviesSet, watermarkAttr, knownIds):
        self.polls += 1
        if self.polls > MAX_POLLS:
            raise AssertionError(
                "The input generator polled %d times with nowhere to put "
                "its results: every batch it yields is compressed for "
                "nothing." % self.polls
            )

        closed = (self._closeAfter is not None
                  and self.polls >= self._closeAfter)

        return [_Mov(self.polls)], closed, closed

    def _loadLogicalSetItemsByIds(self, moviesSet, ids):
        return []

    def drain(self):
        # A tiny wait rather than zero: the protocol deliberately falls
        # back to a one-second sleep when asked for none, so that nothing
        # spins a core, and these tests would inherit it.
        return list(self._iterInputMovies(None, 'movies', waitSecs=0.001))


class TestRelionCompressPollingStopsWhenItCannotFinish(TestCase):
    """Compression is expensive and this loop is what feeds it.

    Once the output can no longer be written, or the run has been
    aborted, everything it yields from then on is compressed and thrown
    away, and the error only surfaces when the producer closes.
    """

    def test_PollingStopsOnceABatchHasFailed(self):
        harness = _PollHarness()
        realDiscover = harness._discoverNewInputItems

        def failOnThirdPoll(*args, **kwargs):
            result = realDiscover(*args, **kwargs)
            if harness.polls == 3:
                harness._failedBatches.append({'id': 'batch-1'})
            return result

        harness._discoverNewInputItems = failOnThirdPoll

        harness.drain()

        self.assertEqual(harness.polls, 3)

    def test_PollingStopsOnceTheRunIsAborted(self):
        harness = _PollHarness()
        realDiscover = harness._discoverNewInputItems

        def abortOnSecondPoll(*args, **kwargs):
            result = realDiscover(*args, **kwargs)
            if harness.polls == 2:
                harness.status = _Status(STATUS_ABORTED)
            return result

        harness._discoverNewInputItems = abortOnSecondPoll

        harness.drain()

        self.assertEqual(harness.polls, 2)

    def test_AHealthyStreamIsPolledUntilItCloses(self):
        harness = _PollHarness(closeAfter=4)

        yielded = harness.drain()

        self.assertEqual(harness.polls, 4)
        self.assertEqual(
            len(yielded),
            4,
            "Every movie a poll found must still be yielded.",
        )
