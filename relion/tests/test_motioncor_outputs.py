# **************************************************************************
# *
# * Relion motion correction in streaming: a movie is published to every
# * output or to none, a failed one is never published, read once per run,
# * and aligned again on Continue.
# *
# **************************************************************************

import os
import tempfile
import unittest

from relion.protocols.protocol_motioncor import ProtRelionMotioncor


class _Movie:
    def __init__(self, objId):
        self._objId = objId

    def getObjId(self):
        return self._objId


class _Output:
    def __init__(self, ids):
        self._ids = set(ids)

    def getIdSet(self):
        return set(self._ids)


class _Flag:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value

    def __bool__(self):
        return bool(self._value)


class _Motioncor(ProtRelionMotioncor):
    def __init__(self, root, doDW=True, saveNonDW=True):
        self._root = root
        self.doDW = _Flag(doDW)
        self.saveNonDW = _Flag(saveNonDW)
        self.updates = []
        self.failed = []
        self.listOfMovies = []
        self.streamClosed = False
        os.makedirs(os.path.join(root, 'extra'), exist_ok=True)

    def _getExtraPath(self, *paths):
        return os.path.join(self._root, 'extra', *paths)

    def _getOutputMicName(self, movie):
        return 'movie_%06d_mic.mrc' % movie.getObjId()

    def _getOutputMicWtName(self, movie):
        return 'movie_%06d_mic_DW.mrc' % movie.getObjId()

    def _writeResults(self, movie):
        for name in (self._getOutputMicName(movie), self._getOutputMicWtName(movie)):
            open(self._getExtraPath(name), 'w').close()

    def _updateOutputMovieSet(self, newDone, streamMode):
        self.updates.append(('outputMovies', [m.getObjId() for m in newDone]))

    def _updateOutputMicSet(self, newDone, sqliteFn, getOutputMicName,
                            outputName, streamMode):
        self.updates.append((outputName, [m.getObjId() for m in newDone]))

    def _writeFailedList(self, movieList):
        self.failed.extend(m.getObjId() for m in movieList)

    def _getFirstJoinStep(self):
        return None

    def debug(self, message):
        pass


class TestEveryOutputGetsEveryMovie(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmp.cleanup()

    def test_AMovieIsPublishedOnlyWhenEveryOutputHasIt(self):
        protocol = _Motioncor(self._tmp.name)
        protocol.outputMovies = _Output([1, 2])
        protocol.outputMicrographs = _Output([1, 2])
        protocol.outputMicrographsDoseWeighted = _Output([1])

        self.assertEqual({1}, ProtRelionMotioncor._getPublishedMovieIds(protocol))

    def test_EachOutputGetsOnlyTheMoviesItLacks(self):
        protocol = _Motioncor(self._tmp.name)
        protocol.outputMovies = _Output([1, 2])
        protocol.outputMicrographs = _Output([1, 2])
        protocol.outputMicrographsDoseWeighted = _Output([1])

        protocol._updateOutputSets([_Movie(2), _Movie(3)], 'open')

        self.assertEqual([('outputMovies', [3]),
                          ('outputMicrographs', [3]),
                          ('outputMicrographsDoseWeighted', [2, 3])],
                         protocol.updates)


class _Step:
    def __init__(self, movieId):
        self.funcName = 'processMovieStep'
        self.argsStr = '[{"object.id": %d}, false]' % movieId

    def isFinished(self):
        return True


class TestAFailedMovieIsNeverPublished(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmp.cleanup()

    def _protocol(self):
        protocol = _Motioncor(self._tmp.name)
        protocol.listOfMovies = [_Movie(1), _Movie(2)]
        protocol._steps = [_Step(1), _Step(2)]
        protocol._writeResults(_Movie(1))
        return protocol

    def test_OnlyMoviesWithTheirResultsArePublished(self):
        protocol = self._protocol()

        ProtRelionMotioncor._checkNewOutput(protocol)

        self.assertEqual(('outputMovies', [1]), protocol.updates[0])
        self.assertEqual([2], protocol.failed)

    def test_AFailedMovieIsLookedAtOnce(self):
        protocol = self._protocol()

        ProtRelionMotioncor._checkNewOutput(protocol)
        protocol.updates = []
        ProtRelionMotioncor._checkNewOutput(protocol)

        self.assertEqual([], protocol.updates)
        self.assertEqual([2], protocol.failed)

    def test_TheRunEndsWhenEveryMovieWasLookedAt(self):
        protocol = self._protocol()
        protocol.streamClosed = True

        ProtRelionMotioncor._checkNewOutput(protocol)

        self.assertTrue(protocol.finished)

    def test_AFailedMovieLosesItsDoneMarker(self):
        protocol = self._protocol()
        markers = {}
        for movieId in (1, 2):
            markers[movieId] = protocol._getExtraPath('DONE', 'movie_%06d.TXT' % movieId)
            os.makedirs(os.path.dirname(markers[movieId]), exist_ok=True)
            open(markers[movieId], 'w').close()
        protocol._getMovieDone = lambda movie: markers[movie.getObjId()]

        for movieId in (1, 2):
            protocol._forgetDoneMarkerIfFailed(_Movie(movieId))

        self.assertTrue(os.path.exists(markers[1]))
        self.assertFalse(os.path.exists(markers[2]))


if __name__ == '__main__':
    unittest.main()
