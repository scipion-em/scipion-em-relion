# **************************************************************************
# *
# * Relion movie compression in streaming: compressed movies keep their
# * file names, a failed batch does not stop the others, and the run fails
# * only when no movie could be compressed (or the output can't be saved).
# *
# **************************************************************************

import os
import tempfile
import unittest

from relion.protocols.protocol_compress_movies_tasks import ProtRelionCompressMoviesTasks


class _Movie:
    def __init__(self, objId, fileName):
        self._objId, self._fileName = objId, fileName

    def getObjId(self):
        return self._objId

    def getFileName(self):
        return self._fileName


class _Output:
    def __init__(self, fileNames=(), size=0):
        self._fileNames, self._size = list(fileNames), size

    def getUniqueValues(self, attribute, where=None):
        return list(self._fileNames)

    def getSize(self):
        return self._size


class _Compress(ProtRelionCompressMoviesTasks):
    def __init__(self, root):
        self._root = root
        self.errors = []
        self.reported = []
        self.published = []

    def _getExtraPath(self, *paths):
        return os.path.join(self._root, 'extra', *paths)

    def _outputFromBatch(self, batch):
        if batch.get('failWrite'):
            raise RuntimeError("database is gone")
        self.published.append(batch['id'])

    def error(self, message):
        self.errors.append(message)

    def warning(self, message):
        self.reported.append(message)


class TestCompressedMoviesKeepTheirNames(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmp.cleanup()

    def test_ACompressedMovieKeepsItsName(self):
        protocol = _Compress(self._tmp.name)

        name = protocol._getOutputMovieFn(_Movie(1, '/data/a/movie.mrcs'))

        self.assertEqual('movie.tif', os.path.basename(name))

    def test_OnlyASecondMovieWithTheSameNameGetsItsId(self):
        protocol = _Compress(self._tmp.name)

        first = protocol._getOutputMovieFn(_Movie(1, '/data/a/movie.mrcs'))
        second = protocol._getOutputMovieFn(_Movie(2, '/data/b/movie.mrcs'))

        self.assertEqual('movie.tif', os.path.basename(first))
        self.assertEqual('000002__movie.tif', os.path.basename(second))

    def test_ANamePublishedBeforeAContinueStaysTaken(self):
        protocol = _Compress(self._tmp.name)
        protocol.outputMovies = _Output([protocol._getExtraPath('movie.tif')])

        name = protocol._getOutputMovieFn(_Movie(2, '/data/b/movie.mrcs'))

        self.assertEqual('000002__movie.tif', os.path.basename(name))


class TestAFailedBatchDoesNotStopTheOthers(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.protocol = _Compress(self._tmp.name)
        self.protocol._startBatchBookkeeping()

    def tearDown(self):
        self._tmp.cleanup()

    def test_ABatchThatFailedToCompressIsLeftOut(self):
        self.protocol._handleProcessedBatch({'id': 1, 'error': 'tiff failed'})
        self.protocol._handleProcessedBatch({'id': 2})

        self.assertEqual([2], self.protocol.published)
        self.assertFalse(self.protocol._mustStopFeeding())

    def test_AnOutputThatCannotBeSavedStopsTheRun(self):
        self.protocol._handleProcessedBatch({'id': 1, 'failWrite': True})

        self.assertTrue(self.protocol._mustStopFeeding())
        with self.assertRaisesRegex(RuntimeError, "could not be saved"):
            self.protocol._finishCompression(_Output(size=3))

    def test_SomeFailedBatchesOnlyWarn(self):
        self.protocol._handleProcessedBatch({'id': 1, 'error': 'tiff failed'})

        self.protocol._finishCompression(_Output(size=3))

        self.assertTrue(self.protocol.reported)

    def test_NoMovieCompressedFailsTheRun(self):
        self.protocol._handleProcessedBatch({'id': 1, 'error': 'tiff failed'})

        with self.assertRaisesRegex(RuntimeError, "every movie"):
            self.protocol._finishCompression(_Output(size=0))


if __name__ == '__main__':
    unittest.main()
