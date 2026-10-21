# **************************************************************************
# *
# * EER options are decided from the movies being processed: _validate runs
# * in the process that launches the protocol, and flags it sets are gone
# * in the one that executes the steps.
# *
# **************************************************************************

import os
import tempfile
import unittest
from unittest import mock

from pwem.protocols import ProtAlignMovies

from relion.protocols.protocol_compress_movies_tasks import ProtRelionCompressMoviesTasks
from relion.protocols.protocol_motioncor import ProtRelionMotioncor


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value

    def __int__(self):
        return int(self._value)

    def __index__(self):
        return int(self._value)


class _Movie:
    def __init__(self, objId, fileName):
        self._objId, self._fileName = objId, fileName

    def getObjId(self):
        return self._objId

    def getFileName(self):
        return self._fileName

    def setFileName(self, fileName):
        self._fileName = fileName


class _Compress(ProtRelionCompressMoviesTasks):
    def __init__(self, root):
        # As in a fresh process: _validate never ran here.
        self.isEER = False
        self._root = root
        self._gainFile = None
        self.eerGroup = _Value(32)
        self.eerSampling = _Value(0)
        self.deflateLevel = _Value(1)
        self.commands = []

    def getEnumText(self, name):
        return 'lzw'

    def _getExtraPath(self, *paths):
        return os.path.join(self._root, 'extra', *paths)

    def _runProgram(self, cmd, **kwargs):
        self.commands.append(cmd)

    def info(self, message):
        pass

    def error(self, message):
        pass


class TestCompressionKnowsItsMoviesAreEer(unittest.TestCase):
    def _compress(self, fileName):
        with tempfile.TemporaryDirectory() as root:
            protocol = _Compress(root)
            batchPath = os.path.join(root, 'batch')
            os.makedirs(batchPath)
            protocol._processBatch({'id': 1, 'path': batchPath,
                                    'items': [_Movie(1, '/data/' + fileName)]})
            return protocol.commands[0]

    def test_EerMoviesAreCompressedWithTheirGrouping(self):
        cmd = self._compress('movie.eer')

        self.assertIn('--eer_grouping 32', cmd)
        self.assertIn('--eer_upsampling 1', cmd)

    def test_OtherMoviesAreNot(self):
        self.assertNotIn('--eer_grouping', self._compress('movie.tif'))


class _Movies:
    def __init__(self, first):
        self._first = first

    def getFirstItem(self):
        return self._first


class _Pointer:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class TestMotionCorrectionKnowsItsMoviesAreEer(unittest.TestCase):
    def test_TheRunKnowsItWhenItStarts(self):
        protocol = ProtRelionMotioncor.__new__(ProtRelionMotioncor)
        protocol.isEER = False
        protocol.inputMovies = _Pointer(_Movies(_Movie(1, '/data/movie.eer')))

        with mock.patch.object(ProtAlignMovies, '_insertAllSteps'), \
                mock.patch.object(ProtRelionMotioncor, '_prepareFailedForRetry'):
            ProtRelionMotioncor._insertAllSteps(protocol)

        self.assertTrue(protocol.isEER)

    def test_EachMovieSaysItself(self):
        self.assertTrue(ProtRelionMotioncor._isEERMovie(_Movie(1, '/d/m.EER')))
        self.assertFalse(ProtRelionMotioncor._isEERMovie(_Movie(1, '/d/m.tif')))
        self.assertFalse(ProtRelionMotioncor._isEERMovie(None))


if __name__ == '__main__':
    unittest.main()
