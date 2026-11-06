# **************************************************************************
# *
# * The gain is estimated with the program and options the input movies
# * need, decided where the step runs: _validate runs in the process that
# * launches the protocol, not in the one that executes its steps.
# *
# **************************************************************************

import unittest
from unittest import mock

import pwem.emlib as emlib

from relion.protocols.protocol_estimate_gain import ProtRelionCompressEstimateGain


class _Value:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value

    def __str__(self):
        return str(self._value)

    def __bool__(self):
        return bool(self._value)


class _Movie:
    def __init__(self, fileName):
        self._fileName = fileName

    def getFileName(self):
        return self._fileName


class _Movies:
    def __init__(self, fileName):
        self._first = _Movie(fileName)

    def getFirstItem(self):
        return self._first


class _Pointer:
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _EstimateGain(ProtRelionCompressEstimateGain):
    def __init__(self, fileName):
        # As in a fresh process: _validate never ran here.
        self.isFloat32 = False
        self.isEER = False
        self.inputMovies = _Pointer(_Movies(fileName))
        self.numberOfThreads = _Value(4)
        self.threshold = _Value(50)
        self.maxFrames = _Value(0)
        self.eerSampling = _Value(0)
        self.random = _Value(False)
        self.jobs = []

    def _getFileName(self, key):
        return 'input_movies.star'

    def getPath(self, *paths):
        return '/run/' + '/'.join(paths)

    def runJob(self, program, args, **kwargs):
        self.jobs.append((program, args))


class TestTheInputDecidesHowTheGainIsEstimated(unittest.TestCase):
    def test_EerMoviesAreNotInverted(self):
        protocol = _EstimateGain('/data/movie.eer')

        protocol.estimateGainStep()

        program, args = protocol.jobs[0]
        self.assertEqual('relion_estimate_gain', program)
        self.assertIn('--dont_invert', args)

    def test_Float32MrcMoviesAreEstimatedWhileConverting(self):
        protocol = _EstimateGain('/data/movie.mrcs')

        with mock.patch.object(emlib.image.ImageHandler, 'getDataType',
                               return_value=emlib.DT_FLOAT):
            protocol.estimateGainStep()

        program, args = protocol.jobs[0]
        self.assertEqual('relion_convert_to_tiff', program)
        self.assertIn('--estimate_gain', args)

    def test_IntegerMoviesUseTheGainEstimator(self):
        protocol = _EstimateGain('/data/movie.tif')

        protocol.estimateGainStep()

        program, args = protocol.jobs[0]
        self.assertEqual('relion_estimate_gain', program)
        self.assertNotIn('--dont_invert', args)


if __name__ == '__main__':
    unittest.main()
