# **************************************************************************
# *
# * Building RELION picks a CUDA host compiler nvcc accepts, and builds
# * without the GUI when the Xft files it needs are missing.
# *
# **************************************************************************

import subprocess
import unittest
from unittest import mock

from relion import build_toolchain as tc


class TestTheGccEachCudaAccepts(unittest.TestCase):
    def test_KnownReleases(self):
        self.assertEqual(11, tc.maxGccForCuda((11, 8)))
        self.assertEqual(10, tc.maxGccForCuda((11, 2)))
        self.assertEqual(12, tc.maxGccForCuda((12, 0)))
        self.assertEqual(13, tc.maxGccForCuda((12, 6)))
        self.assertEqual(14, tc.maxGccForCuda((12, 9)))

    def test_UnknownOrMissingCuda(self):
        self.assertIsNone(tc.maxGccForCuda(None))
        self.assertIsNone(tc.maxGccForCuda((9, 2)))


def _fakeRun(outputs):
    def run(cmd):
        return outputs.get(cmd[0], '')
    return run


class TestTheHostCompilerChoice(unittest.TestCase):
    def _choose(self, cudaRelease, defaultGcc, installed):
        outputs = {'/cuda/bin/nvcc': 'Cuda compilation tools, release %s, V0' % cudaRelease,
                   'cc': defaultGcc}
        which = lambda name: '/usr/bin/' + name if name in installed else None
        with mock.patch.object(tc, '_run', _fakeRun(outputs)):
            return tc.cudaHostCompiler('/cuda/bin/nvcc', 'cc', which=which)

    def test_ATooNewDefaultGccIsReplacedByTheNewestOneCudaAccepts(self):
        # Ubuntu 24.04 with CUDA 11.8: the RELION build failed with
        # "unsupported GNU version! gcc versions later than 11".
        self.assertEqual('/usr/bin/gcc-11',
                         self._choose('11.8', '13', {'gcc-11', 'gcc-12', 'gcc-9'}))

    def test_AnAcceptedDefaultGccIsLeftAlone(self):
        self.assertIsNone(self._choose('12.4', '13', {'gcc-11'}))

    def test_WithoutASuitableGccNothingChanges(self):
        self.assertIsNone(self._choose('11.8', '13', {'gcc-12'}))

    def test_WithoutCudaNothingChanges(self):
        with mock.patch.object(tc, '_run', _fakeRun({})):
            self.assertIsNone(tc.cudaHostCompiler(None, 'cc'))
            self.assertEqual('', tc.cmakeToolchainFlags(None, which=lambda n: None))


class TestTheGuiFlag(unittest.TestCase):
    def _flag(self, env):
        return subprocess.run(['/bin/sh', '-c', 'echo ' + tc.GUI_FLAG_SHELL],
                              capture_output=True, text=True, env=env).stdout.strip()

    def test_WithoutPkgConfigTheGuiIsLeftOut(self):
        self.assertEqual('-DGUI=OFF', self._flag({'PATH': '/nonexistent'}))


if __name__ == '__main__':
    unittest.main()
