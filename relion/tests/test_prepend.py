#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Mikel Iceta (miceta@cnb.csic.es)                     [1]
# *
# * [1] Unidad de  Biocomputacion, Centro Nacional de Biotecnologia, CSIC (CNB-CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

""" What gets prepended to the RELION commands.

RELION_PREPEND and RELION_PREPEND_CPU are what a cluster needs to load
its modules before any relion command. They used to be applied in
ProtRelionBase, which left out every protocol that does not inherit from
it, so the tests below walk those too. No RELION binaries needed.
"""

from pyworkflow.plugin import Domain
from pyworkflow.tests import BaseTest, setupTestOutput

from relion import Plugin
from ..constants import RELION_PREPEND, RELION_PREPEND_CPU
from ..protocols import (ProtRelionClassify2D, ProtRelionResizeVolume,
                         ProtRelionCreateMask3D, ProtRelionExpandSymmetry,
                         ProtRelionSymmetrizeVolume,
                         ProtRelionCompressEstimateGain,
                         ProtRelionCompressMoviesTasks)


GPU_CMD = 'module load relion-gpu &&'
CPU_CMD = 'module load relion-cpu &&'


class TestRelionPrepend(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        # Protocol._getEnviron() reaches the plugin through _package, which
        # is only filled in once the protocols have been registered, as a
        # real Scipion session does on startup
        Domain.getProtocols()

    def setUp(self):
        self._saved = {name: Plugin._vars.get(name, None)
                       for name in (RELION_PREPEND, RELION_PREPEND_CPU)}
        Plugin._vars[RELION_PREPEND] = GPU_CMD
        Plugin._vars[RELION_PREPEND_CPU] = CPU_CMD

    def tearDown(self):
        Plugin._vars.update(self._saved)

    def test_theVariablesAreDeclaredByThePlugin(self):
        """ #426: reading them straight out of os.environ kept them out of
        the plugin's variables, so nothing listed or documented them. """
        self.assertIn(RELION_PREPEND, Plugin._vars)
        self.assertIn(RELION_PREPEND_CPU, Plugin._vars)

    def test_gpuAndCpuAreTakenFromTheirOwnVariable(self):
        self.assertEqual(GPU_CMD, Plugin.getEnviron(useGpu=True).getPrepend())
        self.assertEqual(CPU_CMD, Plugin.getEnviron(useGpu=False).getPrepend())

    def test_theCpuOneIsTheDefault(self):
        """ Whoever does not say anything is not using the GPU. """
        self.assertEqual(CPU_CMD, Plugin.getEnviron().getPrepend())

    def test_nothingSetMeansNothingPrepended(self):
        Plugin._vars[RELION_PREPEND] = ''
        Plugin._vars[RELION_PREPEND_CPU] = ''

        self.assertEqual('', Plugin.getEnviron(useGpu=True).getPrepend())
        self.assertEqual('', Plugin.getEnviron().getPrepend())

    def test_anUnsetVariableIsNotTheStringNone(self):
        Plugin._vars[RELION_PREPEND_CPU] = None

        self.assertEqual('', Plugin.getEnviron().getPrepend())

    def test_everyProtocolGetsIt(self):
        """ #307: these do not inherit from ProtRelionBase, which is where
        the prepend used to be applied, and they all run relion binaries.
        """
        for cls in (ProtRelionResizeVolume, ProtRelionCreateMask3D,
                    ProtRelionExpandSymmetry, ProtRelionSymmetrizeVolume,
                    ProtRelionCompressEstimateGain,
                    ProtRelionCompressMoviesTasks):
            self.assertEqual(CPU_CMD, cls()._getEnviron().getPrepend(),
                             "%s runs without the prepend" % cls.__name__)

    def test_aGpuProtocolGetsTheGpuOne(self):
        self.assertEqual(GPU_CMD,
                         ProtRelionClassify2D()._getEnviron().getPrepend())
