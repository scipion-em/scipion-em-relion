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

""" What crop/resize volumes asks relion_image_handler to do.

The arguments are what this protocol is, so they are checked directly:
running the binary would only prove that relion_image_handler works. No
RELION binaries and no GPU needed.
"""

from pwem.objects import Volume
from pyworkflow.tests import BaseTest, setupTestOutput

from ..protocols import ProtRelionResizeVolume
from ..protocols.protocol_crop_resize import FLIP_NO, FLIP_AXIS


class TestRelionResizeVolumeArgs(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _protocol(self, **kwargs):
        vol = Volume()
        vol.setSamplingRate(1.35)

        prot = ProtRelionResizeVolume()
        prot.inputVolumes.set(vol)
        for name, value in kwargs.items():
            getattr(prot, name).set(value)

        return prot

    def test_theSamplingRateIsAlwaysPassed(self):
        args = self._protocol()._getArgDict()

        self.assertEqual({' --angpix': 1.35}, args)

    def test_mirroringAddsTheFlagOfThatAxis(self):
        for axis in ('X', 'Y', 'Z'):
            args = self._protocol(flipAxis=FLIP_AXIS.index(axis))._getArgDict()

            self.assertIn(' --flip%s' % axis, args)
            # a flag, so it must not be given a value
            self.assertEqual('', args[' --flip%s' % axis])

    def test_noMirroringLeavesNoFlipFlag(self):
        args = self._protocol(flipAxis=FLIP_NO)._getArgDict()

        self.assertEqual([], [k for k in args if 'flip' in k])

    def test_mirroringGoesTogetherWithTheRest(self):
        args = self._protocol(flipAxis=FLIP_AXIS.index('Z'),
                              doRescale=True, rescaleSamplingRate=2.7,
                              doResize=True, resizeSize=128)._getArgDict()

        self.assertEqual({' --angpix': 1.35,
                          ' --rescale_angpix': 2.7,
                          ' --new_box': 128,
                          ' --flipZ': ''}, args)

    def test_mirroringAloneIsEnoughToRun(self):
        """ It used to take a rescale or a resize to get past _validate. """
        prot = self._protocol(flipAxis=FLIP_AXIS.index('Z'))

        self.assertEqual([], prot._validate())

    def test_doingNothingIsStillRefused(self):
        prot = self._protocol()

        self.assertEqual(["You have to select at least one option!"],
                         prot._validate())
