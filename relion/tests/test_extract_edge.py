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

""" Discarding the particles at the edge of the micrograph.

relion has no option for this: a box that sticks out of the micrograph is
filled with the nearest pixels inside it, which smears lines along that
edge of the particle. The protocol can leave those coordinates out before
handing them over, and these tests pin down which ones go. They need no
RELION binaries and no GPU.
"""

from pwem.objects import Micrograph, Coordinate
from pwem.emlib.image import ImageHandler
from pyworkflow.tests import BaseTest, setupTestOutput

from ..protocols import ProtRelionExtractParticles


MIC_X, MIC_Y = 200, 160
BOX = 20


class TestRelionExtractDiscardEdge(BaseTest):
    @classmethod
    def setUpClass(cls):
        """ One real micrograph: getDim() reads the file header, there is
        no setter for it. """
        setupTestOutput(cls)
        cls.micFn = cls.getOutputPath('mic_01.mrc')
        ImageHandler().createEmptyImage(cls.micFn, xDim=MIC_X, yDim=MIC_Y)

    def _micrograph(self):
        mic = Micrograph(location=self.micFn)
        mic.setMicName('mic_01.mrc')

        return mic

    def _protocol(self, boxSize=BOX):
        prot = ProtRelionExtractParticles()
        prot.boxSize.set(boxSize)
        # set by _setupBasicProperties() in a real run; here the
        # coordinates are already in the pixels of this micrograph
        prot._getPos = lambda coord: coord.getPosition()
        self.messages = []
        prot.info = self.messages.append

        return prot

    def _discard(self, positions, boxSize=BOX):
        coords = []
        for i, (x, y) in enumerate(positions):
            coord = Coordinate()
            coord.setObjId(i + 1)
            coord.setPosition(x, y)
            coords.append(coord)

        kept = self._protocol(boxSize)._discardEdgeCoords(self._micrograph(),
                                                          coords)

        return [coord.getPosition() for coord in kept]

    def test_keepsTheOnesThatFit(self):
        positions = [(100, 80), (50, 40), (150, 120)]

        self.assertEqual(positions, self._discard(positions))
        self.assertEqual([], self.messages)

    def test_dropsTheOnesOverEachEdge(self):
        kept = self._discard([(100, 80),    # middle
                              (9, 80),      # over the left edge
                              (191, 80),    # over the right edge
                              (100, 9),     # over the bottom edge
                              (100, 151)])  # over the top edge

        self.assertEqual([(100, 80)], kept)
        self.assertEqual(1, len(self.messages))
        self.assertIn("4 of 5 coordinates", self.messages[0])
        self.assertIn("10 px", self.messages[0])

    def test_theBoxIsAllowedToTouchTheEdge(self):
        """ Half a box away is exactly where it still fits. """
        kept = self._discard([(10, 10), (MIC_X - 10, MIC_Y - 10)])

        self.assertEqual([(10, 10), (190, 150)], kept)

    def test_anOddBoxSizeIsNotRoundedIntoTheEdge(self):
        """ half is 10.5 px for a box of 21, so 10 is out and 11 is in. """
        kept = self._discard([(10, 80), (11, 80)], boxSize=21)

        self.assertEqual([(11, 80)], kept)

    def test_aBoxBiggerThanTheMicrographDropsEverything(self):
        """ Nothing fits, and _validate() is the one that stops this. """
        self.assertEqual([], self._discard([(100, 80)], boxSize=400))

    def test_validateRefusesABoxBiggerThanTheMicrographs(self):
        prot = self._protocol(boxSize=400)
        prot.doDiscardEdge.set(True)
        prot.getInputMicrographs = lambda: _MicSetStub()

        errors = [e for e in prot._validate() if 'larger than' in e]

        self.assertEqual(1, len(errors))
        self.assertIn("400 px", errors[0])
        self.assertIn("200 x 160", errors[0])

    def test_validateIsQuietWhenTheBoxFits(self):
        prot = self._protocol(boxSize=20)
        prot.doDiscardEdge.set(True)
        prot.getInputMicrographs = lambda: _MicSetStub()

        self.assertEqual([], [e for e in prot._validate()
                              if 'larger than' in e])


class _MicSetStub:
    """ Only getDim() is reached by _validate(). """
    def getDim(self):
        return MIC_X, MIC_Y, 1
