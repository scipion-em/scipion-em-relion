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

""" Reading back the coordinates that relion_autopick writes.

One STAR file per micrograph, named after the micrograph id. These used to
be read through convert_deprecated, with xmipp rows. No RELION binaries,
no GPU and no datasets needed.
"""

import os

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestOutput
from pwem.objects import SetOfCoordinates, SetOfMicrographs, Micrograph

import relion.convert as convert
from ..protocols import ProtRelionAutopickLoG


AUTOPICK_STAR = """
# version 50001

data_

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2
_rlnClassNumber #3
_rlnAutopickFigureOfMerit #4
_rlnAnglePsi #5
   93.090909    31.030303            1     0.442270   300.000000
  473.917355    31.030303            1     0.431193   115.000000
  194.644628    39.493113            2     1.021585   145.000000
"""

# picking without references writes neither a class nor an angle
LOG_STAR = """
data_

loop_
_rlnCoordinateX #1
_rlnCoordinateY #2
_rlnAutopickFigureOfMerit #3
  100.000000   200.000000     0.500000
  300.000000   400.000000     0.600000
"""


class TestRelionReadCoordinates(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _micrographs(self, count):
        fn = self.getOutputPath('%s_mics.sqlite' % self._testMethodName)
        pwutils.cleanPath(fn)
        mics = SetOfMicrographs(filename=fn)
        mics.setSamplingRate(7.08)
        micList = []
        for i in range(count):
            mic = Micrograph(location='mic_%06d.mrc' % (i + 1))
            mic.setObjId(i + 1)
            mic.setMicName('mic_%06d.mrc' % (i + 1))
            mics.append(mic)
            micList.append(mic)
        mics.write()

        return mics, micList

    def _coordinates(self, mics):
        fn = self.getOutputPath('%s_coords.sqlite' % self._testMethodName)
        pwutils.cleanPath(fn)
        coordSet = SetOfCoordinates(filename=fn)
        coordSet.setMicrographs(mics)

        return coordSet

    def _starFiles(self, contents):
        paths = []
        for i, content in enumerate(contents):
            fn = self.getOutputPath('%s_mic_%06d_autopick.star'
                                    % (self._testMethodName, i + 1))
            if content is not None:
                with open(fn, 'w') as f:
                    f.write(content)
            paths.append(fn)

        return paths

    def test_theFunctionIsThereUnderItsOldName(self):
        """ protocol_autopick calls the module level function, not the
        reader's method of the same name. """
        self.assertTrue(callable(convert.readSetOfCoordinates))
        self.assertTrue(callable(convert.readCoordinates))

    def test_readsPositionsAndExtraLabels(self):
        mics, micList = self._micrographs(1)
        coordSet = self._coordinates(mics)

        convert.readSetOfCoordinates(coordSet, self._starFiles([AUTOPICK_STAR]),
                                     micList)

        coordSet.write()
        self.assertEqual(3, coordSet.getSize())
        # iterating a Set hands back the same instance over and over, so
        # the values have to be taken out of it as it goes
        rows = [(c.getPosition(), c.getMicId(), c.getMicName(),
                 round(c._rlnAutopickFigureOfMerit.get(), 5),
                 c._rlnClassNumber.get(), round(c._rlnAnglePsi.get(), 5))
                for c in coordSet]

        self.assertEqual(((93, 31), 1, 'mic_000001.mrc', 0.44227, 1, 300.0),
                         rows[0])
        self.assertEqual([(93, 31), (473, 31), (194, 39)],
                         [r[0] for r in rows])

    def test_eachFileGoesWithItsMicrograph(self):
        mics, micList = self._micrographs(2)
        coordSet = self._coordinates(mics)

        convert.readSetOfCoordinates(
            coordSet, self._starFiles([AUTOPICK_STAR, AUTOPICK_STAR]), micList)

        self.assertEqual(6, coordSet.getSize())
        self.assertEqual([1, 1, 1, 2, 2, 2],
                         sorted(c.getMicId() for c in coordSet))

    def test_aLabelThatIsNotThereIsNotInvented(self):
        mics, micList = self._micrographs(1)
        coordSet = self._coordinates(mics)

        convert.readSetOfCoordinates(coordSet, self._starFiles([LOG_STAR]),
                                     micList)

        coordSet.write()
        rows = [(c.getPosition(), round(c._rlnAutopickFigureOfMerit.get(), 5),
                 c.hasAttribute('_rlnClassNumber'),
                 c.hasAttribute('_rlnAnglePsi')) for c in coordSet]

        self.assertEqual([((100, 200), 0.5, False, False),
                          ((300, 400), 0.6, False, False)], rows)

    def test_aMissingFileDoesNotStopTheRest(self):
        """ Picking is run per batch and a micrograph can come out empty. """
        mics, micList = self._micrographs(2)
        coordSet = self._coordinates(mics)

        convert.readSetOfCoordinates(
            coordSet, self._starFiles([None, AUTOPICK_STAR]), micList)

        self.assertEqual(3, coordSet.getSize())
        self.assertEqual([2, 2, 2], [c.getMicId() for c in coordSet])

    def test_aBrokenFileDoesNotStopTheRest(self):
        mics, micList = self._micrographs(2)
        coordSet = self._coordinates(mics)

        convert.readSetOfCoordinates(
            coordSet, self._starFiles(['not a star file at all\n',
                                       AUTOPICK_STAR]), micList)

        self.assertEqual(3, coordSet.getSize())

    def test_withoutAMicListItUsesTheOnesOfTheSet(self):
        mics, _ = self._micrographs(1)
        coordSet = self._coordinates(mics)

        convert.readSetOfCoordinates(coordSet, self._starFiles([AUTOPICK_STAR]))

        coordSet.write()
        self.assertEqual(3, coordSet.getSize())
        self.assertEqual([1, 1, 1], [c.getMicId() for c in coordSet])


class TestRelionAutopickReadsItsOutput(BaseTest):
    """ The protocol's own hook, which is what broke when the module level
    readSetOfCoordinates went away with convert_deprecated: nothing in the
    suite called it without launching relion_autopick first.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def test_readCoordsFromMics(self):
        workingDir = self.getOutputPath('autopick')
        pwutils.makePath(os.path.join(workingDir, 'tmp'))

        prot = ProtRelionAutopickLoG()
        prot.setWorkingDir(workingDir)

        micList = []
        for i in range(2):
            mic = Micrograph(location='mic_%06d.mrc' % (i + 1))
            mic.setObjId(i + 1)
            mic.setMicName('mic_%06d.mrc' % (i + 1))
            micList.append(mic)
            with open(prot._getTmpPath('mic_%06d_autopick.star'
                                       % mic.getObjId()), 'w') as f:
                f.write(AUTOPICK_STAR)

        fn = self.getOutputPath('autopick_coords.sqlite')
        pwutils.cleanPath(fn)
        coordSet = SetOfCoordinates(filename=fn)

        prot.readCoordsFromMics(workingDir, micList, coordSet)

        self.assertEqual(6, coordSet.getSize())
        self.assertEqual([1, 1, 1, 2, 2, 2],
                         sorted(c.getMicId() for c in coordSet))
