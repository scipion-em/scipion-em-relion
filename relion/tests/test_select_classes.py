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

""" Matching the 2D class ranker output back to the input classes.

The lookup goes through the class average locations, because the width of
the index RELION writes in 'index@stack' is not stable: it went from 6 to
8 digits in RELION 5.0.1 (FILENAMENUMBERLENGTH). These tests pin both
widths down. They need neither RELION binaries nor a GPU.
"""

import os

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestOutput
from pwem.objects import Class2D, Particle

from ..protocols import ProtRelionSelectClasses2D


CLASS_STACK = 'Class2D/job008/run_it025_classes.mrcs'
# Digits relion_class_ranker pads the class average index to
RELION_50_DIGITS = 6
RELION_51_DIGITS = 8


def classAveragesStar(classes, digits):
    """ A class_averages.star like the one relion_class_ranker writes, with
    one row per (class number, predicted score, estimated resolution). """
    rows = "\n".join(
        "%0*d@%s  %f  0.154511  0.0  0.0  %f" % (digits, number, CLASS_STACK,
                                                 score, resolution)
        for number, score, resolution in classes)
    return """# version 30001

data_

loop_
_rlnReferenceImage #1
_rlnPredictedClassScore #2
_rlnClassDistribution #3
_rlnAccuracyRotations #4
_rlnAccuracyTranslationsAngst #5
_rlnEstimatedResolution #6
""" + rows + "\n"


def class2d(number):
    """ A class whose representative sits at the given position of the
    class averages stack, the way the classification protocol sets it. """
    representative = Particle()
    representative.setLocation(number, CLASS_STACK)
    cls2d = Class2D()
    cls2d.setRepresentative(representative)

    return cls2d


class TestRelionSelectClassesLookup(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _protocol(self, classes, digits, minResolution=25):
        """ A protocol whose extra/ holds the given class_averages.star. """
        workingDir = self.getOutputPath(self._testMethodName)
        pwutils.makePath(os.path.join(workingDir, 'extra'))

        prot = ProtRelionSelectClasses2D()
        prot.setWorkingDir(workingDir)
        prot.minResolution.set(minResolution)

        with open(prot._getExtraPath('class_averages.star'), 'w') as f:
            f.write(classAveragesStar(classes, digits))

        return prot

    def _rowOf(self, prot, number):
        return prot._getClassRow(prot._loadClassRows(), class2d(number))

    def test_findsTheClassesRelion51Wrote(self):
        """ From 5.0.1 on the index is padded to 8 digits. """
        prot = self._protocol([(23, 0.68, 16.18), (39, 0.64, 18.88)],
                              RELION_51_DIGITS)

        for number, score in [(23, 0.68), (39, 0.64)]:
            row = self._rowOf(prot, number)
            self.assertIsNotNone(row, "Class %d was not matched" % number)
            self.assertAlmostEqual(score, row.rlnPredictedClassScore, places=5)

    def test_findsTheClassesRelion50Wrote(self):
        """ Up to 5.0.0 it was padded to 6, which must keep working. """
        prot = self._protocol([(23, 0.68, 16.18), (39, 0.64, 18.88)],
                              RELION_50_DIGITS)

        for number in (23, 39):
            self.assertIsNotNone(self._rowOf(prot, number),
                                 "Class %d was not matched" % number)

    def test_skipsClassesTheRankerDidNotScore(self):
        """ Only the selected classes make it into class_averages.star. """
        prot = self._protocol([(23, 0.68, 16.18)], RELION_51_DIGITS)

        self.assertIsNone(self._rowOf(prot, 40))

    def test_skipsClassesBelowTheResolutionLimit(self):
        prot = self._protocol([(23, 0.68, 16.18), (39, 0.64, 30.0)],
                              RELION_51_DIGITS, minResolution=20)

        self.assertIsNotNone(self._rowOf(prot, 23))
        self.assertIsNone(self._rowOf(prot, 39))
