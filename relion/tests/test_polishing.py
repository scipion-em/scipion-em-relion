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

""" CTF pre-multiplication in bayesian polishing.

relion_motion_refine can multiply the polished particles by their CTF,
and marks the result with rlnCtfDataAreCtfPremultiplied so that the CTF is
not applied twice. These tests cover the arguments and the rules around
the option. They need no RELION binaries and no GPU.
"""

import os

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestOutput
from pwem.objects import SetOfParticles, Acquisition, Particle

from ..convert.convert31 import OpticsGroups
from ..protocols import ProtRelionBayesianPolishing


SHINY_OPTICS = """
data_optics

loop_
_rlnOpticsGroupName #1
_rlnOpticsGroup #2
_rlnVoltage #3
_rlnSphericalAberration #4
_rlnAmplitudeContrast #5
_rlnImagePixelSize #6
_rlnImageSize #7
_rlnImageDimensionality #8
_rlnCtfDataAreCtfPremultiplied #9
opticsGroup1  1  300.0  2.7  0.1  1.35  256  2  1
"""


class _Holder:
    """ Stands in for a pointer param. """
    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _VolumeStub:
    def getSamplingRate(self):
        return 1.35


class _PostprocessStub:
    """ Only what _getPolishArgs() reaches of the postprocess protocol. """
    outputVolume = _VolumeStub()

    def _getExtraPath(self, *paths):
        return os.path.join('/postprocess/extra', *paths)


class TestRelionPolishingCtfPremultiply(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _protocol(self, **kwargs):
        workingDir = self.getOutputPath(self._testMethodName)
        pwutils.makePath(workingDir)

        prot = ProtRelionBayesianPolishing()
        prot.setWorkingDir(workingDir)
        prot._createFilenameTemplates()
        prot.inputPostprocess = _Holder(_PostprocessStub())
        for name, value in kwargs.items():
            getattr(prot, name).set(value)

        return prot

    def test_nothingIsPremultipliedByDefault(self):
        args = self._protocol()._getPolishArgs()

        self.assertNotIn('--ctf_multiply', args)
        self.assertNotIn('--crop', args)

    def test_premultiplyingIsPassedOn(self):
        args = self._protocol(doCtfMultiply=True)._getPolishArgs()

        self.assertIn('--ctf_multiply', args)

    def test_noCropUnlessASizeIsGiven(self):
        args = self._protocol(doCtfMultiply=True, cropSize=-1)._getPolishArgs()

        self.assertNotIn('--crop', args)

    def test_cropGoesWithItsSize(self):
        args = self._protocol(doCtfMultiply=True,
                              cropSize=256)._getPolishArgs()

        self.assertIn('--crop 256', args)

    def test_trainingDoesNotPremultiply(self):
        """ There are no particles to write out when training. """
        prot = self._protocol(doCtfMultiply=True, cropSize=256)
        prot.operation.set(prot.OP_TRAIN)

        args = prot._getPolishArgs()

        self.assertNotIn('--ctf_multiply', args)
        self.assertNotIn('--crop', args)

    def test_anOddCropSizeIsRefused(self):
        prot = self._protocol(doCtfMultiply=True, cropSize=255,
                              extrSize=512, rescaledSize=256)

        errors = [e for e in prot._validate() if 'crop size' in e]

        self.assertEqual(1, len(errors))
        self.assertIn("even number", errors[0])

    def test_aCropLargerThanTheBoxIsRefused(self):
        prot = self._protocol(doCtfMultiply=True, cropSize=300,
                              extrSize=512, rescaledSize=256)

        errors = [e for e in prot._validate() if 'crop size' in e]

        self.assertEqual(1, len(errors))
        self.assertIn("larger than", errors[0])

    def test_aCropThatFitsIsAccepted(self):
        prot = self._protocol(doCtfMultiply=True, cropSize=128,
                              extrSize=512, rescaledSize=256)

        self.assertEqual([], [e for e in prot._validate()
                              if 'crop size' in e])

    def test_theCropRulesOnlyApplyWhenPremultiplying(self):
        prot = self._protocol(doCtfMultiply=False, cropSize=255,
                              extrSize=512, rescaledSize=256)

        self.assertEqual([], [e for e in prot._validate()
                              if 'crop size' in e])

    def test_warnsAboutWhatTheLabelDoesNotReach(self):
        quiet = self._protocol()
        loud = self._protocol(doCtfMultiply=True)

        self.assertEqual(1, len(quiet._warnings()))
        self.assertEqual(2, len(loud._warnings()))
        self.assertIn('rlnCtfDataAreCtfPremultiplied', loud._warnings()[1])


class TestRelionPolishingOutputOptics(BaseTest):
    """ The output has to carry the optics table relion wrote, not the one
    of the input particles: with --ctf_multiply the two differ by exactly
    the label that says the CTF has already been applied.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _outputSet(self, name, opticsStar):
        """ Run the protocol's own step over a shiny.star holding the given
        optics table, and give back the set it set the optics on. """
        workingDir = self.getOutputPath(name)
        pwutils.makePath(os.path.join(workingDir, 'extra'))

        prot = ProtRelionBayesianPolishing()
        prot.setWorkingDir(workingDir)
        prot._createFilenameTemplates()
        with open(prot._getFileName('shiny'), 'w') as f:
            f.write(opticsStar)

        fn = self.getOutputPath('%s.sqlite' % name)
        pwutils.cleanPath(fn)
        parts = SetOfParticles(filename=fn)
        parts.setSamplingRate(1.35)
        parts.setAcquisition(Acquisition(voltage=300, sphericalAberration=2.7,
                                         amplitudeContrast=0.1,
                                         magnification=50000))

        prot._setOutputOptics(parts)

        return parts

    def test_thePremultipliedFlagReachesTheOutputSet(self):
        parts = self._outputSet('shiny', SHINY_OPTICS)

        # read back what the set now carries, the way a writer would
        stored = OpticsGroups.fromString(
            parts.getAcquisition().opticsGroupInfo.get())

        self.assertEqual(1, stored.first().rlnCtfDataAreCtfPremultiplied)

    def test_anInputWithoutTheFlagDoesNotGainOne(self):
        parts = self._outputSet('plain', SHINY_OPTICS.replace(
            "_rlnCtfDataAreCtfPremultiplied #9\n", "").replace(
            "2  1\n", "2\n"))
        stored = OpticsGroups.fromString(
            parts.getAcquisition().opticsGroupInfo.get())

        self.assertFalse(stored.hasColumn('rlnCtfDataAreCtfPremultiplied'))


SHINY_PARTICLES = """
data_optics

loop_
_rlnOpticsGroupName #1
_rlnOpticsGroup #2
_rlnVoltage #3
_rlnSphericalAberration #4
_rlnAmplitudeContrast #5
_rlnImagePixelSize #6
_rlnImageSize #7
_rlnImageDimensionality #8
opticsGroup1  1  300.0  2.7  0.1  1.35  256  2

data_particles

loop_
_rlnImageId #1
_rlnImageName #2
%s
"""


class TestRelionPolishingOutputParticles(BaseTest):
    """ Which particles come out of polishing and where they point.

    shiny.star can hold fewer particles than the input: RELION leaves out
    the ones whose movie, or whose motion model, it did not have. The
    input items without a row have to be left out of the output.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _protocol(self, rows):
        workingDir = self.getOutputPath(self._testMethodName)
        pwutils.makePath(os.path.join(workingDir, 'extra'))

        prot = ProtRelionBayesianPolishing()
        prot.setWorkingDir(workingDir)
        prot._createFilenameTemplates()
        with open(prot._getFileName('shiny'), 'w') as f:
            f.write(SHINY_PARTICLES % '\n'.join(
                '%d 000%03d@shiny.mrcs' % (objId, i + 1)
                for i, objId in enumerate(rows)))

        return prot

    def _inputSet(self, size=4):
        fn = self.getOutputPath('%s_in.sqlite' % self._testMethodName)
        pwutils.cleanPath(fn)
        parts = SetOfParticles(filename=fn)
        parts.setSamplingRate(1.35)
        for i in range(size):
            part = Particle()
            part.setObjId(i + 1)
            part.setLocation(i + 1, 'input.mrcs')
            parts.append(part)
        parts.write()

        return parts

    def _outputSet(self, rows, size=4):
        prot = self._protocol(rows)
        inputSet = self._inputSet(size)

        fn = self.getOutputPath('%s_out.sqlite' % self._testMethodName)
        pwutils.cleanPath(fn)
        outSet = SetOfParticles(filename=fn)
        outSet.copyInfo(inputSet)

        prot._newLocations = prot._readNewLocations()
        outSet.copyItems(inputSet, updateItemCallback=prot._updatePtcl)

        return outSet

    def test_everyParticlePointsAtItsPolishedImage(self):
        outSet = self._outputSet([1, 2, 3, 4])

        self.assertEqual(4, outSet.getSize())
        self.assertEqual([(1, 'shiny.mrcs'), (2, 'shiny.mrcs'),
                          (3, 'shiny.mrcs'), (4, 'shiny.mrcs')],
                         [p.getLocation() for p in outSet])

    def test_theParticlesWithoutARowAreLeftOut(self):
        outSet = self._outputSet([2, 4])

        self.assertEqual(2, outSet.getSize())
        self.assertEqual([2, 4], [p.getObjId() for p in outSet])

    def test_theRowsDoNotHaveToBeSorted(self):
        """ The xmipp iterator this replaced walked both sides in id order,
        which meant sorting the whole table in memory first. """
        outSet = self._outputSet([4, 1, 3, 2])

        self.assertEqual(4, outSet.getSize())
        self.assertEqual([1, 2, 3, 4], [p.getObjId() for p in outSet])
        # the second row of the file is the one with id 1
        self.assertEqual((2, 'shiny.mrcs'),
                         outSet.getFirstItem().getLocation())

    def test_anIdThatIsNotInTheInputDoesNotShiftTheRest(self):
        outSet = self._outputSet([1, 2, 99, 3])

        self.assertEqual([1, 2, 3], [p.getObjId() for p in outSet])
        self.assertEqual([(1, 'shiny.mrcs'), (2, 'shiny.mrcs'),
                          (4, 'shiny.mrcs')],
                         [p.getLocation() for p in outSet])

    def test_anEmptyShinyLeavesAnEmptyOutput(self):
        outSet = self._outputSet([])

        self.assertEqual(0, outSet.getSize())
