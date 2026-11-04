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

""" Optics group bookkeeping.

These tests cover the places where the optics group assignment can be
silently lost or silently changed. They need neither RELION binaries nor a
GPU, since they exercise the protocol helpers directly.
"""

from emtable import Table

from pyworkflow.object import Integer
from pyworkflow.tests import BaseTest, setupTestOutput
from pyworkflow.utils import cleanPath
from pwem.objects import (Movie, Micrograph, SetOfMicrographs,
                          SetOfCoordinates, SetOfParticles, Particle,
                          Acquisition)

from ..convert import createWriter, readSetOfParticles

from ..constants import SAME_AS_PICKING, OTHER
from ..convert.convert31 import OpticsGroups
from ..protocols import (ProtRelionMotioncor, ProtRelionAssignOpticsGroup,
                         ProtRelionExtractParticles)


def opticsStar(nGroups, first=1):
    """ An optics table with the given number of groups, numbered from
    'first' on: a single group does not have to be group 1. """
    rows = "\n".join("afis_group_%d  %d  200.0  1.4  0.1  0.885  3710  2"
                     % (i, i) for i in range(first, first + nGroups))
    return """
data_optics

loop_
_rlnOpticsGroupName #1
_rlnOpticsGroup #2
_rlnVoltage #3
_rlnSphericalAberration #4
_rlnAmplitudeContrast #5
_rlnMicrographPixelSize #6
_rlnImageSize #7
_rlnImageDimensionality #8
""" + rows + "\n"


class TestRelionMotioncorOptics(BaseTest):
    """ The base class only copies objId, micName and the filename when it
    builds the output micrographs, so motion correction has to carry the
    optics group over itself. Otherwise everything downstream falls back to
    the first optics group without any warning.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def test_micKeepsTheMovieOpticsGroup(self):
        movie = Movie(location='movie_0001.tiff')
        movie._rlnOpticsGroup = Integer(3)
        mic = Micrograph(location='movie_0001_aligned_mic.mrc')

        ProtRelionMotioncor._setOpticsGroup(movie, mic)

        self.assertEqual(3, mic.getAttributeValue('_rlnOpticsGroup'))

    def test_ungroupedMovieLeavesMicUntouched(self):
        """ Without any assignment the micrograph must not get a made up
        group: the readers already default to the first one. """
        movie = Movie(location='movie_0001.tiff')
        mic = Micrograph(location='movie_0001_aligned_mic.mrc')

        ProtRelionMotioncor._setOpticsGroup(movie, mic)

        self.assertIsNone(mic.getAttributeValue('_rlnOpticsGroup', None))
        self.assertEqual(1, mic.getAttributeValue('_rlnOpticsGroup', 1))


class TestRelionAssignOpticsWarnings(BaseTest):
    """ Items that do not match the input STAR file are left out of the
    output set, so a fully unmatched input produces an empty output. That
    has to be reported.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _warningsOf(self, matched, skipped):
        prot = ProtRelionAssignOpticsGroup()
        messages = []
        prot.warning = messages.append
        prot._warnAboutSkipped({'matched': matched, 'skipped': skipped},
                               '/data/afis.star')
        return messages

    def test_quietWhenEverythingMatches(self):
        self.assertEqual([], self._warningsOf(matched=10, skipped=0))

    def test_warnsAboutSkippedItems(self):
        messages = self._warningsOf(matched=7, skipped=3)

        self.assertEqual(1, len(messages))
        self.assertIn("3 out of 10", messages[0])

    def test_warnsAboutEmptyOutput(self):
        messages = self._warningsOf(matched=0, skipped=10)

        self.assertEqual(1, len(messages))
        self.assertIn("EMPTY", messages[0])
        self.assertIn("/data/afis.star", messages[0])


class TestRelionExtractOpticsWarnings(BaseTest):
    """ When extracting from a different set of micrographs, the particles
    take the optics groups of that set. Warn when it does not agree with
    the set the coordinates were picked on.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _outputDb(self, name):
        """ A fresh sqlite file per test: the mapper takes the column
        layout from an existing database, and the sets built here differ in
        whether the items carry _rlnOpticsGroup or not. """
        fn = self.getOutputPath('%s_%s.sqlite' % (self._testMethodName, name))
        cleanPath(fn)
        return fn

    def _createSetOfMics(self, name, groups, nOpticsGroups):
        """ A set of micrographs with the given per-item assignment, and an
        optics table with nOpticsGroups groups (0 to leave the set without
        any optics information at all). """
        mics = SetOfMicrographs(filename=self._outputDb(name))
        mics.setSamplingRate(0.885)
        mics.setDim((3710, 3838, 1))
        mics.setAcquisition(Acquisition(voltage=200, sphericalAberration=1.4,
                                        amplitudeContrast=0.1,
                                        magnification=50000))

        for i, ogNumber in enumerate(groups):
            mic = mics.ITEM_TYPE()
            # the micName is what matches the two sets, the file name is not
            mic.setFileName(self.getOutputPath('%s_%02d.mrc' % (name, i + 1)))
            mic.setMicName('mic_%02d.mrc' % (i + 1))
            mic.setSamplingRate(0.885)
            if ogNumber is not None:
                mic._rlnOpticsGroup = Integer(ogNumber)
            mics.append(mic)

        if nOpticsGroups:
            OpticsGroups.fromString(opticsStar(nOpticsGroups)).toImages(mics)
        mics.write()

        return mics

    def _warningsOf(self, pickedGroups, pickedN, otherGroups, otherN,
                    micsSource=OTHER):
        label = 'other' if micsSource == OTHER else 'same'
        picked = self._createSetOfMics('picked_%s' % label,
                                       pickedGroups, pickedN)
        other = self._createSetOfMics('other_%s' % label, otherGroups, otherN)

        coords = SetOfCoordinates(filename=self._outputDb('coords'))
        coords.setMicrographs(picked)
        coords.write()

        prot = ProtRelionExtractParticles()
        prot.downsampleType.set(micsSource)
        prot.inputCoordinates.set(coords)
        prot.inputMicrographs.set(other)

        return prot._warnings()

    def test_quietWhenMicsAreTheOnesUsedForPicking(self):
        """ Nothing to compare, the particles follow the picked micrographs. """
        self.assertEqual([], self._warningsOf([1, 2, 3], 3, [None] * 3, 0,
                                              micsSource=SAME_AS_PICKING))

    def test_quietWhenBothSetsAgree(self):
        self.assertEqual([], self._warningsOf([1, 2, 3], 3, [1, 2, 3], 3))

    def test_quietWhenNeitherSetIsGrouped(self):
        self.assertEqual([], self._warningsOf([None] * 3, 0, [None] * 3, 0))

    def test_warnsWhenOtherSetWasNeverGrouped(self):
        """ The usual accident: the micrographs to extract from did not go
        through assign optics groups. """
        messages = self._warningsOf([1, 2, 3], 3, [None] * 3, 0)

        self.assertEqual(3, len(messages))
        self.assertIn("do not declare the same optics groups", messages[0])
        self.assertIn("afis_group_3", messages[0])
        self.assertIn("2 out of 3 micrographs", messages[1])
        self.assertIn("will take the optics groups", messages[2])

    def test_warnsWhenAssignmentDiffers(self):
        """ Same optics table on both sides, shuffled assignment. """
        messages = self._warningsOf([1, 2, 3], 3, [3, 1, 2], 3)

        self.assertEqual(2, len(messages))
        self.assertIn("3 out of 3 micrographs", messages[0])
        self.assertIn("will take the optics groups", messages[1])


class TestRelionWriterOpticsGroup(BaseTest):
    """ Items without an optics group of their own have to be written into
    a group that the optics table actually declares. A single group can be
    numbered anything, so the old hardcoded 1 pointed deposited sets at a
    group that was not there and RELION refused the file.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _writeParticles(self, nOpticsGroups, groups, first=1):
        """ Write a set with the given optics table and the given per-item
        assignment, and give back the rows of the output STAR file. """
        fn = self.getOutputPath('%s.sqlite' % self._testMethodName)
        cleanPath(fn)
        parts = SetOfParticles(filename=fn)
        parts.setSamplingRate(0.885)
        parts.setAcquisition(Acquisition(voltage=200, sphericalAberration=1.4,
                                         amplitudeContrast=0.1,
                                         magnification=50000))

        for i, ogNumber in enumerate(groups):
            part = Particle()
            part.setObjId(i + 1)
            part.setLocation(i + 1, 'particles.mrcs')
            part.setSamplingRate(0.885)
            if ogNumber is not None:
                part._rlnOpticsGroup = Integer(ogNumber)
            parts.append(part)

        OpticsGroups.fromString(opticsStar(nOpticsGroups,
                                           first=first)).toImages(parts)
        parts.write()

        starFile = self.getOutputPath('%s.star' % self._testMethodName)
        createWriter().writeSetOfParticles(parts, starFile)

        return list(Table.iterRows('particles@' + starFile))

    def test_theOnlyGroupIsUsedWhateverItsNumber(self):
        """ One group, numbered 2, and no per-item column: every particle
        has to end up in group 2. """
        rows = self._writeParticles(1, [None] * 3, first=2)

        self.assertEqual([2, 2, 2], [row.rlnOpticsGroup for row in rows])

    def test_thePerItemGroupStillWins(self):
        rows = self._writeParticles(3, [1, 3, 2])

        self.assertEqual([1, 3, 2], [row.rlnOpticsGroup for row in rows])

    def test_withSeveralGroupsUngroupedItemsFallBackToOne(self):
        """ Nothing better to go on, so keep the old answer.

        Note a set cannot mix the two cases: the sqlite layout is taken
        from the first item, so an item carrying _rlnOpticsGroup after one
        that does not would lose it.
        """
        rows = self._writeParticles(3, [None] * 3)

        self.assertEqual([1, 1, 1], [row.rlnOpticsGroup for row in rows])


class TestRelionOpticsGroupRoundTrip(BaseTest):
    """ The optics group of each particle has to survive the trip: read it
    from a STAR file, keep it on the item, and write the same value back.
    Until rlnOpticsGroup was added to PARTICLE_EXTRA_LABELS it was read and
    dropped, so an imported set could only ever be written into one group.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _particlesStar(self, groups):
        """ A particles STAR file with a per-item rlnOpticsGroup column. """
        fn = self.getOutputPath('%s_in.star' % self._testMethodName)
        rows = "\n".join("particles.mrcs  %d  %d" % (i + 1, og)
                          for i, og in enumerate(groups))
        with open(fn, 'w') as f:
            f.write(opticsStar(max(groups)))
            f.write("""
data_particles

loop_
_rlnImageName #1
_rlnImageId #2
_rlnOpticsGroup #3
""" + rows + "\n")

        return fn

    def _read(self, groups):
        starFile = self._particlesStar(groups)
        fn = self.getOutputPath('%s.sqlite' % self._testMethodName)
        cleanPath(fn)
        parts = SetOfParticles(filename=fn)
        readSetOfParticles(starFile, parts)
        parts.write()

        return parts

    def test_readKeepsTheGroupOfEveryParticle(self):
        parts = self._read([3, 1, 2])

        self.assertEqual([3, 1, 2], [p.getAttributeValue('_rlnOpticsGroup')
                                     for p in parts.iterItems()])

    def test_writeGivesBackTheSameGroups(self):
        parts = self._read([3, 1, 2])
        starFile = self.getOutputPath('%s_out.star' % self._testMethodName)
        createWriter().writeSetOfParticles(parts, starFile)

        rows = list(Table.iterRows('particles@' + starFile))

        self.assertEqual([3, 1, 2], [row.rlnOpticsGroup for row in rows])
