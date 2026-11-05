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

""" Reading star files written before RELION 3.1.

Until now this went through convert_deprecated, which built an xmipp
MetaData out of the star file. These tests pin down the behaviour of the
emtable reader that replaces it, including the details the xmipp one got
by accident. They need no RELION binaries, no GPU and no datasets.
"""

import os

import numpy as np

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestOutput
from pwem.constants import ALIGN_2D, ALIGN_NONE, ALIGN_PROJ
from pwem.objects import SetOfParticles

from ..convert import createReader


# A 3.0 particles table: shifts in pixels (rlnOriginX/Y, not Angstrom),
# the acquisition repeated on every row and no optics table at all
COLUMNS = ['rlnImageName', 'rlnImageId', 'rlnEnabled',
           'rlnCoordinateX', 'rlnCoordinateY', 'rlnMicrographName',
           'rlnDefocusU', 'rlnDefocusV', 'rlnDefocusAngle',
           'rlnVoltage', 'rlnSphericalAberration', 'rlnAmplitudeContrast',
           'rlnOriginX', 'rlnOriginY', 'rlnAngleRot', 'rlnAngleTilt',
           'rlnAnglePsi', 'rlnClassNumber', 'rlnNrOfSignificantSamples']

ROWS = [
    ['000001@Particles/parts.mrcs', 1, 1,
     299.0, 48.0, 'Micrographs/mic_01.mrcs',
     35136.06, 33578.89, 100.27,
     300.0, 2.0, 0.1,
     -1.956394, 7.043606, 107.801792, 105.980954,
     -131.130655, 1, 38],
    ['000002@Particles/parts.mrcs', 2, 1,
     610.0, 54.0, 'Micrographs/mic_01.mrcs',
     35136.06, 33578.89, 100.27,
     300.0, 2.0, 0.1,
     0.043606, -1.187644, 78.074471, 117.122663,
     78.769820, 3, 16],
]


def writeStar(fileName, tableName='images', columns=COLUMNS, rows=ROWS):
    """ Write a 3.0 style star file: one table, no version comment. """
    with open(fileName, 'w') as f:
        f.write('\ndata_%s\n\nloop_\n' % tableName)
        for i, col in enumerate(columns):
            f.write('_%s #%d\n' % (col, i + 1))
        for row in rows:
            f.write(' '.join(str(v) for v in row) + '\n')

    return fileName


# convert_deprecated.rowToAlignment() on the first row above, before it
# was removed: the shifts are read as pixels, never scaled by the pixel size
PROJ_MATRIX = [
    [0.66178144, 0.68968279, -0.29390981, 3.56314856],
    [0.40270275, 0.00365581, 0.91532351, -0.76209512],
    [0.63235735, -0.72410240, -0.27531780, -6.33743216],
    [0., 0., 0., 1.]]

ALIGN2D_MATRIX = [
    [-0.65777833, 0.75321157, 0., -1.956394],
    [-0.75321157, -0.65777833, 0., 7.043606],
    [0., 0., 1., 0.],
    [0., 0., 0., 1.]]


class TestRelionReader30(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _read(self, alignType=ALIGN_PROJ, starName='particles.star',
              tableName='images', columns=COLUMNS, rows=ROWS, **kwargs):
        starFile = writeStar(self.getOutputPath(starName), tableName,
                             columns, rows)
        sqliteFn = self.getOutputPath(starName.replace('.star', '.sqlite'))
        pwutils.cleanPath(sqliteFn)
        partsSet = SetOfParticles(filename=sqliteFn)
        partsSet.setSamplingRate(1.0)

        reader = createReader(format='30', alignType=alignType, **kwargs)
        reader.readSetOfParticles(starFile, partsSet,
                                  alignType=alignType, **kwargs)

        return partsSet

    def test_theTableNameIsIgnored(self):
        """ relion wrote data_images, Scipion data_Particles and the
        tutorial data_noname: the reader takes the first table, whatever
        it is called, the way md.MetaData(fn) used to. """
        for i, name in enumerate(['images', 'Particles', 'noname']):
            partsSet = self._read(tableName=name,
                                  starName='table_%d.star' % i)

            self.assertEqual(2, partsSet.getSize())

    def test_readsEveryRow(self):
        partsSet = self._read()
        first = partsSet.getFirstItem()

        self.assertEqual(2, partsSet.getSize())
        self.assertEqual(1, first.getObjId())
        self.assertEqual(1, first.getIndex())
        self.assertEqual('Particles/parts.mrcs', first.getFileName())

    def test_readsTheClassAndTheCoordinate(self):
        partsSet = self._read()
        first = partsSet.getFirstItem()
        coord = first.getCoordinate()

        self.assertEqual(1, first.getClassId())
        self.assertEqual((299, 48), coord.getPosition())
        self.assertEqual('mic_01.mrcs', coord.getMicName())

    def test_readsTheMicrographId(self):
        """ The old reader copied rlnMicrographId onto the particle, and
        the import protocol relies on it to group by micrograph. """
        columns = COLUMNS + ['rlnMicrographId']
        rows = [list(row) + [4] for row in ROWS]

        partsSet = self._read(starName='micid.star', columns=columns,
                              rows=rows)

        self.assertEqual(4, partsSet.getFirstItem().getMicId())

    def test_noMicrographIdWhenTheColumnIsMissing(self):
        self.assertIsNone(self._read().getFirstItem().getMicId())

    def test_readsTheCtf(self):
        partsSet = self._read()
        ctf = partsSet.getFirstItem().getCTF()

        self.assertTrue(partsSet.hasCTF())
        self.assertAlmostEqual(35136.06, ctf.getDefocusU(), places=2)
        self.assertAlmostEqual(33578.89, ctf.getDefocusV(), places=2)
        self.assertAlmostEqual(100.27, ctf.getDefocusAngle(), places=2)

    def test_noCtfWithoutAllThreeDefocusValues(self):
        """ rowToCtfModel asked for U, V and the angle together, and
        setupCTF never filled the missing ones in: it looked for the xmipp
        label ctfDefocusU, which a relion row does not have. """
        keep = [c for c in COLUMNS if c not in ('rlnDefocusV',
                                                'rlnDefocusAngle')]
        rows = [[row[COLUMNS.index(c)] for c in keep] for row in ROWS]

        partsSet = self._read(starName='noctf.star', columns=keep, rows=rows)

        self.assertFalse(partsSet.hasCTF())
        self.assertFalse(partsSet.getFirstItem().hasCTF())

    def test_theProjMatrixIsTheOneOfTheOldReader(self):
        """ Golden matrix, taken from convert_deprecated.rowToAlignment()
        on these same values before it was removed. """
        partsSet = self._read(pixelSize=3.54)
        matrix = partsSet.getFirstItem().getTransform().getMatrix()

        self.assertTrue(np.allclose(PROJ_MATRIX, matrix, atol=1e-6),
                        "Not the alignment of the old reader:\n%s" % matrix)

    def test_the2DMatrixIsTheOneOfTheOldReader(self):
        partsSet = self._read(alignType=ALIGN_2D, starName='align2d.star',
                              pixelSize=3.54)
        matrix = partsSet.getFirstItem().getTransform().getMatrix()

        self.assertTrue(np.allclose(ALIGN2D_MATRIX, matrix, atol=1e-6),
                        "Not the alignment of the old reader:\n%s" % matrix)

    def test_shiftsAreInPixelsNotAngstrom(self):
        """ 3.1 renamed rlnOriginX/Y to rlnOriginX/YAngst and divides them
        by the pixel size. Reading a 3.0 file with the 3.1 reader would
        scale every shift, so the pixel size must not show up here. """
        oneToOne = self._read(starName='px1.star', pixelSize=1.0)
        scaled = self._read(starName='px2.star', pixelSize=3.54)

        self.assertTrue(np.allclose(
            oneToOne.getFirstItem().getTransform().getMatrix(),
            scaled.getFirstItem().getTransform().getMatrix(), atol=1e-6))

    def test_noAlignmentLeavesNoTransform(self):
        partsSet = self._read(alignType=ALIGN_NONE, starName='alignnone.star')

        self.assertIsNone(partsSet.getFirstItem().getTransform())
        self.assertEqual(ALIGN_NONE, partsSet.getAlignment())

    def test_theAlignmentTypeReachesTheSet(self):
        self.assertEqual(ALIGN_PROJ, self._read().getAlignment())

    def test_extraLabelsAreKept(self):
        """ rlnNrOfSignificantSamples is not one of PARTICLE_EXTRA_LABELS,
        it came in through the old IMAGE_EXTRA_LABELS list. """
        partsSet = self._read(extraLabels=['rlnNrOfSignificantSamples'])
        first = partsSet.getFirstItem()

        self.assertEqual(38, first._rlnNrOfSignificantSamples.get())

    def test_disabledRowsAreDropped(self):
        """ removeDisabled=True was a lie: xmipp's removeDisabled() looks
        at the label 'enabled', and a relion star file carries 'rlnEnabled'
        instead, so nothing was ever dropped. """
        rows = [list(row) for row in ROWS]
        rows[0][COLUMNS.index('rlnEnabled')] = 0

        partsSet = self._read(starName='disabled.star', rows=rows)

        self.assertEqual(1, partsSet.getSize())
        self.assertEqual(3, partsSet.getFirstItem().getClassId())

    def test_disabledRowsAreKeptWhenAsked(self):
        rows = [list(row) for row in ROWS]
        rows[0][COLUMNS.index('rlnEnabled')] = 0

        partsSet = self._read(starName='keepdisabled.star', rows=rows,
                              removeDisabled=False)

        self.assertEqual(2, partsSet.getSize())

    def test_theAcquisitionComesFromTheRows(self):
        """ There is no optics table before 3.1. """
        partsSet = self._read(readAcquisition=True)
        acq = partsSet.getAcquisition()

        self.assertAlmostEqual(300.0, acq.getVoltage(), places=2)
        self.assertAlmostEqual(2.0, acq.getSphericalAberration(), places=2)
        self.assertAlmostEqual(0.1, acq.getAmplitudeContrast(), places=2)

    def test_theAcquisitionOfTheSetIsLeftAloneOnImport(self):
        """ The import protocol fills it from the form instead. """
        starFile = writeStar(self.getOutputPath('acq.star'))
        sqliteFn = self.getOutputPath('acq.sqlite')
        pwutils.cleanPath(sqliteFn)
        partsSet = SetOfParticles(filename=sqliteFn)
        partsSet.setSamplingRate(1.0)
        partsSet.getAcquisition().setVoltage(200.0)

        reader = createReader(format='30', alignType=ALIGN_PROJ)
        reader.readSetOfParticles(starFile, partsSet, alignType=ALIGN_PROJ,
                                  readAcquisition=False)

        self.assertAlmostEqual(200.0, partsSet.getAcquisition().getVoltage(),
                               places=2)

    def test_anEmptyStarFileIsRefused(self):
        starFile = writeStar(self.getOutputPath('empty.star'), rows=[])
        sqliteFn = self.getOutputPath('empty.sqlite')
        pwutils.cleanPath(sqliteFn)
        partsSet = SetOfParticles(filename=sqliteFn)

        with self.assertRaises(ValueError):
            createReader(format='30').readSetOfParticles(
                starFile, partsSet, alignType=ALIGN_NONE)

    def test_thePostprocessHookSeesParticleAndRow(self):
        seen = []

        def hook(img, row):
            seen.append((img.getObjId(), row.rlnImageName))

        self._read(postprocessImageRow=hook)

        self.assertEqual([(1, '000001@Particles/parts.mrcs'),
                          (2, '000002@Particles/parts.mrcs')], seen)


class _ProtocolStub:
    """ Only what the stack linking reaches of the import protocol. """
    def __init__(self, path, copyFiles=False):
        self._path = path
        self.copyFiles = copyFiles

    def _getExtraPath(self, *paths):
        return os.path.join(self._path, *paths)

    def getCopyOrLink(self):
        return pwutils.createLink


class TestRelionImport30Stacks(BaseTest):
    """ Bringing the particle stacks of a pre-3.1 run next to the project.

    The 3.1 path symlinks the folder and prefixes it to the name in the
    star file; a 3.0 name can be relative to a parent folder, so this one
    links file by file instead.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _importer(self, copyFiles=False):
        from ..convert.dataimport import RelionImport

        srcDir = self.getOutputPath('data', 'extra')
        extraDir = self.getOutputPath('run_%s' % copyFiles, 'extra')
        pwutils.makePath(srcDir, extraDir)
        # the stack lives one folder above the star file, as in
        # relion_tutorial/import/refine3d/extra/relion_data.star
        with open(os.path.join(srcDir, '..', 'parts.mrcs'), 'w') as f:
            f.write('not really an mrcs')

        importer = RelionImport(_ProtocolStub(extraDir, copyFiles),
                                os.path.join(srcDir, 'relion_data.star'))
        importer._imgPath = srcDir

        return importer, extraDir

    def _particle(self, location=(3, '../parts.mrcs')):
        from pwem.objects import Particle

        part = Particle()
        part.setLocation(*location)

        return part

    def test_aStackAboveTheStarFileIsReachable(self):
        """ os.path.join('extra/Particles', '../parts.mrcs') would only
        work as long as nobody normalises the path. """
        importer, extraDir = self._importer()
        part = self._particle()

        importer._linkOrCopyStack(part)

        self.assertEqual(os.path.join(extraDir, 'parts.mrcs'),
                         part.getFileName())
        self.assertEqual(os.path.normpath(part.getFileName()),
                         part.getFileName())
        self.assertTrue(os.path.exists(part.getFileName()))

    def test_theStackLandsFlatInExtra(self):
        """ Whatever the folder layout of the old run, the stack is
        brought in by its base name, so the path stays normalised. """
        importer, extraDir = self._importer()
        pwutils.makePath(os.path.join(importer._imgPath, 'Particles', 'sub'))
        deep = os.path.join('Particles', 'sub', 'parts.mrcs')
        with open(os.path.join(importer._imgPath, deep), 'w') as f:
            f.write('nor this one')
        part = self._particle((1, deep))

        importer._linkOrCopyStack(part)

        self.assertEqual(os.path.join(extraDir, 'parts.mrcs'),
                         part.getFileName())
        self.assertTrue(os.path.exists(part.getFileName()))

    def test_theIndexIsKept(self):
        importer, _ = self._importer()
        part = self._particle()

        importer._linkOrCopyStack(part)

        self.assertEqual(3, part.getIndex())

    def test_linkedByDefault(self):
        importer, _ = self._importer()
        part = self._particle()

        importer._linkOrCopyStack(part)

        self.assertTrue(os.path.islink(part.getFileName()))

    def test_copiedWhenTheFormSaysSo(self):
        """ The 3.0 path used to link whatever the form said, because
        copyOrLinkFileName was called without its copyFiles argument. """
        importer, _ = self._importer(copyFiles=True)
        part = self._particle()

        importer._linkOrCopyStack(part)

        self.assertFalse(os.path.islink(part.getFileName()))
        self.assertTrue(os.path.isfile(part.getFileName()))

    def test_theSameStackIsOnlyBroughtOnce(self):
        importer, _ = self._importer()
        first, second = self._particle(), self._particle((7, '../parts.mrcs'))

        importer._linkOrCopyStack(first)
        importer._linkOrCopyStack(second)

        self.assertEqual(first.getFileName(), second.getFileName())
        self.assertEqual(7, second.getIndex())
