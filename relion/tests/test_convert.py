#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Roberto Marabini       (roberto@cnb.csic.es) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os
import subprocess
import numpy as np

from pyworkflow.tests import BaseTest, setupTestOutput, DataSet
from pyworkflow.utils import cleanPath, magentaStr, createLink, replaceExt
from pwem.objects import (SetOfParticles, CTFModel, Acquisition,
                          SetOfMicrographs, Coordinate, Particle,
                          SetOfVolumes, Transform)
from pwem.emlib.image import ImageHandler
import pwem.emlib.metadata as md
from pwem.constants import ALIGN_PROJ, ALIGN_2D, ALIGN_3D

from relion import Plugin
import relion.convert as convert
from relion.convert.convert31 import OpticsGroups
from emtable import Table


class TestConvertBinaryFiles(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.ds = DataSet.getDataSet('xmipp_tutorial')  
        cls.dsEmx = DataSet.getDataSet('emx')
        
    def test_hdfToStk(self):
        """ In this case the hdf stack files should be converted
        to .stk spider files for Relion.
        """
        print(magentaStr("\n==> Testing relion - convert hdf files to mrcs:"))
        stackFiles = ['BPV_1386_ptcls.hdf',
                      'BPV_1387_ptcls.hdf',
                      'BPV_1388_ptcls.hdf']
        
        partSet = SetOfParticles(filename=':memory:')
        
        for fn in stackFiles:
            particle = Particle()
            particle.setLocation(1, self.ds.getFile('particles/%s' % fn))
            partSet.append(particle)
            
        outputDir = self.getOutputPath()
        filesDict = convert.convertBinaryFiles(partSet, outputDir)
        partSet.close()
        print(filesDict)
        
    def test_mrcsLink(self):
        """ In this case just a link with .mrcs extension 
        should be created """
        print(magentaStr("\n==> Testing relion - link mrc stack to mrcs:"))
        stackFile = self.dsEmx.getFile('particles/particles.mrc')
        partSet = SetOfParticles(filename=':memory:')
        
        for i in range(1, 10):
            particle = Particle()
            particle.setLocation(i, stackFile)
            partSet.append(particle)
            
        outputDir = self.getOutputPath()
        filesDict = convert.convertBinaryFiles(partSet, outputDir)
        print(filesDict)


SHOW_IMAGES = False  # Launch xmipp_showj to open intermediate results
CLEAN_IMAGES = True  # Remove the output temporary files
PRINT_MATRIX = True
PRINT_FILES = True


def runRelionProgram(cmd):
    print(">>>", cmd)
    cmd = cmd.split()
    p = subprocess.Popen(cmd, env=Plugin.getEnviron())
    return p.wait()


class TestConvertAnglesBase(BaseTest):
    """ Base class to launch both Alignment and Reconstruction tests."""
    IS_ALIGNMENT = None
    CMD = None

    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.dataset = DataSet.getDataSet('emx')

    def launchTest(self, fileKey, mList, alignType=None, **kwargs):
        """ Helper function to launch similar alignment tests
        given the EMX transformation matrix.
        Params:
            fileKey: the file where to grab the input stack images.
            mList: the matrix list of transformations
                (should be the same length of the stack of images)
        """
        print("\n")
        print("*" * 80)
        print("* Launching test: ", fileKey)
        print("*" * 80)

        is2D = alignType == ALIGN_2D

        if fileKey == 'alignShiftRotExp':
            # relion requires mrcs stacks
            origFn = self.dataset.getFile(fileKey)
            stackFn = replaceExt(origFn, ".mrcs")
            createLink(origFn, stackFn)
        else:
            stackFn = self.dataset.getFile(fileKey)

        partFn1 = self.getOutputPath(fileKey + "_particles1.sqlite")
        mdFn = self.getOutputPath(fileKey + "_particles.star")
        partFn2 = self.getOutputPath(fileKey + "_particles2.sqlite")

        if self.IS_ALIGNMENT:
            outputFn = self.getOutputPath(fileKey + "_output.mrcs")
            outputFnRelion = self.getOutputPath(fileKey + "_output")
            goldFn = self.dataset.getFile(fileKey + '_Gold_output_relion.mrcs')
        else:
            outputFn = self.getOutputPath(fileKey + "_output.vol")
            goldFn = self.dataset.getFile("reconstruction/gold/" + fileKey + '_Gold_rln_output.vol')

        if PRINT_FILES:
            print("BINARY DATA: ", stackFn)
            print("SET1:        ", partFn1)
            print("  MD:        ", mdFn)
            print("SET2:        ", partFn2)
            print("OUTPUT:      ", outputFn)
            print("GOLD:        ", goldFn)

        if alignType == ALIGN_2D or alignType == ALIGN_PROJ:
            partSet = SetOfParticles(filename=partFn1)
        else:
            partSet = SetOfVolumes(filename=partFn1)
        partSet.setAlignment(alignType)

        acq = Acquisition(voltage=300,
                          sphericalAberration=2,
                          amplitudeContrast=0.1,
                          magnification=60000)
        og = OpticsGroups.create(rlnMtfFileName="mtfFile1.star",
                                 rlnImageSize=128)
        partSet.setSamplingRate(1.0)
        partSet.setAcquisition(acq)
        og.toImages(partSet)
        # Populate the SetOfParticles with images
        # taken from images.mrc file
        # and setting the previous alignment parameters
        aList = [np.array(m) for m in mList]
        for i, a in enumerate(aList):
            p = Particle()
            p.setLocation(i + 1, stackFn)
            p.setTransform(Transform(a))
            partSet.append(p)
        # Write out the .sqlite file and check that are correctly aligned
        print("Partset", partFn1)
        partSet.printAll()
        partSet.write()
        # Convert to a Xmipp metadata and also check that the images are
        # aligned correctly
        if alignType == ALIGN_2D or alignType == ALIGN_PROJ:
            starWriter = convert.createWriter()
            starWriter.writeSetOfParticles(partSet, mdFn, alignType=alignType)
            partSet2 = SetOfParticles(filename=partFn2)
        else:
            convert.writeSetOfVolumes(partSet, mdFn, alignType=alignType)
            partSet2 = SetOfVolumes(filename=partFn2)
        # Let's create now another SetOfImages reading back the written
        # Xmipp metadata and check one more time.
        partSet2.copyInfo(partSet)
        if alignType == ALIGN_2D or alignType == ALIGN_PROJ:
            convert.readSetOfParticles(mdFn, partSet2,
                                       alignType=alignType)
        else:
            convert.readSetOfParticles(mdFn, partSet2,
                                       rowToFunc=convert.rowToVolume,
                                       alignType=alignType)

        partSet2.write()

        if PRINT_MATRIX:
            for i, img in enumerate(partSet2):
                m1 = aList[i]
                m2 = img.getTransform().getMatrix()
                print("-" * 5)
                print(img.getFileName(), img.getIndex())
                print('m1:\n', m1, convert.geometryFromMatrix(m1, False))

                print('m2:\n', m2, convert.geometryFromMatrix(m2, False))
                self.assertTrue(np.allclose(m1, m2, rtol=1e-2))

        # Launch apply transformation and check result images
        runRelionProgram(self.CMD % locals())

        if SHOW_IMAGES:
            runRelionProgram('scipion show %(outputFn)s' % locals())

        if os.path.exists(goldFn):
            self.assertTrue(
                ImageHandler().compareData(goldFn, outputFn, tolerance=0.001),
                "Different data files:\n>%s\n<%s" % (goldFn, outputFn))

        if CLEAN_IMAGES:
            cleanPath(outputFn)


class TestAlignment(TestConvertAnglesBase):
    IS_ALIGNMENT = True
    CMD = "relion_stack_create  --i %(mdFn)s --o %(outputFnRelion)s --apply_transformation"

    def test_isInverse(self):
        """Consistency between forwards and backwards geometrical
        transformations"""
        def _testInv(matrix):
            matrix = np.array(matrix)
            a = Transform(matrix)

            row1 = md.Row()
            convert.alignmentToRow(a, row1, alignType=ALIGN_2D)

            row2 = None
            
            row3 = md.Row()
            convert.alignmentToRow(a, row3, alignType=ALIGN_PROJ)

            return row1, row2, row3

        row1, row2, row3 = _testInv([[1.0, 0.0, 0.0, 20.0],
                                    [0.0, 1.0, 0.0, 0.0],
                                    [0.0, 0.0, 1.0, 0.0],
                                    [0.0, 0.0, 0.0, 1.0]])

        self.assertAlmostEqual(row1.getValue(md.RLN_ORIENT_ORIGIN_X), 20., 4)
        self.assertAlmostEqual(row1.getValue(md.RLN_ORIENT_ORIGIN_Y), 0., 4)
        self.assertAlmostEqual(row1.getValue(md.RLN_ORIENT_PSI), 0., 4)

        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ORIGIN_X), 20., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ORIGIN_Y), 0., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ORIGIN_Z), 0., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ROT), 0., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_TILT), 0., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_PSI), 0., 4)

        row1, row2, row3 = _testInv([[0.93969, 0.34202, 0.0, -6.8404],
                                    [-0.34202, 0.93969, 0.0, 18.7939],
                                    [0.0, 0.0, 1.0, 0.0],
                                    [0.0, 0.0, 0.0, 1.0]])
        self.assertAlmostEqual(row1.getValue(md.RLN_ORIENT_ORIGIN_X), -6.8404, 4)
        self.assertAlmostEqual(row1.getValue(md.RLN_ORIENT_ORIGIN_Y), 18.7939, 4)
        self.assertAlmostEqual(row1.getValue(md.RLN_ORIENT_PSI), -20., 4)

        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ORIGIN_X), -12.8558097352, 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ORIGIN_Y), 15.3209632479, 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ORIGIN_Z), 0., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_ROT), -20., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_TILT), 0., 4)
        self.assertAlmostEqual(row3.getValue(md.RLN_ORIENT_PSI), 0., 4)

    def test_alignShiftRotExp(self):
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: 2D alignment
        Misalignment: angles, shifts
        """
        mList = [[[1.0, 0.0, 0.0, 20.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[0.86602539, 0.5, 0.0, 20.0],
                  [-0.5, 0.86602539, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[0.86602539, 0.5, 0.0, 27.706396],
                  [-0.5, 0.86602539, 0.0, 0.331312],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[0.5, 0.86602539, 0.0, 3.239186],
                  [-0.86602539, 0.5, 0.0, 20.310715],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[0.0, 1.0, 0.0, 10.010269],
                  [-1.0, 0.0, 0.0, 3.6349521],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]]

        self.launchTest('alignShiftRotExp', mList, alignType=ALIGN_2D)

    def aatest_alignShiftRot3D(self):
        # TODO: 3D alignment not tested since we need a program in Relion
        # apply transformation.
        
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: 3D alignment
        Misalignment: angles, shifts
        0.71461016 0.63371837 -0.29619813         15
        -0.61309201 0.77128059 0.17101008         25
        0.33682409 0.059391174 0.93969262         35
                 0          0          0          1
        """
        mList = [[[0.71461016, 0.63371837, -0.29619813, 15],
                  [-0.61309201, 0.77128059, 0.17101008, 25],
                  [0.33682409, 0.059391174, 0.93969262, 35],
                  [0, 0, 0, 1]],

                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]]

        self.launchTest('alignShiftRot3D', mList, alignType=ALIGN_3D)

    def test_alignRotOnly(self):
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: 2D alignment
        Misalignment: angles
        mList[0] * (0,32,0,1)' -> (32,0,0,1)'  T -> ref
        mList[1] * (-32,0,0,1)' -> (32,0,0,1)'
        """
        mList = [[[0.0, 1.0, 0.0, 0.0],
                  [-1.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[-1.0, 0.0, 0.0, 0.0],
                  [0.0, -1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]]

        self.launchTest('alignRotOnly', mList, alignType=ALIGN_2D)

    def aatest_alignRotOnly3D(self):
        # TODO: 3D alignment not tested since we need a program in Relion
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: 3D alignment
        Misalignment: angles
           mList[0] * (22.8586, -19.6208, 10.7673)' -> (32,0,0,1); T -> ref
           mList[1] * (22.8800, 20.2880, -9.4720) ->(32,0,0,1)'; T -> ref
        """

        mList = [[[0.715, -0.613, 0.337, 0.],
                  [0.634, 0.771, 0.059, 0.],
                  [-0.296, 0.171, 0.94, 0.],
                  [0., 0., 0., 1.]],
                 [[0.71433, 0.63356, -0.29586, -0.00000],
                  [-0.61315, 0.77151, 0.17140, -0.00000],
                  [0.33648, 0.05916, 0.93949, -0.00000],
                  [0.00000, 0.00000, 0.00000, 1.00000]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]]

        self.launchTest('alignRotOnly3D', mList, alignType=ALIGN_3D)

    def test_alignShiftOnly(self):
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: 2D alignment
        Misalignment: shifts
         mList[0] * ( 32,0,0,1)-> (0,0,0,1); T -> ref
         mList[1] * (  0,32,0,1)-> (0,0,0,1); T -> ref

        """
        mList = [[[1.0, 0.0, 0.0, -32.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, -32.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]]

        self.launchTest('alignShiftOnly', mList, alignType=ALIGN_2D)

    def aatest_alignShiftOnly3D(self):
        # TODO: 3D alignment not tested since we need a program in Relion
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: 3D alignment
        Misalignment: shifts
        mList[0] * (0,0,0) -> (32,0,0)
        mList[1] * (32,-16,0) -> (32,0,0)
        mList[1] * (32,0,-8) -> (32,0,0)
        mList[2] reference
        """
        mList = [[[1.0, 0.0, 0.0, 32.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 16.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 8.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]]

        self.launchTest('alignShiftOnly3D', mList, alignType=ALIGN_3D)

    def test_alignShiftRot(self):
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: 2D alignment
        Misalignment: angles and shifts
        note that when rot and shift are involved
        the correct transformation matrix is
        rot matrix + rot*shift
        mList[0] * (16,32,0) -> (32,0,0)
        mList[1] * (-32,16,0) -> (32,0,0)
        mList[2] reference
        """
        mList = [[[0.0, 1.0, 0.0, 0.0],
                  [-1.0, 0.0, 0.0, 16.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[-1.0, 0.0, 0.0, 0.0],
                  [0.0, -1.0, 0.0, 16.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[1.0, 0.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]]]

        self.launchTest('alignShiftRot', mList, alignType=ALIGN_2D)


class TestReconstruct(TestConvertAnglesBase):
    IS_ALIGNMENT = False
    CMD = "relion_reconstruct --i %(mdFn)s --o %(outputFn)s"

    def test_forward_backwards(self):
        """convert transformation matrix to xmipp and back"""

        mList = [[[0.71461016, 0.63371837, -0.29619813,  1.],  # a1
                  [-0.61309201, 0.77128059, 0.17101008,  2.],
                  [0.33682409, 0.059391174, 0.93969262,  3.],
                  [0,          0,          0,            1.]],
                 [[0., 0., -1., 0.],  # a2
                  [0., 1., 0., 0.],
                  [1., 0., 0., 0.],
                  [0., 0., 0., 1.]],
                 [[0., 1., 0., 0.],  # a3
                  [0., 0., 1., 0.],
                  [1., 0., 0., 0.],
                  [0., 0., 0., 1.]],
                 [[0.22612257, 0.82379508, -0.51983678, 0.],  # a4
                  [-0.88564873, 0.39606407, 0.24240388,  0.],
                  [0.40557978, 0.40557978, 0.81915206,  0.],
                  [0.,          0.,          0.,           1.]],
                 [[-0.78850311, -0.24329656, -0.56486255,   0.],  # a5
                  [0.22753462, -0.96866286, 0.099600501,  0.],
                  [-0.57139379, -0.049990479, 0.81915206,  0.],
                  [0.,            0.,           0.,           1.]],
                 [[1.0, 0.0, 0.0, 0.0],  # a6
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[0., 0., -1., 0.],  # a7
                  [-1., 0., 0.,  0.],
                  [0., 1., 0.,  0.],
                  [0., 0., 0., 1.]]
                 ]

        aList = [np.array(m) for m in mList]
        rowa = md.Row()
        rowb = md.Row()
        labelList = [md.RLN_ORIENT_ROT,
                     md.RLN_ORIENT_TILT,
                     md.RLN_ORIENT_PSI,
                     md.RLN_ORIENT_ORIGIN_X,
                     md.RLN_ORIENT_ORIGIN_Y,
                     md.RLN_ORIENT_ORIGIN_Z]

        for i, a in enumerate(aList):
            a = Transform(aList[i])
            convert.alignmentToRow(a, rowa, ALIGN_PROJ)
            b = convert.rowToAlignment(rowa, ALIGN_PROJ)
            convert.alignmentToRow(b, rowb, ALIGN_PROJ)
            # same two matrices
            self.assertTrue(np.allclose(a.getMatrix(), b.getMatrix(),
                                        rtol=1e-2))
            for label in labelList:
                auxBtilt = rowb.getValue(label)
                auxAtilt = rowa.getValue(label)
                # same two rows
                self.assertAlmostEqual(auxBtilt, auxAtilt, places=3)

            b = convert.rowToAlignment(rowa, ALIGN_PROJ)
            convert.alignmentToRow(b, rowb, ALIGN_PROJ)
            aMatrix = a.getMatrix()
            # aMatrix[0,:] *= -1; aMatrix[2,:] *= -1;
            # same two matrices with flip
            print("\naMatrix: \n", aMatrix, "\nbMatrix: \n", b.getMatrix())
            
            self.assertTrue(np.allclose(aMatrix, b.getMatrix(), rtol=1e-2))

    def test_reconstRotOnly(self):
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: reconstruction from projections
        Misalignment: angles
        Volume to be reconstructed is empty except for an small sphere
        at (4,8,16)
        mList[0] * (3.3429    9.6554   15.2184)   -> (4,8,16)
        mList[1] * (14    4    4)   -> (4,8,16)
        mList[2] * (16    4    8)   -> (4,8,16)
        mList[3] * ( 0.30858   12.95297   12.96632)   -> (4,8,16)
        mList[4] * (-10.4760   -9.5223   11.6438)   -> (4,8,16)
        mList[5] * (4    8   16)   -> (4,8,16)
        mList[6] * (-8   16   -4)  -> (4,8,16)
        """
        mList = [[[0.71461016, 0.63371837, -0.29619813,  0.],  # a1
                  [-0.61309201, 0.77128059, 0.17101008,  0.],
                  [0.33682409, 0.059391174, 0.93969262,  0.],
                  [0,          0,          0,            1.]],
                 [[0., 0., -1., 0.],  # a2
                  [0., 1., 0., 0.],
                  [1., 0., 0., 0.],
                  [0., 0., 0., 1.]],
                 [[0., 1., 0., 0.],  # a3
                  [0., 0., 1., 0.],
                  [1., 0., 0., 0.],
                  [0., 0., 0., 1.]],
                 [[0.22612257, 0.82379508, -0.51983678, 0.],  # a4
                  [-0.88564873, 0.39606407, 0.24240388,  0.],
                  [0.40557978, 0.40557978, 0.81915206,  0.],
                  [0.,          0.,          0.,           1.]],
                 [[-0.78850311, -0.24329656, -0.56486255,   0.],  # a5
                  [0.22753462, -0.96866286, 0.099600501,  0.],
                  [-0.57139379, -0.049990479, 0.81915206,  0.],
                  [0.,            0.,           0.,           1.]],
                 [[1.0, 0.0, 0.0, 0.0],  # a6
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[0., 0., -1., 0.],  # a7
                  [-1., 0., 0.,  0.],
                  [0., 1., 0.,  0.],
                  [0., 0., 0., 1.]]
                 ]

        self.launchTest('reconstRotOnly', mList, alignType=ALIGN_PROJ)

    def test_reconstRotandShift(self):
        """ Check that for a given alignment object,
        the corresponding Xmipp metadata row is generated properly.
        Goal: reconstruction from projections
        Misalignment: angles and shifts
        Volume to be reconstructed is empty except for an small sphere
        at (4,8,16)
        mList[0] * 64*(0.0210515   0.0037119   0.0587308)   -> (4,8,16)
        mList[1] * 64*(0.250000   0.125000   0.062500)   -> (4,8,16)
        mList[2] * 64*(0.218750   0.062500   0.062500)   -> (4,8,16)
        mList[3] * 64*(-0.0037446   0.0683733   0.1968962)   -> (4,8,16)
        mList[4] * 64*(-0.109277  -0.176747   0.256331)   -> (4,8,16)
        mList[5] * 64*(4    8   16)   -> (4,8,16)
        mList[6] * 64*(-0.125000   0.250000  -0.062500)  -> (4,8,16)
        """
        mList = [[[0.71461016, 0.63371837, -0.29619813,  4.],  # a1
                  [-0.61309201, 0.77128059, 0.17101008,  8.],
                  [0.33682409, 0.059391174, 0.93969262, 12.],
                  [0,          0,          0,            1.]],
                 [[0., 0., -1., 0.],  # a2
                  [0., 1., 0., 0.],
                  [1., 0., 0., 0.],
                  [0., 0., 0., 1.]],
                 [[0., 1., 0., 0.],  # a3
                  [0., 0., 1., 4.],
                  [1., 0., 0., 2.],
                  [0., 0., 0., 1.]],
                 [[0.22612257, 0.82379508, -0.51983678, 7.],  # a4
                  [-0.88564873, 0.39606407, 0.24240388,  3.],
                  [0.40557978, 0.40557978, 0.81915206,  4.],
                  [0.,          0.,          0.,           1.]],
                 [[-0.78850311, -0.24329656, -0.56486255,   5.],  # a5
                  [0.22753462, -0.96866286, 0.099600501, -3.],
                  [-0.57139379, -0.049990479, 0.81915206, -2.],
                  [0.,            0.,           0.,           1.]],
                 [[1.0, 0.0, 0.0, 0.0],  # a6
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0]],
                 [[0., 0., -1., 0.],  # a7
                  [-1., 0., 0.,  0.],
                  [0., 1., 0.,  0.],
                  [0., 0., 0., 1.]]
                 ]

        self.launchTest('reconstRotandShift', mList, alignType=ALIGN_PROJ)


class TestRelionWriter(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def _createSetOfMics(self, n=10, nOptics=2):
        micName = 'BPV_13%02d.mrc'
        psdName = 'BPV_13%02d_PSD.ctf:mrc'
        ogName = 'opticsGroup%d'
        mtfFile = 'mtfFile%d.star'

        cleanPath(self.getOutputPath('micrographs.sqlite'))
        micsDb = self.getOutputPath('micrographs.sqlite')
        outputMics = SetOfMicrographs(filename=micsDb)
        outputMics.setSamplingRate(1.234)

        mic = SetOfMicrographs.ITEM_TYPE()
        acq = Acquisition(voltage=300,
                          sphericalAberration=2,
                          amplitudeContrast=0.1,
                          magnification=60000)

        og = OpticsGroups.create(rlnMtfFileName='')

        fog = og.first()

        ctf = CTFModel(defocusU=10000, defocusV=15000, defocusAngle=15)
        outputMics.setAcquisition(acq)
        mic.setAcquisition(acq)
        mic.setCTF(ctf)

        itemsPerOptics = n // nOptics

        for i in range(1, n+1):
            mic.setFileName(micName % i)
            ctf = mic.getCTF()
            ctf.setPsdFile(psdName % i)
            ctf.setFitQuality(np.random.uniform())
            ctf.setResolution(np.random.uniform(3, 15))
            ogNumber = (i-1) // itemsPerOptics + 1

            ogDict = {
                'rlnOpticsGroup': ogNumber,
                'rlnOpticsGroupName': ogName % ogNumber,
                'rlnMtfFileName': mtfFile % ogNumber
            }

            if ogNumber in og:
                og.update(ogNumber, **ogDict)
            else:
                og.add(fog._replace(**ogDict))

            mic.rlnOpticsGroup = ogNumber
            mic.setObjId(None)
            outputMics.append(mic)

        print(">>> Writing micrograph set to: ", micsDb)
        outputMics.write()

        return outputMics

    def _createSetOfParts(self, nMics=10, nOptics=2, partsPerMic=10):
        micSet = self._createSetOfMics(nMics, nOptics)
        outputSqlite = self.getOutputPath('particles.sqlite')
        cleanPath(outputSqlite)
        print(">>> Writing to particles db: %s" % outputSqlite)
        outputParts = SetOfParticles(filename=outputSqlite)
        outputParts.setSamplingRate(1.234)
        outputParts.setAcquisition(micSet.getAcquisition())

        part = SetOfParticles.ITEM_TYPE()
        coord = Coordinate()

        for mic in micSet:
            for i in range(1, partsPerMic+1):
                part.setLocation(i, mic.getFileName().replace('mrc', 'mrcs'))
                coord.setPosition(x=np.random.randint(0, 1000),
                                  y=np.random.randint(0, 1000))
                coord.setMicrograph(mic)
                part.setObjId(None)
                part.setCoordinate(coord)
                part.setAcquisition(mic.getAcquisition())
                outputParts.append(part)

        outputParts.write()

        return outputParts

    def test_micrographsToStar(self):
        """ Write a SetOfParticles to Relion star input file. """
        outputStar = self.getOutputPath("micrographs.star")
        outputMics = self._createSetOfMics(10)
        print(">>> Writing to micrographs: %s" % outputStar)
        starWriter = convert.createWriter()
        starWriter.writeSetOfMicrographs(outputMics, outputStar)

    def test_particlesToStar(self):
        """ Write a SetOfParticles to Relion star input file. """
        outputStar = self.getOutputPath("particles.star")
        outputParts = self._createSetOfParts(10, 2, 10)
        print(">>> Writing to particles star: %s" % outputStar)
        starWriter = convert.createWriter()
        starWriter.writeSetOfParticles(outputParts, outputStar)

    def test_particlesImportToStar(self):
        sqliteFn = self.ds.getFile("import/case2/particles.sqlite")
        partsSet = SetOfParticles(filename=sqliteFn)
        partsSet.loadAllProperties()
        outputStar = self.getOutputPath("particles.star")
        print(">>> Writing to particles star: %s" % outputStar)
        starWriter = convert.createWriter()
        starWriter.writeSetOfParticles(partsSet, outputStar)


class TestRelionReader(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.ds = DataSet.getDataSet('relion31_tutorial_precalculated')
        cls.ds2 = DataSet.getDataSet('relion_tutorial')
        cls.starOldFormat = cls.ds2.getFile('import/refine3d_case2/relion_data.star')

    def __readParticles(self, partsStar, outputSqlite=None, **kwargs):
        outputSqlite = outputSqlite or self.getOutputPath('particles.sqlite')
        print("<<< Reading star file: \n   %s\n" % partsStar)
        cleanPath(outputSqlite)
        print(">>> Writing to particles db: \n   %s\n" % outputSqlite)
        partsSet = SetOfParticles(filename=outputSqlite)
        convert.readSetOfParticles(partsStar, partsSet, **kwargs)
        return partsSet

    def test_readSetOfParticles(self):
        partsSet = self.__readParticles(
            self.ds.getFile("Extract/job018/particles.star"),
            extraLabels=['rlnNrOfSignificantSamples']
        )
        partsSet.write()

        first = partsSet.getFirstItem()
        first.printAll()

        self.assertAlmostEqual(first.getSamplingRate(), 1.244531)
        self.assertEqual(first.getClassId(), 4)
        self.assertTrue(hasattr(first, '_rlnNrOfSignificantSamples'))

        fog = OpticsGroups.fromImages(partsSet).first()
        self.assertEqual(fog.rlnMtfFileName, 'mtf_k2_200kV.star')
        self.assertEqual(fog.rlnOpticsGroupName, 'opticsGroup1')

    def test_readSetOfParticlesAfterCtf(self):
        starFile = self.ds.getFile("CtfRefine/job023/particles_ctf_refine.star")
        partsReader = Table.Reader(starFile, tableName='particles',
                                   types=convert.LABELS_DICT)
        firstRow = partsReader.getRow()

        partsSet = self.__readParticles(starFile)
        first = partsSet.getFirstItem()

        ogLabels = ['rlnBeamTiltX', 'rlnBeamTiltY']
        extraLabels = ['rlnCtfBfactor', 'rlnCtfScalefactor', 'rlnPhaseShift']
        for l in extraLabels:
            value = getattr(first, '_%s' % l)
            self.assertIsNotNone(value, "Missing label: %s" % l)
            self.assertAlmostEqual(getattr(firstRow, l), value)

        fog = OpticsGroups.fromImages(partsSet).first()
        self.assertTrue(all(hasattr(fog, l) for l in ogLabels))

        # Also test writing and preserving extra labels
        outputStar = self.getOutputPath('particles.star')
        print(">>> Writing to particles star: %s" % outputStar)
        starWriter = convert.createWriter()
        starWriter.writeSetOfParticles(partsSet, outputStar)

        fog = OpticsGroups.fromStar(outputStar).first()
        self.assertTrue(all(hasattr(fog, l) for l in ogLabels))

        partsReader = Table.Reader(outputStar, tableName='particles',
                                   types=convert.LABELS_DICT)
        firstRow = partsReader.getRow()
        for l in extraLabels:
            value = getattr(first, '_%s' % l)
            self.assertIsNotNone(value, "Missing label: %s" % l)
            self.assertAlmostEqual(getattr(firstRow, l), value)

    def test_readOldStarFormat(self):
        partsSet = self.__readParticles(
            self.starOldFormat,
            alignType=ALIGN_PROJ,
            format='30')
        partsSet.write()
        first = partsSet.getFirstItem()
        first.printAll()

        coord = first.getCoordinate()
        x, y = coord.getPosition()
        self.assertEqual(first.getClassId(), 1)
        self.assertEqual(x, 299)
        self.assertEqual(coord.getMicName(), 'Falcon_2012_06_12-14_33_35_0_movie.mrcs')


class TestRelionOpticsGroups(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls.ds = DataSet.getDataSet('relion31_tutorial_precalculated')

    def test_fromStar(self):
        partsStar = self.ds.getFile("Extract/job018/particles.star")

        print("<<< Reading optics groups from file: \n   %s\n" % partsStar)
        og = OpticsGroups.fromStar(partsStar)
        fog = og.first()

        # test hasColumn method
        for colName in ['rlnMtfFileName', 'rlnOpticsGroupName']:
            self.assertTrue(og.hasColumn(colName))

        self.assertEqual(fog.rlnMtfFileName, 'mtf_k2_200kV.star')
        self.assertEqual(fog.rlnOpticsGroupName, 'opticsGroup1')
        self.assertEqual(og['opticsGroup1'], fog)

    def test_string(self):
        og = OpticsGroups.create(rlnMtfFileName='mtf_k2_200kV.star')
        fog = og.first()

        self.assertEqual(fog.rlnMtfFileName, 'mtf_k2_200kV.star')
        self.assertEqual(fog.rlnOpticsGroupName, 'opticsGroup1')
        self.assertEqual(og['opticsGroup1'], fog)

        # try update by id
        og.update(1, rlnMtfFileName="new_mtf_k2.star")
        # try update by name
        og.update('opticsGroup1', rlnImageSize=512)

        fog = og.first()
        self.assertEqual(fog.rlnMtfFileName, 'new_mtf_k2.star')
        self.assertEqual(fog.rlnImageSize, 512)
        self.assertEqual(fog.rlnOpticsGroupName, 'opticsGroup1')
        self.assertEqual(og['opticsGroup1'], fog)

    def test_add(self):
        """ Testing adding more groups or columns. """
        og = OpticsGroups.create()

        og1 = og.first()

        self.assertEqual(len(og), 1)
        self.assertAlmostEqual(og1.rlnVoltage, 300.)

        og2 = og1._replace(rlnOpticsGroup=2,
                           rlnOpticsGroupName='opticsGroup2',
                           rlnVoltage=200.)
        og.add(og2)
        og3 = og1._replace(rlnOpticsGroup=3,
                           rlnOpticsGroupName='opticsGroup3',
                           rlnVoltage=100.)
        og.add(og3)

        self.assertEqual(len(og), 3)

        og.addColumns(rlnMtfFileName='mtf_k2_200kV.star')

        self.assertTrue(all(hasattr(ogx, 'rlnMtfFileName') for ogx in og))

        og.update(2, rlnVoltage=300.)
        og.update(3, rlnVoltage=300.)

        for ogx in og:
            self.assertAlmostEqual(ogx.rlnVoltage, 300.)

        og.updateAll(rlnVoltage=200.)

        for ogx in og:
            self.assertAlmostEqual(ogx.rlnVoltage, 200.)
