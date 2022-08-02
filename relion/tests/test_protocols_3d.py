# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es) [1]
# *             Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *             J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [3]
# *
# * [1] Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
# * [3] SciLifeLab, Stockholm University
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

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.plugin import Domain
from pyworkflow.protocol.constants import STATUS_FINISHED
import pyworkflow.utils as pwutils
from pwem.objects import SetOfParticles, Volume
from pwem.protocols import ProtImportParticles, ProtImportPdb
from pwem.emlib.image import ImageHandler

from ..protocols import *
from .test_protocols_base import TestRelionBase, USE_GPU, RUN_CPU


class TestRelionCalculateFSC(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.half1 = cls.ds.getFile('volumes/emd_3488_Noisy_half1.vol')
        cls.half2 = cls.ds.getFile('volumes/emd_3488_Noisy_half2.vol')
        cls.map = cls.ds.getFile('volumes/emd_3488.map')
        cls.model = cls.ds.getFile('PDBx_mmCIF/5ni1.pdb')

        cls.protImport1 = cls.runImportVolumes(cls.half1, 1.05)
        cls.protImport2 = cls.runImportVolumes(cls.half2, 1.05)
        cls.protImportMap = cls.runImportVolumes(cls.map, 1.05)

        print(pwutils.magentaStr("\n==> Importing data - pdb:"))
        cls.protImportModel = cls.newProtocol(ProtImportPdb,
                                              inputPdbData=1,
                                              pdbFile=cls.model)
        cls.launchProtocol(cls.protImportModel)

    def test_fsc_overall(self):
        print(pwutils.magentaStr("\n==> Testing relion - calculate fsc:"))
        calcRun = self.newProtocol(
            ProtRelionCalculateFSC,
            objLabel='FSC overall',
            fscType=0,
            half1=self.protImport1.outputVolume,
            half2=self.protImport2.outputVolume
        )
        self.launchProtocol(calcRun)
        self.assertIsNotNone(calcRun.outputFSC)

    def test_fsc_model_map(self):
        print(pwutils.magentaStr("\n==> Testing relion - calculate fsc:"))
        calcRun = self.newProtocol(
            ProtRelionCalculateFSC,
            objLabel='FSC model-map',
            fscType=1,
            model=self.protImportModel.outputPdb,
            map=self.protImportMap.outputVolume
        )
        self.launchProtocol(calcRun)
        self.assertIsNotNone(calcRun.outputFSC)

    def test_fsc_work_free(self):
        print(pwutils.magentaStr("\n==> Testing relion - calculate fsc:"))
        calcRun = self.newProtocol(
            ProtRelionCalculateFSC,
            objLabel='FSC work/free',
            fscType=2,
            model_half1=self.protImportModel.outputPdb,
            half1=self.protImport1.outputVolume,
            half2=self.protImport2.outputVolume
        )
        self.launchProtocol(calcRun)
        self.assertIsNotNone(calcRun.outputSetOfFSCs)


class TestRelionClassify3D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testProtRelionClassify3D(self):
        def _runRelionClassify3D(doGpu=False, label=''):
            relion3DClass = self.newProtocol(ProtRelionClassify3D,
                                             numberOfClasses=3,
                                             numberOfIterations=4,
                                             doCTF=False, runMode=1,
                                             maskDiameterA=320,
                                             numberOfMpi=2, numberOfThreads=2)

            relion3DClass.setObjLabel(label)
            relion3DClass.inputParticles.set(self.protNormalize.outputParticles)
            relion3DClass.referenceVolume.set(self.protImportVol.outputVolume)
            relion3DClass.doGpu.set(doGpu)
            self.launchProtocol(relion3DClass)
            return relion3DClass

        def _checkAsserts(relionProt):
            self.assertIsNotNone(relionProt.outputClasses, "There was a "
                                                           "problem with "
                                                           "Relion 3D classify")

            for class3D in relionProt.outputClasses:
                self.assertTrue(class3D.hasAlignmentProj())

        if RUN_CPU:
            print(pwutils.magentaStr("\n==> Testing relion - classify 3d on CPU:"))
            relionProt = _runRelionClassify3D(doGpu=False,
                                              label="Run Relion classify3D CPU")
            _checkAsserts(relionProt)

        if USE_GPU:
            print(pwutils.magentaStr("\n==> Testing relion - classify 3d on GPU:"))
            relionGpu = _runRelionClassify3D(doGpu=True,
                                             label="Relion classify3D GPU")
            _checkAsserts(relionGpu)


class TestRelionCreate3dMask(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.volFn = cls.ds.getFile('import/refine3d/extra/relion_class001.mrc')
        cls.protImportVol = cls.runImportVolumes(cls.volFn, 3.0)

    def _validations(self, mask, dims, pxSize, prot):
        self.assertIsNotNone(mask, "There was a problem with mask 3d protocol, "
                                   "using %s protocol as input" % prot)
        xDim = mask.getXDim()
        sr = mask.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your volume is (%d)^3 "
                                     "and must be (%d)^3" % (xDim, dims))

        self.assertAlmostEqual(sr, pxSize, delta=0.0001,
                               msg="Pixel size of your volume is %0.5f and"
                                   " must be %0.5f" % (sr, pxSize))

    def test_createMask(self):
        maskProt = self.newProtocol(ProtRelionCreateMask3D,
                                    initialLowPassFilterA=10)  # filter at 10 A
        vol = self.protImportVol.outputVolume
        maskProt.inputVolume.set(vol)
        print(pwutils.magentaStr("\n==> Testing relion - create mask 3d:"))
        self.launchProtocol(maskProt)

        self._validations(maskProt.outputMask, vol.getXDim(),
                          vol.getSamplingRate(), maskProt)

        ih = ImageHandler()
        img = ih.read(maskProt.outputMask)
        mean, std, _min, _max = img.computeStats()
        # Check the mask is non empty and between 0 and 1
        self.assertAlmostEqual(_min, 0)
        self.assertAlmostEqual(_max, 1)
        self.assertTrue(mean > 0)


class TestRelionExpandSymmetry(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partRef3dFn = cls.dataset.getFile('import/refine3d/extra/relion_data.star')
        cls.protImport = cls.runImportParticlesStar(cls.partRef3dFn, 7.08)

    def test_ExpandSymmetry(self):
        print(pwutils.magentaStr("\n==> Testing relion - expand symmetry:"))
        prot = self.newProtocol(ProtRelionExpandSymmetry)
        prot.inputParticles.set(self.protImport.outputParticles)
        prot.symmetryGroup.set("D2")
        self.launchProtocol(prot)

        self.assertIsNotNone(prot.outputParticles,
                             "There was a problem with expand symmetry protocol")
        sizeIn = self.protImport.outputParticles.getSize()
        sizeOut = prot.outputParticles.getSize()
        self.assertAlmostEqual(sizeIn * 4, sizeOut, delta=0.0001,
                               msg="Number of output particles is %d and"
                                   " must be %d" % (sizeOut, sizeIn * 4))


class TestRelionInitialModel(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.dataset.getFile('import/classify2d/extra/relion_it015_data.star')
        cls.protImport = cls.runImportParticlesStar(cls.partFn, 7.08)

    def testProtRelionIniModel(self):
        def _runRelionIniModel(doGpu=True, label=''):
            print(pwutils.magentaStr("\n==> Testing relion - initial model:"))
            relionIniModel = self.newProtocol(ProtRelionInitialModel,
                                              doCTF=False,
                                              doGpu=doGpu,
                                              maskDiameterA=340,
                                              symmetryGroup='C1',
                                              allParticlesRam=True,
                                              numberOfMpi=1,
                                              numberOfThreads=8)
            relionIniModel.setObjLabel(label)
            relionIniModel.inputParticles.set(self.protImport.outputParticles)
            self.launchProtocol(relionIniModel)

            return relionIniModel

        def _checkAsserts(relionProt):
            relionProt._initialize()  # Load filename templates
            dataSqlite = relionProt._getIterData(relionProt._lastIter())
            outImgSet = SetOfParticles(filename=dataSqlite)

            self.assertIsNotNone(relionProt.outputVolume,
                                 "There was a problem with Relion initial model")
            self.assertAlmostEqual(outImgSet[1].getSamplingRate(),
                                   self.protImport.outputParticles[1].getSamplingRate(),
                                   msg="The sampling rate is wrong", delta=0.00001)

        relionProt = _runRelionIniModel(
            doGpu=USE_GPU, label="Relion initial model %s"
                                 % ('GPU' if USE_GPU else 'CPU'))
        _checkAsserts(relionProt)


class TestRelionLocalRes(TestRelionBase):
    """ This class helps to test local resolution protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        pathFns = 'import/refine3d/extra'
        cls.partFn = cls.ds.getFile(os.path.join(pathFns, 'relion_it025_data.star'))
        cls.protImport = cls.runImportParticlesStar(cls.partFn, 7.08)
        cls.volFn = cls.ds.getFile(os.path.join(pathFns, 'relion_class001.mrc'))
        cls.half1Fn = cls.ds.getFile(os.path.join(pathFns, 'relion_it025_half1_class001.mrc'))
        cls.half2Fn = cls.ds.getFile(os.path.join(pathFns, 'relion_it025_half2_class001.mrc'))
        cls.modelFn = cls.ds.getFile(os.path.join(pathFns, 'relion_model.star'))
        cls.protImportVol = cls.runImportVolumes(cls.volFn, 7.08)

    def _createRef3DProtBox(self, label):
        prot = self.newProtocol(ProtRelionRefine3D)
        self.saveProtocol(prot)

        prot.setObjLabel(label)
        pwutils.makePath(prot._getPath())
        pwutils.makePath(prot._getExtraPath())
        pwutils.makePath(prot._getTmpPath())

        prot.inputParticles.set(self.protImport.outputParticles)
        prot.referenceVolume.set(self.protImportVol.outputVolume)

        volume = Volume()
        volume.setFileName(prot._getExtraPath('test.mrc'))
        pxSize = prot.inputParticles.get().getSamplingRate()
        volume.setSamplingRate(pxSize)

        prot._defineOutputs(outputVolume=volume)
        prot.setStatus(STATUS_FINISHED)

        return prot

    def _validations(self, vol, dims, pxSize):
        self.assertIsNotNone(vol, "There was a problem with localres protocol ")
        xDim = vol.getXDim()
        sr = vol.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your volume is (%d)^3 "
                                     "and must be (%d)^3" % (xDim, dims))
        self.assertAlmostEqual(sr, pxSize, delta=0.0001,
                               msg="Pixel size of your volume is %0.5f and"
                                   " must be %0.5f" % (sr, pxSize))

    def test_runRelionLocalRes(self):
        protRef = self._createRef3DProtBox("auto-refine")

        protRef._createFilenameTemplates()
        volPath = protRef._getFileName('finalvolume', ref3d=1).split(':')[0]
        volHalf1 = protRef._getFileName('final_half1_volume', ref3d=1).split(':')[0]
        volHalf2 = protRef._getFileName('final_half2_volume', ref3d=1).split(':')[0]

        pwutils.copyFile(self.volFn, volPath)
        pwutils.copyFile(self.half1Fn, volHalf1)
        pwutils.copyFile(self.half2Fn, volHalf2)
        pwutils.copyFile(self.modelFn,
                         protRef._getExtraPath('relion_model.star'))

        protRef.outputVolume.setFileName(volPath)
        protRef.outputVolume.setHalfMaps([volHalf1, volHalf2])
        project = protRef.getProject()
        project._storeProtocol(protRef)

        restProt = self.newProtocol(ProtRelionLocalRes, protRefine=protRef)
        print(pwutils.magentaStr("\n==> Testing relion - local resolution:"))
        restProt.setObjLabel('Relion local resolution')

        self.launchProtocol(restProt)
        self._validations(restProt.outputVolume, 60, 7.08)

        # Add also a test case for the mask
        print(pwutils.magentaStr("\n==> Importing data - mask 3D:"))
        protMask = self.newProtocol(ProtRelionCreateMask3D,
                                    inputVolume=protRef.outputVolume,
                                    initialLowPassFilterA=30)
        self.launchProtocol(protMask)

        print(pwutils.magentaStr("\n==> Testing relion - local resolution (with mask):"))
        restProt = self.newProtocol(ProtRelionLocalRes,
                                    objLabel='relion localres (with mask)',
                                    protRefine=protRef,
                                    solventMask=protMask.outputMask)
        self.launchProtocol(restProt)


class TestRelionPostprocess(TestRelionBase):
    """ This class helps to test postprocess protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        pathFns = 'import/refine3d/extra'
        cls.volFn = cls.ds.getFile(os.path.join(pathFns, 'relion_class001.mrc'))
        cls.half1Fn = cls.ds.getFile(os.path.join(pathFns, 'relion_it025_half1_class001.mrc'))
        cls.half2Fn = cls.ds.getFile(os.path.join(pathFns, 'relion_it025_half2_class001.mrc'))
        cls.protImportVol = cls.runImportVolumes(cls.volFn, 3.0)

    def importPartsFromScipion(self):
        partFn = self.ds.getFile('import/particles.sqlite')
        protPart = self.newProtocol(ProtImportParticles,
                                    objLabel='Import Particles',
                                    importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                    sqliteFile=partFn,
                                    samplingRate=3,
                                    haveDataBeenPhaseFlipped=True
                                    )
        self.launchProtocol(protPart)
        return protPart

    def _createRef3DProtBox(self, label, protocol,
                            storeIter=False, iterN=0):
        from pyworkflow.protocol.constants import STATUS_FINISHED

        prot = self.newProtocol(protocol)
        self.saveProtocol(prot)

        prot.setObjLabel(label)
        pwutils.makePath(prot._getPath())
        pwutils.makePath(prot._getExtraPath())
        pwutils.makePath(prot._getTmpPath())

        prot.inputParticles.set(self.importPartsFromScipion().outputParticles)

        protClassName = prot.getClassName()
        outputVol = self.protImportVol.outputVolume
        if protClassName.startswith('ProtRelionRefine3D'):
            prot.referenceVolume.set(outputVol)
        elif protClassName.startswith('XmippProtProjMatch'):
            prot.input3DReferences.set(outputVol)
        elif protClassName.startswith('EmanProtRefine'):
            prot.input3DReference.set(outputVol)

        volume = Volume()
        volume.setFileName(prot._getExtraPath('test.mrc'))
        pxSize = prot.inputParticles.get().getSamplingRate()
        volume.setSamplingRate(pxSize)
        if storeIter:
            prot._setLastIter(iterN)
        prot._defineOutputs(outputVolume=volume)

        prot.setStatus(STATUS_FINISHED)

        # Create a mask protocol, because now it is not part of post-process
        protMask = self.newProtocol(ProtRelionCreateMask3D)
        protMask.inputVolume.set(outputVol)
        self.launchProtocol(protMask)

        return prot, protMask

    def _validations(self, vol, dims, pxSize, prot=""):
        self.assertIsNotNone(vol, "There was a problem with postprocess "
                                  "protocol, using %s protocol as input" % prot)
        xDim = vol.getXDim()
        sr = vol.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your volume is (%d)^3 "
                                     "and must be (%d)^3" % (xDim, dims))

        self.assertAlmostEqual(sr, pxSize, delta=0.0001,
                               msg="Pixel size of your volume is %0.5f and"
                                   " must be %0.5f" % (sr, pxSize))

    def test_postProcess_from_autorefine(self):
        print(pwutils.magentaStr("\n==> Testing relion - postprocess after refine 3d:"))
        protRef, protMask = self._createRef3DProtBox("auto-refine",
                                                     ProtRelionRefine3D)
        protRef._createFilenameTemplates()
        volPath = protRef._getFileName('finalvolume', ref3d=1).split(':')[0]
        volHalf1 = protRef._getFileName('final_half1_volume', ref3d=1).split(':')[0]
        volHalf2 = protRef._getFileName('final_half2_volume', ref3d=1).split(':')[0]

        pwutils.copyFile(self.volFn, volPath)
        pwutils.copyFile(self.half1Fn, volHalf1)
        pwutils.copyFile(self.half2Fn, volHalf2)

        protRef.outputVolume.setFileName(volPath)
        protRef.outputVolume.setHalfMaps([volHalf1, volHalf2])
        project = protRef.getProject()
        project._storeProtocol(protRef)

        postProt = self.newProtocol(ProtRelionPostprocess,
                                    protRefine=protRef,
                                    solventMask=protMask.outputMask)
        postProt.setObjLabel('post process Auto-refine')

        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 3, "Relion auto-refine")

    def test_postProcess_from_projMatch(self):
        print(pwutils.magentaStr("\n==> Testing relion - postprocess after xmipp proj. match.:"))
        XmippProtProjMatch = Domain.importFromPlugin('xmipp3.protocols',
                                                     'XmippProtProjMatch')

        if XmippProtProjMatch is None:
            print("WARNING: Can not load xmipp3.protocols.XmippProtProjMatch."
                  "Skipping the tests.")
            return

        protRef, protMask = self._createRef3DProtBox(
            "Proj Match", XmippProtProjMatch, storeIter=True, iterN=2)

        pwutils.makePath(os.path.join(protRef._getExtraPath(), 'iter_002'))
        protRef._initialize()
        volXmipp = protRef._getFileName('reconstructedFileNamesIters',
                                        iter=2, ref=1)
        half1Xmipp = protRef._getFileName('reconstructedFileNamesItersSplit1',
                                          iter=2, ref=1)
        half2Xmipp = protRef._getFileName('reconstructedFileNamesItersSplit2',
                                          iter=2, ref=1)

        ih = ImageHandler()
        ih.convert(ih.getVolFileName(self.volFn), volXmipp)
        ih.convert(ih.getVolFileName(self.half1Fn), half1Xmipp)
        ih.convert(ih.getVolFileName(self.half2Fn), half2Xmipp)

        protRef.outputVolume.setFileName(volXmipp)
        protRef.outputVolume.setHalfMaps([half1Xmipp, half2Xmipp])
        project = protRef.getProject()
        project._storeProtocol(protRef)

        postProt = self.newProtocol(ProtRelionPostprocess,
                                    protRefine=protRef,
                                    solventMask=protMask.outputMask)
        postProt.setObjLabel('post process Xmipp Projection Matching')
        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 3, "Projection Matching")

    def test_postProcess_from_eman_refineEasy(self):
        print(pwutils.magentaStr("\n==> Testing relion - postprocess after eman2 refine easy:"))
        EmanProtRefine = Domain.importFromPlugin('eman2.protocols', 'EmanProtRefine')
        convertImage = Domain.importFromPlugin('eman2.convert', 'convertImage')

        protRef, protMask = self._createRef3DProtBox(
            "Eman refine Easy", EmanProtRefine)

        pwutils.makePath(os.path.join(protRef._getExtraPath(), 'refine_01'))
        protRef._createFilenameTemplates()
        volEman = protRef._getFileName("mapFull", run=1, iter=2)
        half1Eman = protRef._getFileName("mapEvenUnmasked", run=1)
        half2Eman = protRef._getFileName("mapOddUnmasked", run=1)

        convertImage(self.volFn, volEman)
        convertImage(self.half1Fn, half1Eman)
        convertImage(self.half2Fn, half2Eman)

        protRef.outputVolume.setFileName(volEman)
        protRef.outputVolume.setHalfMaps([half1Eman, half2Eman])
        project = protRef.getProject()
        project._storeProtocol(protRef)

        postProt = self.newProtocol(ProtRelionPostprocess,
                                    protRefine=protRef,
                                    solventMask=protMask.outputMask)
        postProt.setObjLabel('post process Eman2 refine-easy')
        self.launchProtocol(postProt)
        self._validations(postProt.outputVolume, 60, 3, "Eman refine easy")


class TestRelionRefine(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)

    def testProtRelionRefine(self):
        def _runRelionRefine(doGpu=False, label=''):
            relionRefine = self.newProtocol(ProtRelionRefine3D,
                                            doCTF=False, runMode=1,
                                            maskDiameterA=340,
                                            symmetryGroup="d6",
                                            numberOfMpi=3, numberOfThreads=2)
            relionRefine.setObjLabel(label)
            relionRefine.inputParticles.set(self.protNormalize.outputParticles)
            relionRefine.referenceVolume.set(self.protImportVol.outputVolume)
            relionRefine.doGpu.set(doGpu)
            self.launchProtocol(relionRefine)
            return relionRefine

        def _checkAsserts(relionRefine):
            relionRefine._initialize()  # Load filename templates
            dataSqlite = relionRefine._getIterData(3)
            outImgSet = SetOfParticles(filename=dataSqlite)

            self.assertIsNotNone(relionRefine.outputVolume,
                                 "There was a problem with Relion autorefine")
            self.assertAlmostEqual(outImgSet[1].getSamplingRate(),
                                   self.protNormalize.outputParticles[1].getSamplingRate(),
                                   msg="The sampling rate is wrong", delta=0.00001)

            self.assertAlmostEqual(outImgSet[1].getFileName(),
                                   self.protNormalize.outputParticles[1].getFileName(),
                                   msg="The particles filenames are wrong")

        if RUN_CPU:
            print(pwutils.magentaStr("\n==> Testing relion - refine 3d on CPU:"))
            relionProt = _runRelionRefine(doGpu=False,
                                          label="Run Relion auto-refine CPU")
            _checkAsserts(relionProt)

        if USE_GPU:
            print(pwutils.magentaStr("\n==> Testing relion - refine 3d on GPU:"))
            relionGpu = _runRelionRefine(doGpu=True,
                                         label="Run Relion auto-refine GPU")
            _checkAsserts(relionGpu)


class TestRelionResizeVolume(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.vol = cls.ds.getFile('import/refine3d/extra/relion_class001.mrc')
        cls.protImport = cls.runImportVolumes(cls.vol, 3.0)

    def _validations(self, vol, dims, pxSize):
        self.assertIsNotNone(vol, "There was a problem with crop/resize "
                                  "volumes protocol")
        xDim = vol.getXDim()
        sr = vol.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your volume is (%d)^3 "
                                     "and must be (%d)^3" % (xDim, dims))

        self.assertAlmostEqual(sr, pxSize, delta=0.0001,
                               msg="Pixel size of your volume is %0.5f and"
                                   " must be %0.5f" % (sr, pxSize))

    def test_resizeVol(self):
        resizeProt = self.newProtocol(ProtRelionResizeVolume,
                                      doRescale=True, rescaleSamplingRate=1.5,
                                      doResize=True, resizeSize=128)
        vol = self.protImport.outputVolume
        resizeProt.inputVolumes.set(vol)
        print(pwutils.magentaStr("\n==> Testing relion - crop/resize volumes:"))
        self.launchProtocol(resizeProt)

        self._validations(resizeProt.outputVol, 128, 1.5)


class TestRelionSubtract(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)

    def test_subtract(self):
        print(pwutils.magentaStr("\n==> Running relion - refine 3d:"))
        relionRefine = self.newProtocol(ProtRelionRefine3D,
                                        doCTF=False, runMode=1,
                                        maskDiameterA=340,
                                        symmetryGroup="d6",
                                        numberOfMpi=3, numberOfThreads=2)
        relionRefine.inputParticles.set(self.protImport.outputParticles)
        relionRefine.referenceVolume.set(self.protImportVol.outputVolume)
        relionRefine.doGpu.set(False)
        self.launchProtocol(relionRefine)

        print(pwutils.magentaStr("\n==> Running relion - create mask 3d:"))
        protMask = self.newProtocol(ProtRelionCreateMask3D, threshold=0.045)
        protMask.inputVolume.set(relionRefine.outputVolume)
        self.launchProtocol(protMask)

        print(pwutils.magentaStr("\n==> Testing relion - subtract projection:"))
        protSubtract = self.newProtocol(ProtRelionSubtract,
                                        refMask=protMask.outputMask,
                                        numberOfMpi=2)
        protSubtract.inputProtocol.set(relionRefine)
        self.launchProtocol(protSubtract)
        self.assertIsNotNone(protSubtract.outputParticles,
                             "There was a problem with subtract projection")


class TestRelionSymmetrizeVolume(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('resmap')
        cls.protImportVol = cls.runImportVolumes(cls.ds.getFile('betaGal.mrc'), 3.54)

    def test_symmetrizeVolume(self):
        print(pwutils.magentaStr("\n==> Testing relion - symmetrize volume:"))
        symmVol = self.newProtocol(
            ProtRelionSymmetrizeVolume,
            objLabel='symmetryze d2',
            inputVolume=self.protImportVol.outputVolume,
            symmetryGroup='d2'
        )
        self.launchProtocol(symmVol)

        self.assertIsNotNone(getattr(symmVol, 'outputVolumeAligned'))
        self.assertIsNotNone(getattr(symmVol, 'outputVolumeSymmetrized'))
