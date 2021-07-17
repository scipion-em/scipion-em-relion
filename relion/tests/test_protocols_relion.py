# **************************************************************************
# *
# * Authors:    Laura del Cano (ldelcano@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.plugin import Domain
from pyworkflow.protocol.constants import STATUS_FINISHED
from pyworkflow.utils import magentaStr
from pwem import Config
from pwem.objects import SetOfParticles
from pwem.protocols import *

from ..protocols import *
from ..convert import *
from ..constants import *


def useGpu():
    """ Helper function to determine if GPU can be used.
    Return a boolean and a label to be used in protocol's label. """
    cudaPath = Plugin.getVar('RELION_CUDA_LIB', Config.CUDA_LIB)

    if cudaPath and os.path.exists(cudaPath):
        return True, 'GPU'
    else:
        return False, 'CPU'


USE_GPU = useGpu()[0]
ONLY_GPU = int(os.environ.get('SCIPION_TEST_RELION_ONLY_GPU', 0))
RUN_CPU = not USE_GPU or ONLY_GPU


class TestRelionBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='mda'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.particlesFn = cls.dataset.getFile('particles')
        cls.vol = cls.dataset.getFile('volumes')

    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an output was generated and
        the condition is valid.
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)
    
    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False):
        """ Run an Import particles protocol. """
        print(magentaStr("\n==> Importing data - particles from files:"))
        protImport = cls.newProtocol(ProtImportParticles, 
                                     filesPath=pattern,
                                     samplingRate=samplingRate,
                                     checkStack=checkStack)
        cls.launchProtocol(protImport)
        cls.assertIsNotNone(protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return protImport

    @classmethod
    def runImportParticlesStar(cls, partStar, mag, samplingRate):
        """ Import particles from Relion star file. """
        print(magentaStr("\n==> Importing data - particles from star:"))
        protImport = cls.newProtocol(ProtImportParticles,
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=partStar,
                                     magnification=mag,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=True
                                     )
        cls.launchProtocol(protImport)
        cls.assertIsNotNone(protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return protImport

    @classmethod
    def runNormalizeParticles(cls, particles):
        """ Run normalize particles protocol """
        print(magentaStr("\n==> Running relion - preprocess particles:"))
        protPreproc = cls.newProtocol(ProtRelionPreprocessParticles,
                                      doNormalize=True)
        protPreproc.inputParticles.set(particles)
        cls.launchProtocol(protPreproc)
        return protPreproc
    
    @classmethod
    def runImportVolumes(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        print(magentaStr("\n==> Importing data - volumes:"))
        protImport = cls.newProtocol(ProtImportVolumes, 
                                     filesPath=pattern,
                                     samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        return protImport

    @classmethod
    def runImportMovies(cls, pattern, mag, samplingRate, dose):
        """ Run an Import movies protocol. """
        print(magentaStr("\n==> Importing data - movies:"))
        protImport = cls.newProtocol(ProtImportMovies,
                                     filesPath=pattern,
                                     magnification=mag,
                                     samplingRate=samplingRate,
                                     dosePerFrame=dose)
        cls.launchProtocol(protImport)
        return protImport


class TestRelionPicking(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.ds.getFile('import/classify2d/extra/relion_it015_data.star')

        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImportMics = cls.newProtocol(
            ProtImportMicrographs,
            samplingRateMode=0,
            filesPath='%s/*.mrc' % cls.ds.getFile('micrographs'),
            samplingRate=7.08,
            magnification=50000,
            voltage=300,
            sphericalAberration=0.1)
        cls.launchProtocol(cls.protImportMics)

    def _checkOutput(self, pickProt, minCoords, maxCoords):
        """ Check that the outputCoordinates is not None
         and that it should be between 300 and 400 coordinates
         per micrograph.
        """
        coordSet = getattr(pickProt, 'outputCoordinates', None)
        self.assertIsNotNone(coordSet)
        for micAgg in coordSet.aggregate(["count"], "_micId", ["_micId"]):
            self.assertGreaterEqual(micAgg['count'], minCoords)
            self.assertLessEqual(micAgg['count'], maxCoords)

    def testPickingLog(self):
        print(magentaStr("\n==> Testing relion - autopick LoG:"))
        protPickLog = self.newProtocol(
            ProtRelionAutopickLoG,
            objLabel='autopick LoG',
            inputMicrographs=self.protImportMics.outputMicrographs,
            boxSize=64,
            minDiameter=260,
            maxDiameter=360,
            streamingBatchSize=5,
        )
        self.launchProtocol(protPickLog)
        self._checkOutput(protPickLog, 300, 400)

    def testPickingRef(self):
        # Create a subset with a few good averages
        ih = ImageHandler()
        avgsFn = self.ds.getFile('import/classify2d/extra/'
                                 'relion_it015_classes.mrcs')
        outAvgsFn = os.path.abspath(self.proj.getTmpPath('averages.mrcs'))

        for i, index in enumerate([5, 16, 17, 18, 31]):
            ih.convert((index, avgsFn), (i + 1, outAvgsFn))

        print(magentaStr("\n==> Importing data - averages:"))
        protAvg = self.newProtocol(ProtImportAverages,
                                   importFrom=ProtImportParticles.IMPORT_FROM_FILES,
                                   filesPath=outAvgsFn,
                                   samplingRate=7.08
                                   )
        self.launchProtocol(protAvg)

        # We need CTF estimation for picking ref with Relion
        # Now estimate CTF on the micrographs with ctffind
        ProtCTFFind = Domain.importFromPlugin(
            'cistem.protocols', 'CistemProtCTFFind', doRaise=True)

        print(magentaStr("\n==> Running cistem - ctffind:"))
        protCtf = self.newProtocol(
            ProtCTFFind,
            inputMicrographs=self.protImportMics.outputMicrographs,
            minDefocus=12000, maxDefocus=30000,
            slowSearch=False,
            resamplePix=False,
        )
        self.launchProtocol(protCtf)

        print(magentaStr("\n==> Testing relion - autopick ref:"))
        protPickRef = self.newProtocol(
            ProtRelion2Autopick,
            inputMicrographs=self.protImportMics.outputMicrographs,
            ctfRelations=protCtf.outputCTF,
            inputReferences=protAvg.outputAverages,
            streamingBatchSize=5,
        )
        self.launchProtocol(protPickRef)
        self._checkOutput(protPickRef, 240, 330)


class TestRelionClassify2D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protNormalize = cls.runNormalizeParticles(cls.protImport.outputParticles)
    
    def testRelion2D(self):
        def _runRelionClassify2D(doGpu=False, label=''):
            print(magentaStr("\n==> Testing relion - classify 2d:"))
            prot2D = self.newProtocol(ProtRelionClassify2D,
                                      doCTF=False, maskDiameterA=340,
                                      numberOfMpi=4, numberOfThreads=1)
            prot2D.numberOfClasses.set(4)
            prot2D.numberOfIterations.set(3)
            prot2D.inputParticles.set(self.protNormalize.outputParticles)
            prot2D.setObjLabel(label)
            prot2D.doGpu.set(doGpu)
            self.launchProtocol(prot2D)
            return prot2D

        def _checkAsserts(relionProt):

            self.assertIsNotNone(relionProt.outputClasses,
                                 "There was a problem with Relion 2D classify")
        
            partsPixSize = self.protNormalize.outputParticles.getSamplingRate()
            classsesPixSize = relionProt.outputClasses.getImages().getSamplingRate()
            self.assertAlmostEquals(partsPixSize, classsesPixSize,
                                    msg="There was a problem with the sampling rate "
                                    "of the particles", delta=0.001)

            for class2D in relionProt.outputClasses:
                self.assertTrue(class2D.hasAlignment2D())

        if RUN_CPU:
            relionProt = _runRelionClassify2D(doGpu=False,
                                              label="Run Relion classify2D CPU")
            _checkAsserts(relionProt)

        if USE_GPU:
            relionGpu = _runRelionClassify2D(doGpu=True,
                                             label="Relion classify2D GPU")
            _checkAsserts(relionGpu)


class TestRelionClassify3D(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
    
    def testProtRelionClassify3D(self):
        print(magentaStr("\n==> Running relion - preprocess particles:"))
        relionNormalize = self.newProtocol(ProtRelionPreprocessParticles)
        relionNormalize.inputParticles.set(self.protImport.outputParticles)
        relionNormalize.doNormalize.set(True)
        self.launchProtocol(relionNormalize)

        def _runRelionClassify3D(doGpu=False, label=''):
            relion3DClass = self.newProtocol(ProtRelionClassify3D,
                                             numberOfClasses=3,
                                             numberOfIterations=4,
                                             doCTF=False, runMode=1,
                                             maskDiameterA=320,
                                             numberOfMpi=2, numberOfThreads=2)

            relion3DClass.setObjLabel(label)
            relion3DClass.inputParticles.set(relionNormalize.outputParticles)
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
            print(magentaStr("\n==> Testing relion - classify 3d on CPU:"))
            relionProt = _runRelionClassify3D(doGpu=False,
                                              label="Run Relion classify3D CPU")
            _checkAsserts(relionProt)

        if USE_GPU:
            print(magentaStr("\n==> Testing relion - classify 3d on GPU:"))
            relionGpu = _runRelionClassify3D(doGpu=True,
                                             label="Relion classify3D GPU")
            _checkAsserts(relionGpu)


class TestRelionRefine(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)
    
    def testProtRelionRefine(self):
        print(magentaStr("\n==> Running relion - preprocess particles:"))
        relNorm = self.newProtocol(ProtRelionPreprocessParticles)
        relNorm.inputParticles.set(self.protImport.outputParticles)
        relNorm.doNormalize.set(True)
        self.launchProtocol(relNorm)
        
        def _runRelionRefine(doGpu=False, label=''):
            relionRefine = self.newProtocol(ProtRelionRefine3D,
                                            doCTF=False, runMode=1,
                                            maskDiameterA=340,
                                            symmetryGroup="d6",
                                            numberOfMpi=3, numberOfThreads=2)
            relionRefine.setObjLabel(label)
            relionRefine.inputParticles.set(relNorm.outputParticles)
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
                                   relNorm.outputParticles[1].getSamplingRate(),
                                   msg="The sampling rate is wrong", delta=0.00001)
            
            self.assertAlmostEqual(outImgSet[1].getFileName(),
                                   relNorm.outputParticles[1].getFileName(),
                                   msg="The particles filenames are wrong")
        
        if RUN_CPU:
            print(magentaStr("\n==> Testing relion - refine 3d on CPU:"))
            relionProt = _runRelionRefine(doGpu=False,
                                          label="Run Relion auto-refine CPU")
            _checkAsserts(relionProt)

        if USE_GPU:
            print(magentaStr("\n==> Testing relion - refine 3d on GPU:"))
            relionGpu = _runRelionRefine(doGpu=True,
                                         label="Run Relion auto-refine GPU")
            _checkAsserts(relionGpu)


class TestRelionInitialModel(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.dataset.getFile('import/classify2d/extra/relion_it015_data.star')
        cls.protImport = cls.runImportParticlesStar(cls.partFn, 50000, 7.08)

    def testProtRelionIniModel(self):
        def _runRelionIniModel(doGpu=True, label=''):
            kwargs = {
                'doCTF': False,
                'doGpu': doGpu,
                'maskDiameterA': 340,
                'symmetryGroup': 'C1',
                'allParticlesRam': True,
                'numberOfMpi': 3,
                'numberOfThreads': 2
            }
            kwargs.update({'numberOfIterInitial': 10,
                           'numberOfIterInBetween': 30,
                           'numberOfIterFinal': 10})

            print(magentaStr("\n==> Testing relion - initial model:"))
            relionIniModel = self.newProtocol(ProtRelionInitialModel, **kwargs)
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

        
class TestRelionPreprocess(TestRelionBase):
    """ This class helps to test all different preprocessing particles options
    on Relion. """
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)

    def _validations(self, imgSet, dims, pxSize):
        self.assertIsNotNone(imgSet, "There was a problem with preprocess "
                                     "particles")
        xDim = imgSet.getXDim()
        sr = imgSet.getSamplingRate()
        self.assertEqual(xDim, dims, "The dimension of your particles are %d x "
                                     "%d and must be  %d x %d" % (xDim, xDim,
                                                                  dims, dims))
        self.assertAlmostEqual(sr, pxSize, delta=0.0001,
                               msg="Pixel size of your particles are  %0.5f and"
                               " must be %0.5f" % (sr, pxSize))

    def test_NormalizeAndDust(self):
        print(magentaStr("\n==> Testing relion - preprocess particles (norm, remove dust):"))
        protocol = self.newProtocol(ProtRelionPreprocessParticles,
                                    doNormalize=True, backRadius=40,
                                    doRemoveDust=True, whiteDust=4, blackDust=8)
        protocol.setObjLabel('relion: norm-dust')
        protocol.inputParticles.set(self.protImport.outputParticles)
        self.launchProtocol(protocol)
        
        self._validations(protocol.outputParticles, 100, 3.5)
        
    def test_ScaleAndInvert(self):
        print(magentaStr("\n==> Testing relion - preprocess particles (scale, invert):"))
        protocol = self.newProtocol(ProtRelionPreprocessParticles,
                                    doNormalize=False,
                                    doScale=True, scaleSize=50,
                                    doInvert=True)
        protocol.setObjLabel('relion: scale-invert')
        protocol.inputParticles.set(self.protImport.outputParticles)
        
        self.launchProtocol(protocol)
        self._validations(protocol.outputParticles, 50, 7.0)


class TestRelionSubtract(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestRelionBase.setData('mda')
        cls.protImport = cls.runImportParticles(cls.particlesFn, 3.5)
        cls.protImportVol = cls.runImportVolumes(cls.vol, 3.5)

    def test_subtract(self):
        print(magentaStr("\n==> Running relion - refine 3d:"))
        relionRefine = self.newProtocol(ProtRelionRefine3D,
                                        doCTF=False, runMode=1,
                                        maskDiameterA=340,
                                        symmetryGroup="d6",
                                        numberOfMpi=3, numberOfThreads=2)
        relionRefine.inputParticles.set(self.protImport.outputParticles)
        relionRefine.referenceVolume.set(self.protImportVol.outputVolume)
        relionRefine.doGpu.set(False)
        self.launchProtocol(relionRefine)

        print(magentaStr("\n==> Running relion - create mask 3d:"))
        protMask = self.newProtocol(ProtRelionCreateMask3D, threshold=0.045)
        protMask.inputVolume.set(relionRefine.outputVolume)
        self.launchProtocol(protMask)

        print(magentaStr("\n==> Testing relion - subtract projection:"))
        protSubtract = self.newProtocol(ProtRelionSubtract,
                                        refMask=protMask.outputMask,
                                        numberOfMpi=2)
        protSubtract.inputProtocol.set(relionRefine)
        self.launchProtocol(protSubtract)
        self.assertIsNotNone(protSubtract.outputParticles,
                             "There was a problem with subtract projection")


class TestRelionPostprocess(TestRelionBase):
    """ This class helps to test postprocess protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        pathFns = 'import/refine3d/extra'
        cls.volFn = cls.ds.getFile(join(pathFns, 'relion_class001.mrc'))
        cls.half1Fn = cls.ds.getFile(join(pathFns, 'relion_it025_half1_class001.mrc'))
        cls.half2Fn = cls.ds.getFile(join(pathFns, 'relion_it025_half2_class001.mrc'))

    def importVolume(self):

        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=self.volFn,
                                   samplingRate=3)
        self.launchProtocol(protVol)
        return protVol

    def importPartsFromScipion(self):
        partFn = self.ds.getFile('import/particles.sqlite')
        protPart = self.newProtocol(ProtImportParticles,
                                    objLabel='Import Particles',
                                    importFrom=ProtImportParticles.IMPORT_FROM_SCIPION,
                                    sqliteFile=partFn,
                                    magnification=10000,
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
        outputVol = self.importVolume().outputVolume
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
        print(magentaStr("\n==> Testing relion - postprocess after refine 3d:"))
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
        print(magentaStr("\n==> Testing relion - postprocess after xmipp proj. match.:"))
        XmippProtProjMatch = Domain.importFromPlugin('xmipp3.protocols',
                                                     'XmippProtProjMatch')

        if XmippProtProjMatch is None:
            print("WARNING: Can not load xmipp3.protocols.XmippProtProjMatch."
                  "Skipping the tests.")
            return

        protRef, protMask = self._createRef3DProtBox(
            "Proj Match", XmippProtProjMatch, storeIter=True, iterN=2)

        pwutils.makePath(join(protRef._getExtraPath(), 'iter_002'))
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
        print(magentaStr("\n==> Testing relion - postprocess after eman2 refine easy:"))
        EmanProtRefine = Domain.importFromPlugin('eman2.protocols', 'EmanProtRefine')
        convertImage = Domain.importFromPlugin('eman2.convert', 'convertImage')

        protRef, protMask = self._createRef3DProtBox(
            "Eman refine Easy", EmanProtRefine)

        pwutils.makePath(join(protRef._getExtraPath(), 'refine_01'))
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


class TestRelionLocalRes(TestRelionBase):
    """ This class helps to test local resolution protocol from Relion. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        pathFns = 'import/refine3d/extra'
        cls.partFn = cls.ds.getFile(join(pathFns, 'relion_it025_data.star'))
        cls.protImport = cls.runImportParticlesStar(cls.partFn, 50000, 7.08)
        cls.volFn = cls.ds.getFile(join(pathFns, 'relion_class001.mrc'))
        cls.half1Fn = cls.ds.getFile(join(pathFns, 'relion_it025_half1_class001.mrc'))
        cls.half2Fn = cls.ds.getFile(join(pathFns, 'relion_it025_half2_class001.mrc'))
        cls.modelFn = cls.ds.getFile(join(pathFns, 'relion_model.star'))

    def importVolume(self):
        print(magentaStr("\n==> Importing data - volume:"))
        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=self.volFn,
                                   samplingRate=7.08)
        self.launchProtocol(protVol)
        return protVol

    def _createRef3DProtBox(self, label):
        prot = self.newProtocol(ProtRelionRefine3D)
        self.saveProtocol(prot)

        prot.setObjLabel(label)
        pwutils.makePath(prot._getPath())
        pwutils.makePath(prot._getExtraPath())
        pwutils.makePath(prot._getTmpPath())

        prot.inputParticles.set(self.protImport.outputParticles)
        prot.referenceVolume.set(self.importVolume().outputVolume)

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
        print(magentaStr("\n==> Testing relion - local resolution:"))
        restProt.setObjLabel('Relion local resolution')

        self.launchProtocol(restProt)
        self._validations(restProt.outputVolume, 60, 7.08)

        # Add also a test case for the mask
        print(magentaStr("\n==> Importing data - mask 3D:"))
        protMask = self.newProtocol(ProtRelionCreateMask3D,
                                    inputVolume=protRef.outputVolume,
                                    initialLowPassFilterA=30)
        self.launchProtocol(protMask)

        print(magentaStr("\n==> Testing relion - local resolution (with mask):"))
        restProt = self.newProtocol(ProtRelionLocalRes,
                                    objLabel='relion localres (with mask)',
                                    protRefine=protRef,
                                    solventMask=protMask.outputMask)
        self.launchProtocol(restProt)


class TestRelionExpandSymmetry(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partRef3dFn = cls.dataset.getFile('import/refine3d/extra/relion_data.star')

    def importParticles(self, partStar):
        """ Import particles from Relion star file. """
        print(magentaStr("\n==> Importing data - particles:"))
        protPart = self.newProtocol(ProtImportParticles,
                                    importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                    starFile=partStar,
                                    magnification=10000,
                                    samplingRate=7.08,
                                    haveDataBeenPhaseFlipped=False
                                    )
        self.launchProtocol(protPart)
        return protPart

    def test_ExpandSymmetry(self):
        print(magentaStr("\n==> Testing relion - expand symmetry:"))
        prot = self.newProtocol(ProtRelionExpandSymmetry)
        importRun = self.importParticles(self.partRef3dFn)
        prot.inputParticles.set(importRun.outputParticles)
        prot.symmetryGroup.set("D2")
        self.launchProtocol(prot)

        self.assertIsNotNone(prot.outputParticles,
                             "There was a problem with expand symmetry protocol")
        sizeIn = importRun.outputParticles.getSize()
        sizeOut = prot.outputParticles.getSize()
        self.assertAlmostEqual(sizeIn * 4, sizeOut, delta=0.0001,
                               msg="Number of output particles is %d and"
                               " must be %d" % (sizeOut, sizeIn * 4))


class TestRelionCreate3dMask(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def importVolume(self):
        print(magentaStr("\n==> Importing data - volume:"))
        volFn = self.ds.getFile('import/refine3d/extra/relion_class001.mrc')
        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=volFn,
                                   samplingRate=3)
        self.launchProtocol(protVol)
        return protVol

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
        importProt = self.importVolume()

        maskProt = self.newProtocol(ProtRelionCreateMask3D,
                                    initialLowPassFilterA=10)  # filter at 10 A
        vol = importProt.outputVolume
        maskProt.inputVolume.set(vol)
        print(magentaStr("\n==> Testing relion - create mask 3d:"))
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


class TestRelionSymmetrizeVolume(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('resmap')

    def test_symmetrizeVolume(self):
        print(magentaStr("\n==> Importing data - volume:"))
        importVol = self.newProtocol(
            ProtImportVolumes,
            objLabel='import volume',
            filesPath=self.ds.getFile('betaGal.mrc'),
            samplingRate=3.54)

        self.launchProtocol(importVol)

        print(magentaStr("\n==> Testing relion - symmetrize volume:"))
        symmVol = self.newProtocol(
            ProtRelionSymmetrizeVolume,
            objLabel='symmetryze d2',
            inputVolume=importVol.outputVolume,
            symmetryGroup='d2'
        )
        self.launchProtocol(symmVol)

        self.assertIsNotNone(getattr(symmVol, 'outputVolumeAligned'))
        self.assertIsNotNone(getattr(symmVol, 'outputVolumeSymmetrized'))


class TestRelionExtractParticles(TestRelionBase):
    """This class check if the protocol to extract particles
    in Relion works properly.
    """

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """

        # We have two options:
        # 1) pass the SamplingRate or
        # 2) the ScannedPixelSize + microscope magnification
        if samplingRate is not None:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=0,
                                             filesPath=pattern,
                                             samplingRate=samplingRate,
                                             magnification=magnification,
                                             voltage=voltage,
                                             sphericalAberration=sphericalAberration)
        else:
            cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                             samplingRateMode=1,
                                             filesPath=pattern,
                                             scannedPixelSize=scannedPixelSize,
                                             voltage=voltage,
                                             magnification=magnification,
                                             sphericalAberration=sphericalAberration)

        cls.protImport.setObjLabel('import mics')
        cls.launchProtocol(cls.protImport)
        cls.assertIsNotNone(cls.protImport.outputMicrographs,
                            "SetOfMicrographs has not been produced.")

        return cls.protImport

    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        print(magentaStr("\n==> Importing data - micrographs:"))
        return cls.runImportMicrograph(pattern, samplingRate=1.237,
                                       voltage=300, sphericalAberration=2,
                                       scannedPixelSize=None, magnification=56000)

    @classmethod
    def runDownsamplingMicrographs(cls, mics, downFactorValue, threads=1):
        # test downsampling a set of micrographs
        XmippProtPreprocessMicrographs = Domain.importFromPlugin(
            'xmipp3.protocols', 'XmippProtPreprocessMicrographs')

        if XmippProtPreprocessMicrographs is None:
            print("WARNING: Can not load xmipp3.protocols.XmippProtPreprocessMicrographs."
                  "Skipping the tests.")
            return

        print(magentaStr("\n==> Running xmipp - preprocess micrographs:"))
        cls.protDown = XmippProtPreprocessMicrographs(doDownsample=True,
                                                      downFactor=downFactorValue,
                                                      numberOfThreads=threads)
        cls.protDown.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protDown, wait=True)
        return cls.protDown

    @classmethod
    def runFakedPicking(cls, mics, pattern):
        """ Run a faked particle picking. Coordinates already existing. """

        # TODO This fake picking depends on Xmipp particle picking
        # TODO: Can we change this to an import coordinates?

        XmippProtParticlePicking = Domain.importFromPlugin(
            'xmipp3.protocols', 'XmippProtParticlePicking')

        if XmippProtParticlePicking is None:
            print("WARNING: Can not load xmipp3.protocols.XmippProtParticlePicking."
                  "Skipping the tests.")
            return None

        print(magentaStr("\n==> Running xmipp - fake particle picking:"))
        cls.protPP = XmippProtParticlePicking(importFolder=pattern, runMode=1)
        cls.protPP.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protPP, wait=True)
        # check that faked picking has run ok
        cls.assertIsNotNone(cls.protPP.outputCoordinates,
                            "SetOfCoordinates has not been produced.")

        return cls.protPP

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('xmipp_tutorial')
        cls.micFn = cls.dataset.getFile('mic1')
        cls.micsFn = cls.dataset.getFile('allMics')
        cls.coordsDir = cls.dataset.getFile('posSupervisedDir')
        cls.allCrdsDir = cls.dataset.getFile('posAllDir')

        cls.DOWNSAMPLING = 5.0
        cls.protImport = cls.runImportMicrographBPV(cls.micsFn)
        cls.protDown = cls.runDownsamplingMicrographs(cls.protImport.outputMicrographs,
                                                      cls.DOWNSAMPLING)

        print(magentaStr("\n==> Importing data - ctfs:"))
        cls.protCTF = cls.newProtocol(ProtImportCTF,
                                      importFrom=ProtImportCTF.IMPORT_FROM_XMIPP3,
                                      filesPath=cls.dataset.getFile('ctfsDir'),
                                      filesPattern='*.ctfparam')
        cls.protCTF.inputMicrographs.set(cls.protImport.outputMicrographs)
        cls.proj.launchProtocol(cls.protCTF, wait=True)

        cls.protPP = cls.runFakedPicking(cls.protDown.outputMicrographs, cls.allCrdsDir)

    def _checkSamplingConsistency(self, outputSet):
        """ Check that the set sampling is the same as item sampling. """
        first = outputSet.getFirstItem()

        self.assertAlmostEqual(outputSet.getSamplingRate(),
                               first.getSamplingRate())

    def testExtractSameAsPicking(self):
        print(magentaStr("\n==> Testing relion - extract particles (no ctf):"))
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=110,
                                       doInvert=False)
        protExtract.setObjLabel("extract-noctf")
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        micSampling = protExtract.inputCoordinates.get().getMicrographs().getSamplingRate()

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertAlmostEqual(outputParts.getSamplingRate() / micSampling,
                               1, 1,
                               "There was a problem generating the output.")
        self._checkSamplingConsistency(outputParts)

        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(83)
        compare(228)

    def testExtractOriginal(self):
        print(magentaStr("\n==> Testing relion - extract particles (other mics):"))
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=550,
                                       downsampleType=OTHER,
                                       doInvert=False)
        protExtract.setObjLabel("extract-other")
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate()
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords

        def compare(objId, delta=1.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx / factor, px, delta=delta)
            self.assertAlmostEquals(cy / factor, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(111)
        compare(7)

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertEqual(outputParts.getSamplingRate(), samplingMics,
                         "Output sampling rate should be equal to input "
                         "sampling rate.")
        self._checkSamplingConsistency(outputParts)

    def testExtractOther(self):
        print(magentaStr("\n==> Testing relion - extract particles (other mics with downsampling):"))
        downFactor = 2.989
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=550, downsampleType=OTHER,
                                       doRescale=True,
                                       rescaledSize=184,
                                       doInvert=False,
                                       doFlip=False)
        # Get all the micrographs ids to validate that all particles
        # has the micId properly set
        micsId = [mic.getObjId() for mic in
                  self.protPP.outputCoordinates.getMicrographs()]

        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.inputMicrographs.set(self.protImport.outputMicrographs)
        protExtract.setObjLabel("extract-other+downsample")
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles
        samplingCoords = self.protPP.outputCoordinates.getMicrographs().getSamplingRate()
        samplingFinal = self.protImport.outputMicrographs.getSamplingRate() * downFactor
        samplingMics = protExtract.inputMicrographs.get().getSamplingRate()
        factor = samplingFinal / samplingCoords
        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")

        def compare(objId, delta=2.0):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx / factor, px, delta=delta)
            self.assertAlmostEquals(cy / factor, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(45)
        compare(229)

        outputSampling = outputParts.getSamplingRate()
        self.assertAlmostEqual(outputSampling / samplingMics,
                               downFactor, 1,
                               "There was a problem generating the output.")
        for particle in outputParts:
            self.assertTrue(particle.getCoordinate().getMicId() in micsId)
            self.assertAlmostEqual(outputSampling, particle.getSamplingRate())

    def testExtractCTF(self):
        print(magentaStr("\n==> Testing relion - extract particles (with CTF):"))
        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=110,
                                       downsampleType=SAME_AS_PICKING,
                                       doInvert=False,
                                       doFlip=True)
        protExtract.inputCoordinates.set(self.protPP.outputCoordinates)
        protExtract.ctfRelations.set(self.protCTF.outputCTF)
        protExtract.setObjLabel("extract-ctf")
        self.launchProtocol(protExtract)

        inputCoords = protExtract.inputCoordinates.get()
        outputParts = protExtract.outputParticles

        def compare(objId, delta=0.001):
            cx, cy = inputCoords[objId].getPosition()
            px, py = outputParts[objId].getCoordinate().getPosition()
            micNameCoord = inputCoords[objId].getMicName()
            micNamePart = outputParts[objId].getCoordinate().getMicName()
            self.assertAlmostEquals(cx, px, delta=delta)
            self.assertAlmostEquals(cy, py, delta=delta)
            self.assertEqual(micNameCoord, micNamePart,
                             "The micName should be %s and its %s"
                             % (micNameCoord, micNamePart))

        compare(228)
        compare(83)

        def compareCTF(partId, ctfId):
            partDefU = outputParts[partId].getCTF().getDefocusU()
            defU = protExtract.ctfRelations.get()[ctfId].getDefocusU()
            self.assertAlmostEquals(partDefU, defU, delta=1)

        compareCTF(1, 1)
        compareCTF(150, 2)
        compareCTF(300, 3)

        self.assertIsNotNone(outputParts,
                             "There was a problem generating the output.")
        self.assertTrue(outputParts.hasCTF(), "Output does not have CTF.")
        self._checkSamplingConsistency(outputParts)


class TestRelionCenterAverages(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('mda')

    def test_basic(self):
        """ Run an Import averages protocol. """
        print(magentaStr("\n==> Importing data - averages:"))
        protImport = self.newProtocol(ProtImportAverages,
                                      filesPath=self.ds.getFile('averages/averages.stk'),
                                      samplingRate=5.04)
        self.launchProtocol(protImport)
        inputAvgs = protImport.outputAverages
        print(magentaStr("\n==> Testing relion - center averages:"))
        protCenter = self.newProtocol(ProtRelionCenterAverages)
        protCenter.inputAverages.set(inputAvgs)
        self.launchProtocol(protCenter)

        conditions = ['outputAverages.getSize()==%d' % inputAvgs.getSize(),
                      'outputAverages.getSamplingRate() - %0.5f < 0.00001'
                      % inputAvgs.getSamplingRate()]
        self.checkOutput(protCenter, 'outputAverages', conditions)


class TestRelionExportParticles(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.ds2 = DataSet.getDataSet('xmipp_tutorial')
        cls.particlesFn = cls.ds.getFile('import/refine3d/extra/relion_data.star')
        cls.particlesFn2 = cls.ds2.getFile('particles')
        cls.starImport = cls.runImportParticlesStar(cls.particlesFn, 10000, 7.08)
        cls.runImportParticles(cls.particlesFn2, 1.237, True)

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False,
                           phaseFlip=False):
        """ Run an Import particles protocol. """
        print(magentaStr("\n==> Importing data - particles (from xmipp):"))
        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         filesPath=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack,
                                         haveDataBeenPhaseFlipped=phaseFlip)
        cls.launchProtocol(cls.protImport)
        cls.assertIsNotNone(cls.protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return cls.protImport

    def run_combinations(self, inputProt, name=''):
        """ Run an Import particles protocol. """
        inputParts = inputProt.outputParticles

        print(magentaStr("\n==> Testing relion - export particles:"))

        def _checkProt(prot, stackType):
            stackFiles = glob(prot._getExportPath('Particles', '*mrcs'))
            print("stackFiles: ", stackFiles)

            n = len(stackFiles)
            if stackType == 0:
                self.assertEqual(n, 0)
            elif stackType == 1:
                self.assertGreaterEqual(n, 1)
            else:
                self.assertEqual(n, 1)

        stackTypes = [0, 1, 2]
        stackNames = ['no', 'multi', 'single']
        alignments = [True, False]
        combinations = [(s, a) for s in stackTypes for a in alignments]

        for s, a in combinations:
            label = 'export %s (stack: %s - align: %s)' % (name, stackNames[s], a)
            exportProt = self.newProtocol(ProtRelionExportParticles,
                                          inputParticles=inputParts,
                                          objLabel=label,
                                          stackType=s, alignmentType=a)
            self.launchProtocol(exportProt)
            _checkProt(exportProt, s)

    def test_basic(self):
        self.run_combinations(self.starImport)

    def test_extra(self):
        self.run_combinations(self.protImport)


class TestRelionExportCtf(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsGrigorieff = DataSet.getDataSet('grigorieff')
        cls.dsEman2 = DataSet.getDataSet('eman')

        print(magentaStr("\n==> Importing data - micrographs:"))
        cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                         filesPath=cls.dsXmipp.getFile('allMics'),
                                         samplingRate=1.237, voltage=300)
        cls.launchProtocol(cls.protImport)

    def runImportXmipp(self):
        print(magentaStr("\n==> Importing data - ctfs (from xmipp)"))
        protCTF = self.newProtocol(ProtImportCTF,
                                   importFrom=ProtImportCTF.IMPORT_FROM_XMIPP3,
                                   filesPath=self.dsXmipp.getFile('ctfsDir'),
                                   filesPattern='*.ctfparam')
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.setObjLabel('import ctfs from xmipp ')
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")
        return protCTF

    def runImportCtffind4(self):
        print(magentaStr("\n==> Importing data - ctfs (from ctffind)"))
        protCTF = self.newProtocol(ProtImportCTF,
                                   importFrom=ProtImportCTF.IMPORT_FROM_GRIGORIEFF,
                                   filesPath=self.dsGrigorieff.getFile('ctffind4'),
                                   filesPattern='BPV*/*txt')
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.setObjLabel('import from ctffind4')
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")
        return protCTF

    def runImportScipion(self):
        print(magentaStr("\n==> Importing data - ctfs (from scipion)"))
        ctfSqlite = self.dsGrigorieff.getFile('ctffind3/ctfs.sqlite')

        protCTF = self.newProtocol(ProtImportCTF,
                                   objLabel='import from scipion',
                                   importFrom=ProtImportCTF.IMPORT_FROM_SCIPION,
                                   filesPath=ctfSqlite)

        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")
        return protCTF

    def runImportEman2(self):
        print(magentaStr("\n==> Importing data - ctfs (from eman2)"))
        protCTF = self.newProtocol(ProtImportCTF,
                                   importFrom=ProtImportCTF.IMPORT_FROM_EMAN2,
                                   filesPath=self.dsEman2.getFile('ctfs'),
                                   filesPattern='BPV*json')
        protCTF.inputMicrographs.set(self.protImport.outputMicrographs)
        protCTF.setObjLabel('import from eman2')
        self.launchProtocol(protCTF)

        self.assertIsNotNone(protCTF.outputCTF,
                             "There was a problem when importing ctfs.")
        return protCTF

    def testExportCtf(self):
        ctfs = [(self.runImportXmipp().outputCTF, 'xmipp'),
                (self.runImportCtffind4().outputCTF, 'ctffind'),
                (self.runImportScipion().outputCTF, 'scipion'),
                (self.runImportEman2().outputCTF, 'eman2')]

        for i in ctfs:
            protExport = self.newProtocol(ProtRelionExportCtf)
            protExport.inputCTF.set(i[0])
            print(magentaStr("\n==> Testing relion - export ctf (from %s)" % i[1]))
            self.launchProtocol(protExport)

            outFn = os.path.exists(protExport._getStarFile()) or None
            self.assertIsNotNone(outFn,
                                 "There was a problem when exporting ctfs.")


class TestRelionRemovePrefViews(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.particlesFn = cls.ds.getFile('import/refine3d/extra/relion_data.star')
        cls.starImport = cls.runImportParticlesStar(cls.particlesFn, 10000, 7.08)

    def test_removePrefViews(self):
        print(magentaStr("\n==> Testing relion - remove preferential views"))
        inputParts = self.starImport.outputParticles
        prot = self.newProtocol(ProtRelionRemovePrefViews,
                                inputParticles=inputParts,
                                numToRemove=50)
        self.launchProtocol(prot)

        self.assertIsNotNone(prot.outputParticles,
                             "There was a problem with remove preferential views protocol.")
        outSize = prot.outputParticles.getSize()
        self.assertEqual(outSize, 4080, "Output size is not 4080!")


class TestRelionResizeVolume(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')

    def importVolume(self):
        print(magentaStr("\n==> Importing data - volume:"))
        volFn = self.ds.getFile('import/refine3d/extra/relion_class001.mrc')
        protVol = self.newProtocol(ProtImportVolumes,
                                   objLabel='import volume',
                                   filesPath=volFn,
                                   samplingRate=3)
        self.launchProtocol(protVol)
        return protVol

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
        importProt = self.importVolume()
        resizeProt = self.newProtocol(ProtRelionResizeVolume,
                                      doRescale=True, rescaleSamplingRate=1.5,
                                      doResize=True, resizeSize=128)
        vol = importProt.outputVolume
        resizeProt.inputVolumes.set(vol)
        print(magentaStr("\n==> Testing relion - crop/resize volumes:"))
        self.launchProtocol(resizeProt)

        self._validations(resizeProt.outputVol, 128, 1.5)
