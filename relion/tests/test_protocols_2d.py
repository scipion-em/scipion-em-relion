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

from glob import glob

from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.plugin import Domain
from pyworkflow.utils import magentaStr
from pwem.objects import SetOfMovies
from pwem.protocols import (ProtImportAverages, ProtImportCTF,
                            ProtImportParticles, ProtImportCoordinates)

from ..protocols import *
from ..convert import *
from ..constants import *
from relion.convert.convert31 import OpticsGroups
from .test_protocols_base import TestRelionBase, USE_GPU, RUN_CPU, CPUS, MTF_FILE


class TestRelionAssignOptics(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion30_tutorial')
        cls.protImport = cls.runImportMovies(filesPath=cls.ds.getFile('Movies/'),
                                             filesPattern='20170629_000?5*tiff',
                                             samplingRate=0.885,
                                             voltage=200,
                                             sphericalAberration=1.4,
                                             dose=1.277,
                                             gain=cls.ds.getFile("Movies/gain.mrc"))

    def _checkOutputMovies(self, prot, size, exists=True,
                           hasAlignment=True):
        movies = getattr(prot, 'outputMovies', None)
        self.assertIsNotNone(movies, "No movies were generated")
        self.assertEqual(size, movies.getSize())

        if hasAlignment:
            self.assertTrue(movies.getFirstItem().hasAlignment())

        if exists:
            for m in movies:
                self.assertTrue(os.path.exists(m.getFileName()))

    def test_single(self):
        def _checkAcq(output):
            og = OpticsGroups.fromImages(output)
            fog = og.first()
            self.assertEqual(fog.rlnOpticsGroupName, "opticsGroupTest")
            self.assertEqual(os.path.basename(fog.rlnMtfFileName),
                             os.path.basename(MTF_FILE))

        print(pwutils.magentaStr("\n==> Testing relion - assign optics groups:"))
        protAssign = self.newProtocol(ProtRelionAssignOpticsGroup,
                                      inputSet=self.protImport.outputMovies,
                                      inputType=0,
                                      opticsGroupName="opticsGroupTest",
                                      mtfFile=MTF_FILE)
        protAssign = self.launchProtocol(protAssign)

        output = protAssign.outputMovies
        self._checkOutputMovies(protAssign, 3, hasAlignment=False)
        _checkAcq(output)
        output.close()

        print("Loading db: %s" % os.path.abspath(output.getFileName()))
        moviesSet = SetOfMovies(filename=output.getFileName())
        moviesSet.loadAllProperties()

    def test_fromStar(self):
        mtfStar = self.getOutputPath('input_mtf.star')
        f = open(mtfStar, 'w')
        f.write("""
data_optics

loop_ 
_rlnOpticsGroupName #1 
_rlnOpticsGroup #2 
_rlnMtfFileName #3 
_rlnMicrographOriginalPixelSize #4 
_rlnVoltage #5 
_rlnSphericalAberration #6 
_rlnAmplitudeContrast #7 
_rlnMicrographPixelSize #8 
opticsGroup1            1 mtf_k2_200kV.star     0.885000   200.000000     1.400000     0.100000     0.885000
opticsGroup2            2 mtf_k2_200kV.star     0.885000   200.000000     1.400000     0.100000     0.885000
opticsGroup3            3 mtf_k2_200kV.star     0.885000   200.000000     1.400000     0.100000     0.885000

data_micrographs

loop_
_rlnMicrographName #1 
_rlnOpticsGroup #2 
20170629_00025_frameImage.tiff 1
20170629_00035_frameImage.tiff 2
20170629_00045_frameImage.tiff 3
        """)
        f.close()

        print(pwutils.magentaStr("\n==> Testing relion - assign optics groups (from star):"))
        protAssign = self.newProtocol(ProtRelionAssignOpticsGroup,
                                      objLabel='assign optics - from star',
                                      inputSet=self.protImport.outputMovies,
                                      inputType=1,  # from star file
                                      inputStar=mtfStar)
        protAssign = self.launchProtocol(protAssign)

        outputMovies = protAssign.outputMovies
        og = OpticsGroups.fromImages(outputMovies)
        self.assertEqual(3, len(og))

        for i, movie in enumerate(outputMovies):
            self.assertEqual(i + 1, movie.getAttributeValue('_rlnOpticsGroup'))


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
                                      numberOfMpi=1, numberOfThreads=8)
            prot2D.numberOfClasses.set(4)
            prot2D.numberOfVDAMBatches.set(100)
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


class TestRelionExportCtf(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsXmipp = DataSet.getDataSet('xmipp_tutorial')
        cls.dsGrigorieff = DataSet.getDataSet('grigorieff')
        cls.dsEman2 = DataSet.getDataSet('eman')
        cls.protImport = cls.runImportMics(cls.dsXmipp.getFile('allMics'), 1.237)

    def runImportXmipp(self):
        print(magentaStr("\n==> Importing data - ctfs (from xmipp):"))
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
        print(magentaStr("\n==> Importing data - ctfs (from ctffind):"))
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
        print(magentaStr("\n==> Importing data - ctfs (from scipion):"))
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
        print(magentaStr("\n==> Importing data - ctfs (from eman2):"))
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
            print(magentaStr("\n==> Testing relion - export ctf (from %s):" % i[1]))
            self.launchProtocol(protExport)

            outFn = os.path.exists(protExport._getStarFile()) or None
            self.assertIsNotNone(outFn,
                                 "There was a problem when exporting ctfs.")


class TestRelionExportParticles(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.ds2 = DataSet.getDataSet('xmipp_tutorial')
        cls.particlesFn = cls.ds.getFile('import/refine3d/extra/relion_data.star')
        cls.particlesFn2 = cls.ds2.getFile('particles')
        cls.starImport = cls.runImportParticlesStar(cls.particlesFn, 7.08)
        cls.protImport = cls.runImportParticles(cls.particlesFn2, 1.237, checkStack=True)

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


class TestRelionExtractParticles(TestRelionBase):
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
    def runFakePicking(cls, mics, pattern):
        """ Run a fake particle picking. Coordinates already exist. """

        # TODO This fake picking depends on Xmipp particle picking
        # TODO: Can we change this to an import coordinates?

        XmippProtParticlePicking = Domain.importFromPlugin(
            'xmipp3.protocols', 'XmippProtParticlePicking')

        if XmippProtParticlePicking is None:
            print("WARNING: Can not load xmipp3.protocols.XmippProtParticlePicking."
                  "Skipping the test.")
            return None

        print(magentaStr("\n==> Running xmipp - fake particle picking:"))
        cls.protPP = XmippProtParticlePicking(importFolder=pattern, runMode=1)
        cls.protPP.inputMicrographs.set(mics)
        cls.proj.launchProtocol(cls.protPP, wait=True)
        # check that fake picking has run ok
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
        cls.protImport = cls.runImportMics(cls.micsFn, 1.237)
        cls.protDown = cls.runDownsamplingMicrographs(cls.protImport.outputMicrographs,
                                                      downFactorValue=5.0)

        print(magentaStr("\n==> Importing data - ctfs:"))
        cls.protCTF = cls.newProtocol(ProtImportCTF,
                                      importFrom=ProtImportCTF.IMPORT_FROM_XMIPP3,
                                      filesPath=cls.dataset.getFile('ctfsDir'),
                                      filesPattern='*.ctfparam')
        cls.protCTF.inputMicrographs.set(cls.protImport.outputMicrographs)
        cls.proj.launchProtocol(cls.protCTF, wait=True)

        cls.protPP = cls.runFakePicking(cls.protDown.outputMicrographs, cls.allCrdsDir)

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


class TestRelionImportCoords(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion31_tutorial_precalculated')
        cls.mics = cls.ds.getFile('MotionCorr/job002/Movies/*frameImage.mrc')
        cls.parts = cls.ds.getFile('Refine3D/job029/run_it018_data.star')
        cls.coords = cls.ds.getFile('AutoPick/job006/Movies')
        cls.protImportMics = cls.runImportMics(cls.mics, 0.885)

    def testImportCoords(self):
        """ Run an Import coords protocol from a particle star file. """
        print(magentaStr("\n==> Importing coordinates (from particles star file):"))
        self.protImport = self.newProtocol(ProtImportCoordinates,
                                           objLabel="from particles star file",
                                           importFrom=2,  # RELION
                                           filesPath=self.parts,
                                           boxSize=128)
        self.protImport.inputMicrographs.set(self.protImportMics.outputMicrographs)
        self.launchProtocol(self.protImport)
        self.assertIsNotNone(self.protImport.outputCoordinates,
                             "SetOfCoordinates has not been produced.")
        self.assertEqual(self.protImport.outputCoordinates.getSize(),
                         4501, "Output size is not 4501!")

    def testImportCoords2(self):
        """ Run an Import coords protocol. """
        print(magentaStr("\n==> Importing coordinates (from coordinates star file):"))
        self.protImport = self.newProtocol(ProtImportCoordinates,
                                           objLabel="from coordinates star file",
                                           importFrom=2,  # RELION
                                           filesPath=self.coords,
                                           filesPattern="*autopick.star",
                                           boxSize=128)
        self.protImport.inputMicrographs.set(self.protImportMics.outputMicrographs)
        self.launchProtocol(self.protImport)
        self.assertIsNotNone(self.protImport.outputCoordinates,
                             "SetOfCoordinates has not been produced.")
        self.assertEqual(self.protImport.outputCoordinates.getSize(),
                         1158, "Output size is not 1158!")


class TestRelionImportParticles(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion31_tutorial_precalculated')

    def checkOutput(self, prot, outputName, conditions=[]):
        """ Check that an ouput was generated and
        the condition is valid.
        """
        o = getattr(prot, outputName, None)
        locals()[outputName] = o
        self.assertIsNotNone(o, "Output: %s is None" % outputName)
        for cond in conditions:
            self.assertTrue(eval(cond), 'Condition failed: ' + cond)

    def test_fromExtract(self):
        """ Import particles.star from Extract job. """
        starFile = self.ds.getFile('Extract/job018/particles.star')
        optics = OpticsGroups.fromStar(starFile).first()
        prot1 = self.runImportParticlesStar(starFile, samplingRate=optics.rlnImagePixelSize,
                                            label="from relion (extract job)")

        self.checkOutput(prot1, 'outputParticles', [])

    def test_fromClassify2D(self):
        """ Import particles from Classify 2d job star file.
        """
        starFile = self.ds.getFile('Class2D/job013/run_it025_data.star')
        optics = OpticsGroups.fromStar(starFile).first()
        prot1 = self.runImportParticlesStar(starFile, samplingRate=optics.rlnImagePixelSize,
                                            label="from relion (classify 2d)")

        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignment2D()'])
        self.checkOutput(prot1, 'outputClasses')

    def test_fromRefine3D(self):
        """ Import particles from Refine3D job star file. """
        starFile = self.ds.getFile('Refine3D/job019/run_it020_data.star')
        optics = OpticsGroups.fromStar(starFile).first()
        prot1 = self.runImportParticlesStar(starFile, samplingRate=optics.rlnImagePixelSize,
                                            label="from relion (refine 3d)")

        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignmentProj()'])

    def test_fromClassify3D(self):
        """ Import particles from Classify3D job star file. """
        starFile = self.ds.getFile('Class3D/job016/run_it025_data.star')
        optics = OpticsGroups.fromStar(starFile).first()
        prot1 = self.runImportParticlesStar(starFile, samplingRate=optics.rlnImagePixelSize,
                                            label="from relion (classify 3d)")

        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignmentProj()'])
        self.checkOutput(prot1, 'outputClasses')


class TestRelionMotioncor(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion30_tutorial')
        cls.protImport = cls.runImportMovies(filesPath=cls.ds.getFile('Movies/'),
                                             filesPattern='20170629_000?5*tiff',
                                             samplingRate=0.885,
                                             voltage=200,
                                             sphericalAberration=1.4,
                                             dose=1.277,
                                             gain=cls.ds.getFile("Movies/gain.mrc"))

    def _runRelionMc(self, protImport, **kwargs):
        protInput = protImport

        args = {
            'objLabel': 'relion - motioncor',
            'patchX': 5,
            'patchY': 5,
            'numberOfThreads': CPUS
        }
        args.update(kwargs)
        protRelionMc = self.newProtocol(ProtRelionMotioncor, **args)
        protRelionMc.inputMovies.set(protInput.outputMovies)
        protRelionMc = self.launchProtocol(protRelionMc)

        return protRelionMc

    def _checkOutputMovies(self, prot, size, exists=True,
                           hasAlignment=True):
        movies = getattr(prot, 'outputMovies', None)
        self.assertIsNotNone(movies, "No movies were generated")
        self.assertEqual(size, movies.getSize())

        if hasAlignment:
            self.assertTrue(movies.getFirstItem().hasAlignment())

        if exists:
            for m in movies:
                self.assertTrue(os.path.exists(m.getFileName()))

    def test_1x1(self):
        print(magentaStr("\n==> Testing relion - motioncor (global):"))
        protRelionMc = self._runRelionMc(self.protImport, objLabel='relion - mc 1x1',
                                         patchX=1, patchY=1)
        self._checkOutputMovies(protRelionMc, 3)

    def test_1x1_PS(self):
        print(magentaStr("\n==> Testing relion - motioncor (global + PS):"))
        protRelionMc = self._runRelionMc(self.protImport, objLabel='relion - mc PS',
                                         patchX=1, patchY=1,
                                         savePSsum=True)
        self._checkOutputMovies(protRelionMc, 3)

    def test_3x3_DW(self):
        print(magentaStr("\n==> Testing relion - motioncor (local + DW):"))
        protRelionMc = self._runRelionMc(self.protImport, objLabel='relion - mc 3x3 DW',
                                         patchX=3, patchY=3, doDW=True)
        self._checkOutputMovies(protRelionMc, 3)


class TestRelionPicking(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.ds.getFile('import/classify2d/extra/relion_it015_data.star')
        cls.protImportMics = cls.runImportMics('%s/*.mrc' % cls.ds.getFile('micrographs'),
                                               samplingRate=7.08)

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
            streamingBatchSize=5
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
            resamplePix=False
        )
        self.launchProtocol(protCtf)

        print(magentaStr("\n==> Testing relion - autopick ref:"))
        protPickRef = self.newProtocol(
            ProtRelion2Autopick,
            inputMicrographs=self.protImportMics.outputMicrographs,
            ctfRelations=protCtf.outputCTF,
            inputReferences=protAvg.outputAverages,
            streamingBatchSize=5
        )
        self.launchProtocol(protPickRef)
        self._checkOutput(protPickRef, 320, 370)


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


class TestRelionRemovePrefViews(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('relion_tutorial')
        cls.particlesFn = cls.ds.getFile('import/refine3d/extra/relion_data.star')
        cls.starImport = cls.runImportParticlesStar(cls.particlesFn, 7.08)

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
