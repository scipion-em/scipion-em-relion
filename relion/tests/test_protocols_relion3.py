# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.tests as pwtests
from pyworkflow.utils import copyTree, magentaStr
import pwem.protocols as emprot
from pwem.tests.workflows import TestWorkflow
from pwem.objects import SetOfMovies

import relion
import relion.convert
from relion.convert.convert31 import OpticsGroups
from relion.protocols import ProtRelionMotioncor, ProtRelionAssignOpticsGroup

CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class Relion3TestProtocolBase(TestWorkflow):
    GROUP_NAME = "opticsGroupTest"
    MTF_FILE = os.path.join(os.path.dirname(relion.convert.__file__), 'mtfs',
                            'mtf_k2_300_ec.star')

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')

    @classmethod
    def _importMovies(cls, **kwargs):
        print(magentaStr("\n==> Importing data - movies:"))
        protImport = cls.newProtocol(
            emprot.ProtImportMovies,
            filesPath=cls.ds.getFile('Movies/'),
            filesPattern=kwargs.get('filesPattern', '20170629_000?5*tiff'),
            samplingRateMode=0,
            samplingRate=0.885,
            magnification=50000,
            scannedPixelSize=7.0,
            voltage=200,
            sphericalAberration=1.4,
            doseInitial=0.0,
            dosePerFrame=1.277,
            gainFile=cls.ds.getFile("Movies/gain.mrc")
        )
        protImport.setObjLabel('import 3 movies')
        protImport.setObjComment('Relion 3 tutorial movies:\n\n'
                                 'Microscope Jeol Cryo-ARM 200\n'
                                 'Data courtesy of Takyuki Kato in the Namba '
                                 'group\n(Osaka University, Japan)')
        protImport = cls.launchProtocol(protImport)

        return protImport

    def _runRelionMc(self, protImport, **kwargs):
        # if not relion.Plugin.IS_30():
        #     protInput = self._runAssignOptics(protImport)
        # else:
        #     protInput = protImport

        protInput = protImport

        args = {
            'objLabel': 'relion - motioncor',
            'patchX': 5,
            'patchY': 5,
            'numberOfThreads': CPUS
        }
        args.update(kwargs)
        print(magentaStr("\nRunning relion - motioncor:"))
        protRelionMc = self.newProtocol(ProtRelionMotioncor, **args)
        protRelionMc.inputMovies.set(protInput.outputMovies)
        protRelionMc = self.launchProtocol(protRelionMc)

        return protRelionMc

    def _runAssignOptics(self, protInput, **kwargs):

        args = {
            'objLabel': 'relion - assign optics',
            'opticsGroupName': self.GROUP_NAME,
            'mtfFile': self.MTF_FILE,
        }
        args.update(kwargs)

        print(magentaStr("\nRunning relion - assign optics groups:"))
        protAssign = self.newProtocol(ProtRelionAssignOpticsGroup, **args)
        protAssign.inputSet.set(protInput.outputMovies)
        return self.launchProtocol(protAssign)


class Relion3TestAssignOptics(Relion3TestProtocolBase):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')
        # Run import only once and with 3 movies
        cls.protImport = cls._importMovies(filesPattern='20170629_000?5*tiff')

    def _checkOutputMovies(self, prot, size, exists=True,
                           hasAlignment=True):
        # Validate output movies
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
            self.assertEqual(fog.rlnOpticsGroupName, self.GROUP_NAME)
            self.assertEqual(os.path.basename(fog.rlnMtfFileName),
                             os.path.basename(self.MTF_FILE))

        print(magentaStr("\n==> Testing relion - assign optics groups:"))
        protRelionAssign = self._runAssignOptics(self.protImport)
        self._checkOutputMovies(protRelionAssign, 3, hasAlignment=False)
        output = protRelionAssign.outputMovies
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
data_micrographs

loop_ 
_rlnMicrographName #1 
_rlnOpticsGroup #2 
20170629_00025_frameImage.tiff 1
20170629_00035_frameImage.tiff 2
20170629_00045_frameImage.tiff 3
        """)
        f.close()

        print(magentaStr("\nRunning relion - assign optics groups:"))
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


class Relion3TestMotioncor(Relion3TestProtocolBase):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')
        # Run import only once and with 3 movies
        cls.protImport = cls._importMovies(filesPattern='20170629_000?5*tiff')

    def _checkOutputMovies(self, prot, size, exists=True,
                           hasAlignment=True):
        # Validate output movies
        movies = getattr(prot, 'outputMovies', None)
        self.assertIsNotNone(movies, "No movies were generated")
        # dims = movies.getDim()
        # cls.assertEqual((3710, 3838, 24), dims)
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


class Relion3TestMultiBody(Relion3TestProtocolBase):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')
        cls.extra = cls.ds.getFile("multibody/extra")
        cls.ref3d = cls.ds.getFile("multibody/ref3d")
        cls.ptcls = cls.ds.getFile("multibody/ref3d/relion_it017_data.star")

    def _setupRefinement(self):
        from pyworkflow.protocol.constants import STATUS_FINISHED
        relionRefine = self.newProtocol(relion.protocols.ProtRelionRefine3D,
                                        objLabel='fake 3D refinement',
                                        referenceMask=None)
        self.saveProtocol(relionRefine)
        relionRefine.setStatus(STATUS_FINISHED)

        # copy ref3d files into protocol dir
        currDir2 = os.path.join(self.proj.getPath(), relionRefine._getExtraPath())
        print("Copying files from %s to %s" % (self.ref3d, currDir2))
        copyTree(self.ref3d, currDir2)

        return relionRefine

    def testMultibody(self):
        print(magentaStr("\n==> Testing relion - multi-body:"))
        relionMbody = self.newProtocol(relion.protocols.ProtRelionMultiBody,
                                       initialOffsetRange=2.0,
                                       initialOffsetStep=0.5,
                                       runFlexAnalysis=False,
                                       pooledParticles=30,
                                       skipPadding=True,
                                       doGpu=True,
                                       gpusToUse='0,1:2,3',
                                       numberOfThreads=12,
                                       numberOfMpis=3)
        protRef = self._setupRefinement()
        relionMbody.protRefine.set(protRef)

        # copy m-body files into protocol dir
        currDir1 = os.path.join(self.proj.getPath(), relionMbody._getPath("Uploads"))
        print("Copying files from %s to %s" % (self.extra, currDir1))
        copyTree(self.extra, currDir1)

        bodyFn = os.path.join(self.proj.getPath(), relionMbody._getPath('Uploads/2-bodies.star'))
        relionMbody.bodyStarFile.set(bodyFn)

        self.saveProtocol(relionMbody)
        self.launchProtocol(relionMbody)
        self.assertIsNotNone(relionMbody.outputVolumes,
                             "There was a problem with Relion multi-body")


class TestRelion31ImportParticles(pwtests.BaseTest):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion31_tutorial_precalculated')

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

        prot1 = self.newProtocol(emprot.ProtImportParticles,
                                 objLabel='from relion (extract job)',
                                 importFrom=emprot.ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=starFile,
                                 magnification=10000,
                                 samplingRate=optics.rlnImagePixelSize,
                                 haveDataBeenPhaseFlipped=False
                                 )
        self.launchProtocol(prot1)
        self.checkOutput(prot1, 'outputParticles', [])

    def test_fromClassify2D(self):
        """ Import particles from Classify 2d job star file.
        """
        starFile = self.ds.getFile('Class2D/job013/run_it025_data.star')
        optics = OpticsGroups.fromStar(starFile).first()

        prot1 = self.newProtocol(emprot.ProtImportParticles,
                                 objLabel='from relion (classify 2d)',
                                 importFrom=emprot.ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=starFile,
                                 magnification=120000,
                                 samplingRate=optics.rlnImagePixelSize,
                                 haveDataBeenPhaseFlipped=False
                                 )
        self.launchProtocol(prot1)
        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignment2D()'])
        self.checkOutput(prot1, 'outputClasses')

    def test_fromRefine3D(self):
        """ Import particles from Refine3D job star file. """
        starFile = self.ds.getFile('Refine3D/job019/run_it020_data.star')
        optics = OpticsGroups.fromStar(starFile).first()

        prot1 = self.newProtocol(emprot.ProtImportParticles,
                                 objLabel='from relion (refine 3d)',
                                 importFrom=emprot.ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=starFile,
                                 magnification=10000,
                                 samplingRate=optics.rlnImagePixelSize,
                                 haveDataBeenPhaseFlipped=False
                                 )
        self.launchProtocol(prot1)
        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignmentProj()'])

    def test_fromClassify3D(self):
        """ Import particles from Classify3D job star file. """
        starFile = self.ds.getFile('Class3D/job016/run_it025_data.star')
        optics = OpticsGroups.fromStar(starFile).first()

        prot1 = self.newProtocol(emprot.ProtImportParticles,
                                 objLabel='from relion (classify 3d)',
                                 importFrom=emprot.ProtImportParticles.IMPORT_FROM_RELION,
                                 starFile=starFile,
                                 magnification=120000,
                                 samplingRate=optics.rlnImagePixelSize,
                                 haveDataBeenPhaseFlipped=False
                                 )
        self.launchProtocol(prot1)
        self.checkOutput(prot1, 'outputParticles', ['outputParticles.hasAlignmentProj()'])
        self.checkOutput(prot1, 'outputClasses')
