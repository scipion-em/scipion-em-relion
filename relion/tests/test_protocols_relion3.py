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
from pwem.tests.workflows import TestWorkflow
from pwem.protocols import ProtImportMovies
from pwem.objects import SetOfMovies
from pyworkflow.utils import copyTree, join, magentaStr

import relion
import relion.convert
from relion.protocols import ProtRelionMotioncor, ProtRelionAssignOpticsGroup

CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class Relion3TestProtocolBase(TestWorkflow):
    GROUP_NAME = "opticsGroupTest"
    MTF_FILE = join(os.path.dirname(relion.convert.__file__), 'mtfs',
                    'mtf_k3_standard_200kV.star')

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')

    @classmethod
    def _importMovies(cls, **kwargs):
        print(magentaStr("\n==> Importing data - movies:"))
        protImport = cls.newProtocol(
            ProtImportMovies,
            filesPath=cls.ds.getFile('Movies/'),
            filesPattern=kwargs.get('filesPattern', '*.tiff'),
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
        if not relion.Plugin.IS_30():
            protInput = self._runAssignOptics(protImport)
        else:
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

    def test_assign(self):
        if relion.Plugin.IS_30():
            print("This test only makes sense for Relion >= 3.1. Exiting...")
            return

        def _checkAcq(obj):
            acq = obj.getAcquisition()
            self.assertEqual(acq.getAttributeValue('opticsGroupName', ''),
                             self.GROUP_NAME)
            self.assertEqual(acq.getAttributeValue('mtfFile', ''),
                             self.MTF_FILE)

        print(magentaStr("\n==> Testing relion - assign optics groups:"))
        protRelionAssign = self._runAssignOptics(self.protImport)
        self._checkOutputMovies(protRelionAssign, 3, hasAlignment=False)
        output = protRelionAssign.outputMovies
        _checkAcq(output)
        _checkAcq(output.getFirstItem())

        output.close()

        print("Loading db: %s" % os.path.abspath(output.getFileName()))
        moviesSet = SetOfMovies(filename=output.getFileName())
        moviesSet.loadAllProperties()


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
        if relion.Plugin.IS_GT30():
            print(magentaStr("\n==> Testing relion - motioncor (global + PS):"))
            protRelionMc = self._runRelionMc(self.protImport, objLabel='relion - mc PS',
                                             patchX=1, patchY=1,
                                             savePSsum=True)
            self._checkOutputMovies(protRelionMc, 3)
        else:
            print("Cannot test motioncorr with PS saving - it's only available for Relion 3.1+")

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
