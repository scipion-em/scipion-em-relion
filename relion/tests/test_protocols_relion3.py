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

import relion
from relion.protocols import ProtRelionMotioncor, ProtRelionAssignOpticsGroup

CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class Relion3TestProtocolBase(TestWorkflow):
    GROUP_NAME = "opticsGroupTest"
    MTF_FILE = "/test/mtf.star"

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')

    @classmethod
    def _importMovies(cls, **kwargs):
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
        if not relion.IS_30:
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
        # # Validate output movies
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

    def test_assign(self):
        if relion.IS_30:
            print("This test only make sense for Relion >= 3.1. Exiting...")
            return

        def _checkAcq(obj):
            acq = obj.getAcquisition()
            self.assertEqual(acq.getAttributeValue('opticsGroupName', ''),
                             self.GROUP_NAME)
            self.assertEqual(acq.getAttributeValue('mtfFile', ''),
                             self.MTF_FILE)

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
        # # Validate output movies
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
        protRelionMc = self._runRelionMc(self.protImport, patchX=1, patchY=1)
        self._checkOutputMovies(protRelionMc, 3)

    def test_1x1_PS(self):
        protRelionMc = self._runRelionMc(self.protImport,
                                         patchX=1, patchY=1,
                                         )
        self._checkOutputMovies(protRelionMc, 3)

    def test_2x2(self):
        protRelionMc = self._runRelionMc(self.protImport, patchX=2, patchY=2)
        self._checkOutputMovies(protRelionMc, 3)
