# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
from pyworkflow.tests.em.workflows import TestWorkflow
import pyworkflow.em as pwem

import relion
import relion.protocols


CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class Relion3TestProtocolBase(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')

    @classmethod
    def _importMovies(cls, **kwargs):
        protImport = cls.newProtocol(
            pwem.ProtImportMovies,
            filesPath=cls.ds.getFile('Movies/'),
            filesPattern=kwargs.get('filesPattern','*.tiff'),
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

    def _runRelionMc(self, protImport, patchX=5, patchY=5):
        protRelionMc = self.newProtocol(
            relion.protocols.ProtRelionMotioncor,
            objLabel='relion - motioncor',
            patchX=patchX, patchY=patchY,
            numberOfThreads=CPUS,
        )

        protRelionMc.inputMovies.set(protImport.outputMovies)
        protRelionMc = self.launchProtocol(protRelionMc)

        return protRelionMc


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

    def test_2x2(self):
        protRelionMc = self._runRelionMc(self.protImport, patchX=2, patchY=2)
        self._checkOutputMovies(protRelionMc, 3)


if __name__ == '__main__':
    import unittest
    unittest.main()
