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
from pyworkflow.plugin import Domain
from pyworkflow.utils import magentaStr
from pwem.tests.workflows import TestWorkflow
from pwem.protocols import ProtImportMovies

from ..protocols import *


CPUS = os.environ.get('SCIPION_TEST_CPUS', 8)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)


class TestWorkflowRelionBetagal(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('relion30_tutorial')

    def _importMovies(self):
        print(magentaStr("\n==> Importing data - movies:"))
        protImport = self.newProtocol(
            ProtImportMovies,
            filesPath=self.ds.getFile('Movies/'),
            filesPattern='*.tiff',
            samplingRateMode=0,
            samplingRate=0.885,
            magnification=50000,
            scannedPixelSize=7.0,
            voltage=200,
            sphericalAberration=1.4,
            doseInitial=0.0,
            dosePerFrame=1.277,
            gainFile=self.ds.getFile("Movies/gain.mrc")
        )
        protImport.setObjLabel('import 24 movies')
        protImport.setObjComment('Relion 3 tutorial movies:\n\n'
                                 'Microscope Jeol Cryo-ARM 200\n'
                                 'Data courtesy of Takyuki Kato in the Namba '
                                 'group\n(Osaka University, Japan)')
        protImport = self.launchProtocol(protImport)

        # Validate output movies
        movies = getattr(protImport, 'outputMovies', None)
        self.assertIsNotNone(movies, "No movies were generated from the import")
        dims = movies.getDim()
        self.assertEqual((3710, 3838, 24), dims)
        self.assertEqual(24, movies.getSize())

        print(magentaStr("\n==> Testing relion - assign optic groups:"))
        protAssign = self.newProtocol(ProtRelionAssignOpticsGroup,
                                      objLabel='assign optics',
                                      opticsGroupName='OpticsGroup1')
        protAssign.inputSet.set(protImport.outputMovies)
        return self.launchProtocol(protAssign)

    def _runRelionMc(self, protImport):
        print(magentaStr("\n==> Testing relion - motioncor:"))
        protRelionMc = self.newProtocol(
            ProtRelionMotioncor,
            objLabel='relion - motioncor',
            patchX=1, patchY=1,
            saveFloat16=False,
            numberOfMpi=CPUS//2)

        protRelionMc.inputMovies.set(protImport.outputMovies)
        protRelionMc = self.launchProtocol(protRelionMc)

        return protRelionMc

    def _runRelionLog(self, protRelionMc):
        print(magentaStr("\n==> Testing relion - autopick LoG:"))
        protRelionLog = self.newProtocol(
            ProtRelionAutopickLoG,
            objLabel='relion - autopick log',
            boxSize=100,
            minDiameter=150, maxDiameter=180,
            threshold2=5.0,
            streamingBatchSize=0,
            numberOfMpi=1)

        protRelionLog.inputMicrographs.set(protRelionMc.outputMicrographsDoseWeighted)
        protRelionLog = self.launchProtocol(protRelionLog)

        return protRelionLog

    def _runCtffind(self, protMc):
        print(magentaStr("\n==> Testing ctffind - estimate ctf:"))
        ProtCtfFind = Domain.importFromPlugin('cistem.protocols', 'CistemProtCTFFind')

        protCtf = self.newProtocol(
            ProtCtfFind,
            objLabel='ctffind',
            lowRes=20, highRes=3,
            astigmatism=100,
            windowSize=512,
            usePowerSpectra=True,
            streamingBatchSize=0
        )

        protCtf.inputMicrographs.set(protMc.outputMicrographsDoseWeighted)
        protCtf = self.launchProtocol(protCtf)

        return protCtf

    def _runCtffind(self, protMc):
        """ Run CTFFind protocol. """
        CistemProtCTFFind = Domain.importFromPlugin('cistem.protocols', 'CistemProtCTFFind')
        protCTF = self.newProtocol(CistemProtCTFFind)
        protCTF.inputMicrographs.set(protMc.outputMicrographsDoseWeighted)
        self.launchProtocol(protCTF)
        return protCTF

    def _runRelionExtract(self, protPicking, protCtf):
        print(magentaStr("\n==> Testing relion - extract particles:"))
        protRelionExtract = self.newProtocol(
            ProtRelionExtractParticles,
            objLabel='relion - extract',
            boxSize=320, doRescale=True, rescaledSize=128,
            doInvert=True, doNormalize=True,
            backDiameter=220,
            numberOfMpi=CPUS,
            streamingBatchSize=0,
            downsamplingType=0  # Micrographs same as picking
        )

        protRelionExtract.ctfRelations.set(protCtf.outputCTF)
        protRelionExtract.inputCoordinates.set(protPicking.outputCoordinates)
        protRelionExtract = self.launchProtocol(protRelionExtract)

        return protRelionExtract

    def _runRelion2D(self, protExtract):
        print(magentaStr("\n==> Testing relion - classify 2D:"))
        protRelion2D = self.newProtocol(
            ProtRelionClassify2D,
            objLabel='relion - 2d',
            inplaneAngularSamplingDeg=11,
            maskDiameterA=220,
            numberOfClasses=20,
            extraParams='--maxsig 25',
            pooledParticles=30,
            doGpu=True,
            numberOfThreads=CPUS,
            numberOfMpi=1,
            allParticlesRam=True
        )

        protRelion2D.inputParticles.set(protExtract.outputParticles)
        return self.launchProtocol(protRelion2D)

    def _runRelion2DSelection(self, prot2D):
        print(magentaStr("\n==> Testing relion - 2D class selection:"))
        prot2DSelect = self.newProtocol(
            ProtRelionSelectClasses2D,
            objLabel='relion - select 2d',
            minThreshold=0.1
        )
        prot2DSelect.inputProtocol.set(prot2D)
        return self.launchProtocol(prot2DSelect)

    def _runInitModel(self, protExtract):
        print(magentaStr("\n==> Testing relion - initial model:"))
        relionIniModel = self.newProtocol(
            ProtRelionInitialModel,
            doCTF=True, doGpu=True,
            maskDiameterA=220,
            numberOfIter=100,
            symmetryGroup='C1',
            allParticlesRam=True,
            pooledParticles=30,
            numberOfThreads=CPUS)
        relionIniModel.inputParticles.set(protExtract.outputParticles)
        return self.launchProtocol(relionIniModel)

    def _runRelionSymmetrize(self, inputVol):
        print(magentaStr("\n==> Testing relion - symmetrize 3D:"))
        relionSym = self.newProtocol(
            ProtRelionSymmetrizeVolume,
            symmetryGroup='D2')
        relionSym.inputVolume.set(inputVol)
        return self.launchProtocol(relionSym)

    def _runRelion3D(self, inputVol, inputPts):
        print(magentaStr("\n==> Testing relion - auto-refine 3D:"))
        protRelion3D = self.newProtocol(
            ProtRelionRefine3D,
            doGpu=True,
            maskDiameterA=220,
            symmetryGroup='D2',
            initialLowPassFilterA=30,
            pooledParticles=30,
            allParticlesRam=True,
            numberOfMpi=3,
            numberOfThreads=4)

        protRelion3D.inputParticles.set(inputPts)
        protRelion3D.referenceVolume.set(inputVol)
        protInitModel = self.launchProtocol(protRelion3D)
        return protInitModel

    def test_workflow(self):
        protImport = self._importMovies()
        protRelionMc = self._runRelionMc(protImport)
        protCtfFind = self._runCtffind(protRelionMc)
        protRelionLog = self._runRelionLog(protRelionMc)
        protRelionExtract = self._runRelionExtract(protRelionLog, protCtfFind)
        protRelion2D = self._runRelion2D(protRelionExtract)
        protRelion2DSel = self._runRelion2DSelection(protRelion2D)
        protInitModel = self._runInitModel(protRelionExtract)
        protSym = self._runRelionSymmetrize(protInitModel.outputVolume)
        protRelion3D = self._runRelion3D(protSym.outputVolumeAligned,
                                         protInitModel.outputParticles)
        self.assertIsNotNone(protRelion2D.outputClasses, protRelion2D.getErrorMessage())
        self.assertIsNotNone(protInitModel.outputVolume, protInitModel.getErrorMessage())
