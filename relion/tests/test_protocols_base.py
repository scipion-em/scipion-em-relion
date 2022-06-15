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

from pyworkflow.tests import BaseTest, DataSet
from pyworkflow.utils import magentaStr
from pwem import Config
from pwem.protocols import *


import relion.convert
from ..convert import *
from ..protocols import ProtRelionPreprocessParticles


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
CPUS = os.environ.get('SCIPION_TEST_CPUS', 4)
GPUS = os.environ.get('SCIPION_TEST_GPUS', 2)
MTF_FILE = os.path.join(os.path.dirname(relion.convert.__file__), 'mtfs',
                        'mtf_k2_300_ec.star')


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
    def runImportParticles(cls, pattern, samplingRate, checkStack=False, phaseFlipped=False):
        """ Run an Import particles protocol. """
        print(magentaStr("\n==> Importing data - particles from files:"))
        protImport = cls.newProtocol(ProtImportParticles,
                                     filesPath=pattern,
                                     samplingRate=samplingRate,
                                     checkStack=checkStack,
                                     haveDataBeenPhaseFlipped=phaseFlipped)
        cls.launchProtocol(protImport)
        cls.assertIsNotNone(protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return protImport

    @classmethod
    def runImportParticlesStar(cls, partStar, samplingRate, phaseFlipped=False,
                               label=None):
        """ Import particles from Relion star file. """
        print(magentaStr("\n==> Importing data - particles from star:"))
        protImport = cls.newProtocol(ProtImportParticles,
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=partStar,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=phaseFlipped)
        if label is not None:
            protImport.setObjLabel(label)
        cls.launchProtocol(protImport)
        cls.assertIsNotNone(protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return protImport

    @classmethod
    def runImportMics(cls, filesPath, samplingRate):
        """ Import micrographs. """
        print(magentaStr("\n==> Importing data - micrographs:"))
        protImport = cls.newProtocol(ProtImportMicrographs,
                                     filesPath=filesPath,
                                     samplingRate=samplingRate)
        cls.launchProtocol(protImport)
        cls.assertIsNotNone(protImport.outputMicrographs,
                            "SetOfMicrographs has not been produced.")

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
    def runImportMovies(cls, filesPath=None, filesPattern=None, samplingRate=1.0,
                        dose=1.0, voltage=300, sphericalAberration=2.7, gain=None):
        """ Run an Import movies protocol. """
        print(magentaStr("\n==> Importing data - movies:"))
        protImport = cls.newProtocol(ProtImportMovies,
                                     filesPath=filesPath,
                                     filesPattern=filesPattern,
                                     voltage=voltage,
                                     sphericalAberration=sphericalAberration,
                                     samplingRate=samplingRate,
                                     dosePerFrame=dose,
                                     gainFile=gain)
        cls.launchProtocol(protImport)
        return protImport
