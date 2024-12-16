# ***************************************************************************
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
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
# ***************************************************************************/

import os
from tempfile import NamedTemporaryFile
import numpy as np

from pyworkflow.tests import setupTestProject
from pyworkflow.constants import SCIPION_DEBUG_NOCLEAN
import pyworkflow.utils as pwutils
from pwem.tests.workflows import TestWorkflow
from pwem.convert import AtomicStructHandler
from pwem.protocols import ProtImportParticles, ProtImportVolumes

from relion.protocols import ProtRelionReconstruct, ProtRelionDynaMight


class TestWorkflowRelionDynamight(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.sampling = 1.35
        cls.size = 128
        cls.gpuList = '0'

        from pwem import Domain
        try:
            cls.xmipp3 = Domain.importFromPlugin('xmipp3', doRaise=True)
            cls.xmippAvailable = True
        except:
            cls.xmippAvailable = False
        cls.absTmpPath = os.path.join(cls.getOutputPath(), cls.proj.getTmpPath())
        pwutils.cleanPath(cls.absTmpPath)

    def runCreateMask(self, pattern, thr):
        """ Create a volume mask. """
        from xmipp3.protocols import XmippProtCreateMask3D

        msk = self.newProtocol(XmippProtCreateMask3D,
                               inputVolume=pattern,
                               volumeOperation=0,  # OPERATION_THRESHOLD,
                               threshold=thr,
                               doSmall=False,
                               smallSize=False,
                               doBig=False,
                               doSymmetrize=False,
                               doMorphological=False,
                               doInvert=False,
                               doSmooth=True,
                               sigmaConvolution=2
                               )
        self.launchProtocol(msk)

        return msk

    def runApplyMask(self, volume, mask):
        """ Apply a mask to a volume. """
        from xmipp3.protocols import XmippProtMaskVolumes

        protApplyMask = self.newProtocol(XmippProtMaskVolumes,
                                         inputVolumes=volume,
                                         source=1,  # SOURCE_VOLUME
                                         inputMask=mask
                                         )
        self.launchProtocol(protApplyMask)

        return protApplyMask

    def runImportVolume(self, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        protImportVol = self.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        self.launchProtocol(protImportVol, wait=True)
        self.assertIsNotNone(protImportVol.outputVolume,
                             "Output volume not created")

        return protImportVol

    def runImportParticles(self, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImportPart = self.newProtocol(ProtImportParticles,
                                          objLabel='from xmipp3 (particles)',
                                          importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                          mdFile=pattern,
                                          samplingRate=samplingRate,
                                          voltage=300,
                                          sphericalAberration=2.0,
                                          amplitudeContrast=0.1
                                          )
        self.launchProtocol(protImportPart, wait=True)
        self.assertIsNotNone(protImportPart.outputParticles,
                             "Output particles not created")

        return protImportPart

    def runRelionReconstruct(self, particles):
        recProt1 = self.newProtocol(
            ProtRelionReconstruct,
            numberOfMpis=4,
            symmetryGroup='C1',
            doCTF=True,
            objLabel='Fourier reconstruction',
            inputParticles=particles)
        self.launchProtocol(recProt1)

        return recProt1

    @classmethod
    def getAtomStructFile(cls, atomStructID):
        aSH = AtomicStructHandler()
        atomStructPath = aSH.readFromPDBDatabase(atomStructID, dir=cls.absTmpPath, type='pdb')
        # filter out HETATM
        tempFile  = os.path.join(cls.absTmpPath, "kk")
        os.system(
            f"grep -v HETATM {atomStructPath} > {tempFile};"
            f"mv {tempFile} {atomStructPath}")

        return atomStructPath

    @classmethod
    def createVolume(cls, atomStructPath, volMapName):
        """ Create a volume from PDB. """
        cls.xmipp3.Plugin.runXmippProgram('xmipp_volume_from_pdb',
                                          f'-i {atomStructPath}'
                                          f' -o {volMapName} --size {cls.size}'
                                          f' --centerPDB --sampling {cls.sampling}')

        return volMapName + '.vol'

    @staticmethod
    def shiftFirstChain(atomStructPath,
                        translation_vector=np.array([0.0, 0.0, 0.0]),
                        nTimes=1):
        # import PDB needs new version of Biopython
        # leave the import here so it is not executed
        # during the test detection
        from Bio.PDB import PDBParser, PDBIO

        def translate_structure(structure, vector):
            """
            Translate the structure by a given vector.

            Args:
                structure: Biopython structure object.
                vector: List or tuple of (x, y, z)
                        translation values.
            """
            vector = vector * nTimes
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            # Translate the atom's coordinates
                            new_coord = atom.coord + vector
                            atom.coord = new_coord
                    break  # Break after the first chain

        parser = PDBParser(QUIET=True)
        print("atomStructPath", atomStructPath)
        structure = parser.get_structure("example_structure", atomStructPath)
        translate_structure(structure, translation_vector)
        io = PDBIO()
        io.set_structure(structure)
        output_path = atomStructPath.replace(".ent", f"_shifted_{nTimes}.pdb")
        io.save(output_path)

        return output_path

    @classmethod
    def createProjections(cls, volumeNames, angular_sampling_rate=15):
        projStackNames = []
        for vol in volumeNames:
            print("Creating projections for volume:",  vol)
            args = f" -i {vol}"
            args += f" -o {vol.replace('.vol', '.mrcs')}"
            args += f" --sampling_rate {angular_sampling_rate}"
            args += " --sym c1"
            args += " --method real_space"
            progname = "xmipp_angular_project_library"
            cls.xmipp3.Plugin.runXmippProgram(progname, args)
            projStackNames.append(vol.replace(".vol", ".doc"))

        return projStackNames

    @classmethod
    def createCTFdata(cls):
        ctfFile = NamedTemporaryFile(delete=False, suffix=".ctfdata")
        command = f"""# XMIPP_STAR_1 *
#  SamplingRate should be the same that the one used in the micrographs
data_fullMicrograph
 _ctfSamplingRate {cls.sampling}
 _ctfVoltage 300
 _ctfDefocusU 6000
 _ctfDefocusV 6000
 _ctfDefocusAngle -140.258
 _ctfSphericalAberration 2
 _ctfQ0 0.1
"""
        ctfFile.write(command.encode('utf8'))
        ctfFile.close()

        return ctfFile.name

    @classmethod
    def unionSets(cls, xmdProjNames):
        for i in range(1, len(xmdProjNames)):
            args = f" -i {xmdProjNames[0]}"
            args += f" --set union {xmdProjNames[i]}"
            args += f" -o {xmdProjNames[0]}"
            progname = "xmipp_metadata_utilities"
            cls.xmipp3.Plugin.runXmippProgram(progname, args)

        return xmdProjNames[0]

    @classmethod
    def addNoiseAndCTF(cls, projDocNames):
        ctfFile = cls.createCTFdata()
        print("ctfFile", ctfFile)
        xmdProjNames = []
        for projDocName in projDocNames:
            args = f" -i {projDocName}"
            args += f" -o {projDocName.replace('.doc', 'ctf.mrcs')}"
            args += f" --ctf {ctfFile}"
            args += " --noise 10 --after_ctf_noise"
            progname = "xmipp_phantom_simulate_microscope"
            cls.xmipp3.Plugin.runXmippProgram(progname, args)
            xmdProjNames.append(projDocName.replace(".doc", "ctf.xmd"))

        xmdProjName = cls.unionSets(xmdProjNames)

        return xmdProjName

    def importData(self, xmdProjName, volName):
        # import particles
        protImportProj = self.runImportParticles(xmdProjName, self.sampling)
        # import volume
        protImportVolume = self.runImportVolume(volName, self.sampling)

        return protImportProj, protImportVolume

    def runDynaMight(self, volume, particles):
        protDynaMight = self.newProtocol(
            ProtRelionDynaMight,
            referenceVolume=volume,
            inputParticles=particles,
            doShow=False,
            gpuList=self.gpuList,
            numberOfGaussians=10000,
            numberOfThreads=4,
            doDeform=False)

        return self.launchProtocol(protDynaMight, wait=True)

    def protDynaMightShow(self, protDynaMight):
        protDynaMight = self.newProtocol(
            ProtRelionDynaMight,
            continueRun=protDynaMight,
            doContinue=True,
            doVisualize=True,
            gpuList=self.gpuList,
            halfSet=0)
        # since this is an interactive protocol
        # do not run it, just save it
        # the user may execute it later
        self.saveProtocol(protDynaMight)

    def testDynamight(self):
        if not self.xmippAvailable:
            # if xmipp is not available, just skip this test
            return
        atomStructID = '3wtg'
        nVolumes = 2  # realistic value 5
        angular_sampling_rate = 30  # realistic value 3
        volMapName = os.path.join(self.absTmpPath, f'{atomStructID}_shifted_%d')
        atomStructPath = self.getAtomStructFile(atomStructID)
        translation_vector = np.array([0.5, 0.0, -0.25])
        volumeNames = []

        for i in range(nVolumes):
            shiftChain = self.shiftFirstChain(
                atomStructPath=atomStructPath,
                translation_vector=translation_vector,
                nTimes=i)
            volumeName = self.createVolume(
                atomStructPath=shiftChain,
                volMapName=volMapName % i)
            volumeNames.append(volumeName)

        projDocNames = self.createProjections(
            volumeNames, angular_sampling_rate)
        xmdProjFile = self.addNoiseAndCTF(projDocNames)
        volName = volumeNames[0]

        protImportProj, protImportVolume = self.importData(xmdProjName=xmdProjFile,
                                                           volName=volName)

        protCreateMask = self.runCreateMask(pattern=protImportVolume.outputVolume, thr=0.1)
        self.assertIsNotNone(protCreateMask.outputMask,
                             "Output mask not created")

        protRelionReconstruct = self.runRelionReconstruct(protImportProj.outputParticles)
        self.assertIsNotNone(protRelionReconstruct.outputVolume,
                             "Failed to reconstruct the volume")

        protApplyMask = self.runApplyMask(protRelionReconstruct.outputVolume,
                                          protCreateMask.outputMask)
        self.assertIsNotNone(protApplyMask.outputVol,
                             "Failed to apply a mask to the volume")

        protDynaMight = self.runDynaMight(protApplyMask.outputVol,
                                          protImportProj.outputParticles)
        protDynaMight._createFilenameTemplates()
        self.assertTrue(os.path.exists(protDynaMight._getFileName('checkpoint_final')),
                        "DynaMight has not produced expected output")

        self.protDynaMightShow(protDynaMight)
        if not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN):
            pwutils.cleanPath(self.absTmpPath)
