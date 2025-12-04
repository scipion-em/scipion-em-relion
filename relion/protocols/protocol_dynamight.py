# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
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
import os.path
from glob import glob
from typing import List

import pyworkflow.protocol.params as params
from pyworkflow.constants import NEW
import pyworkflow.utils as pwutils
from pwem.protocols import ProtAnalysis3D
from pwem.constants import ALIGN_PROJ
from pwem.objects import SetOfVolumes, Volume

from relion import Plugin
import relion.convert as convert
from relion.protocols.protocol_base import ProtRelionBase

TSNE = 0
UMAP = 1
PCA = 2
ICA = 3

class ProtRelionDynaMight(ProtAnalysis3D, ProtRelionBase):
    """
    Relion protocol for continuous flexibility analysis.

    As of release 5.0, Relion comes with a machine-learning approach for
    the analysis of molecular motions and flexibility called DynaMight.
    DynaMight will fit molecular motions for each experimental particle image
    as a 3D deformation field, which is learnt using a variational auto-encoder.
    It also implements functionality to calculate a pseudo-inverse 3D deformation
    field that can then be used in a deformed weighted backprojection algorithm
    to obtain an improved 3D reconstruction of the consensus structure.

    """
    _label = 'DynaMight flexibility'
    _devStatus = NEW
    _possibleOutputs = {"Volumes": SetOfVolumes}
    IS_CLASSIFY = False

    @classmethod
    def isDisabled(cls):
        return not Plugin.IS_GT50()

    def _initialize(self):
        """ This function is meant to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineConstants()
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU ID",
                       help="GPU may have several cores. Set it to zero "
                            "if you do not know what we are talking about. "
                            "First core index is 0, second 1 and so on. "
                            "*DynaMight can use only one GPU*.")

        form.addSection(label='Input')
        form.addParam('doContinue', params.BooleanParam, default=False,
                      label='Analyse a previous run?',
                      help='If you set to *Yes*, you should select a previous '
                           'DynaMight protocol and most of the input parameters '
                           'will be taken from it.')

        form.addParam('continueRun', params.PointerParam,
                      pointerClass='ProtRelionDynaMight',
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to analyse.')

        form.addParam('inputParticles', params.PointerParam,
                      allowsNull=True,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      condition='not doContinue',
                      label="Input particles",
                      help='Select the input images from the project.')

        form.addParam('referenceVolume', params.PointerParam,
                      pointerClass='VolumeMask', allowsNull=True,
                      label="Input consensus mask",
                      condition='not doContinue')

        form.addSection(label='Tasks')
        group = form.addGroup('Deformations', condition='not doContinue')
        group.addParam('numberOfGaussians', params.IntParam, default=5000,
                      condition='not doContinue',
                      label="Number of Gaussians",
                      help="Number of Gaussians to describe the consensus "
                           "map with. Larger structures that one wishes to "
                           "describe at higher resolutions will need more "
                           "Gaussians. As a rule of thumb, you could try and "
                           "use 1-2 Gaussians per amino acid or nucleotide in "
                           "your complex. But note that running DynaMight with "
                           "more than 30,000 Gaussians may be problematic on "
                           "GPUs with a memory of 24 GB.")

        group.addParam('threshold', params.FloatParam, default=-1,
                       condition='not doContinue',
                       label="Initial map threshold (optional)",
                       help="If provided, this threshold will be used to position "
                            "initial Gaussians in the consensus map. If left "
                            "default (-1), an automated procedure will be used to "
                            "estimate the appropriate threshold.")

        group.addParam('regularizeFactor', params.IntParam, default=1,
                       condition='not doContinue',
                       label="Regularization factor",
                       expertLevel=params.LEVEL_ADVANCED,
                       help="This regularization factor defines the relative "
                            "weights between the data over the restraints. "
                            "Values higher than one will put more weights on the "
                            "restraints.")

        group.addParam('latentDim', params.IntParam, default=8,
                       condition='not doContinue',
                       label="Latent dimension",
                       expertLevel=params.LEVEL_ADVANCED,
                       help="Number of latent dimension in encoded deformed latent space.")

        group.addParam('weightDecay', params.FloatParam, default=0.0,
                       condition='not doContinue',
                       label="Weight Decay",
                       expertLevel=params.LEVEL_ADVANCED,
                       help="Weight decay for Adam optimizer.")

        group.addParam('batchSizeD', params.IntParam, default=128,
                       condition='not doContinue',
                       label="Deformations batch size",
                       help="Batch size for processing images.")

        group.addParam('numEpochsD', params.IntParam, default=100,
                       condition='not doContinue',
                       label="Number of epochs",
                       help="Number of epochs for training network.")

        group.addParam('numWorkers', params.IntParam, default=8,
                       condition='not doContinue',
                       label="Number of workers in CPU",
                       help="Number of workers for multiple processes and loading in CPU.")

        group.addParam('allParticlesRam', params.BooleanParam, default=False,
                       label='Pre-read all particles into RAM?',
                       expertLevel=params.LEVEL_ADVANCED,
                       help="If set to Yes, dynamight will preload images into "
                            "memory for learning the forward or inverse deformations "
                            "and for deformed backprojection. This will speed up "
                            "the calculations, but you need to make sure you have "
                            "enough RAM to do so.")

        group = form.addGroup('Latent Space', condition='doContinue')
        group.addParam('doVisualize', params.BooleanParam, default=False,
                       condition='doContinue',
                       label="Do visualization?",
                       help="If set to Yes, dynamight will be run to visualize "
                            "the latent space and deformed models. One can also "
                            "save series of maps to make movies in Chimera, or "
                            "STAR files of particle subsets within this task.")

        group.addParam('halfSet', params.IntParam, default=0,
                       condition='doContinue and doVisualize',
                       label="Half-set to visualize",
                       help="Select halfset 1 or 2 to explore the latent space "
                            "of that halfset. If you select halfset 0, then the "
                            "validation set is being visualised, which will give "
                            "you an estimate of the errors in the deformations.")

        group.addParam('dimRed', params.EnumParam, default=PCA,
                       choices=['TSNE', 'UMAP', 'PCA', 'ICA'],
                       condition='doContinue and doVisualize',
                       label='Dimensionality reduction method',
                       help='Type of DimRed method to use when computing the corresponding latent space.')

        group.addParam('clusterize', params.BooleanParam, default=False,
                       condition='doContinue and doVisualize',
                       label='Clusterize latent space?',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Condition to clusterize latent space with KMEANS to distinguish in the visualization.')

        group = form.addGroup('Inverse Deformations', condition='doContinue')
        group.addParam('doDeform', params.BooleanParam, default=False,
                       condition='doContinue',
                       label="Estimate inverse deformation and backproject?",
                       help="If set to Yes, dynamight will be run to estimate "
                            "inverse-deformations first. These are necessary "
                            "to perform deformed backprojection to calculate "
                            "an improved consensus model.")

        group.addParam('numEpochsI', params.IntParam, default=200,
                       condition='doContinue and doDeform',
                       label="Number of epochs to perform",
                       help="Number of epochs to perform inverse deformations. "
                            "You can monitor the convergence of the loss "
                            "function to assess how many are necessary. "
                            "Often 200 are enough.")

        group.addParam('storeDeforms', params.BooleanParam, default=False,
                       condition='doContinue and doDeform',
                       label="Store deformations in RAM?",
                       expertLevel=params.LEVEL_ADVANCED,
                       help="If set to Yes, dynamight will store deformations "
                            "in the GPU memory, which will speed up the "
                            "calculations, but you need to have enough GPU "
                            "memory to do this.")

        group.addParam('batchSizeI', params.IntParam, default=10,
                       condition='doContinue and doDeform',
                       label="Backprojection batch size",
                       help="Number of images to process in parallel. "
                            "This will speed up the calculation, but will "
                            "cost GPU memory. Try how high you can go on "
                            "your GPU, given your box size and size of the "
                            "neural network.")

        group.addParam('downFactor', params.IntParam, default=2,
                       condition='doContinue and doDeform',
                       label='Downsampling factor for IT',
                       help='Downsampling factor to decrease IT computation to a smaller box. It is then upsampled'
                            ' to its original size.')

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()

        if not self.doContinue:
            self._insertFunctionStep(self.convertInputStep, needsGPU=False)
            self._insertFunctionStep(self.runDynamightStep, needsGPU=True)
        else:
            self.runTasks()
            self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        deform_path = "forward_deformations/checkpoints"
        myDict = {
            'input_particles': self._getExtraPath('input_particles.star'),
            'checkpoint_iter': self._getExtraPath(deform_path, '%(iter)03d.pth'),
            'checkpoint_final': self._getExtraPath(deform_path, 'checkpoint_final.pth')
            }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self, *args):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self.inputParticles.get()
        imgStar = self._getFileName('input_particles')

        # Pass stack file as None to avoid write the images files
        convert.writeSetOfParticles(imgSet, imgStar,
                                    outputDir=self._getExtraPath(),
                                    alignType=ALIGN_PROJ)
        self._convertRef()

    def runDynamightStep(self):
        params = [
            "optimize-deformations",
            f"--refinement-star-file {self._getFileName('input_particles')}",
            f"--output-directory {self._getExtraPath()}",
            f"--initial-model {self._getRefArg()}",
            f"--initial-threshold {self.threshold.get()}",
            f"--mask-file {self.referenceVolume.get()}" if self.referenceVolume else ""
            f"--n-gaussians {self.numberOfGaussians.get()}",
            f"--n-latent-dimensions {self.latentDim.get()}",
            f"--weight-decay {self.weightDecay.get()}",
            f"--regularization-factor {self.regularizeFactor.get()}",
            f"--batch-size {self.batchSizeD.get()}",
            f"--gpu-id {self.gpuList.get()}",
            f"--n-epochs {self.numEpochsD.get()}",
            f"--n-threads {self.numberOfThreads.get()}",
            "--preload-images" if self.allParticlesRam else "",
            f"--n-workers {self.numWorkers.get()}"
        ]

        self.runProgram(params)

    def runTasks(self):
        inputProt = self.continueRun.get()
        pwutils.createLink(inputProt._getExtraPath("forward_deformations"),
                           self._getExtraPath("forward_deformations"))
        checkpoint_file = self._getFileName('checkpoint_final')

        if self.doVisualize:
            params = [
                "explore-latent-space",
                self._getExtraPath(),
                f"--checkpoint-file {checkpoint_file}",
                f"--half-set {self.halfSet.get()}",
                f"--mask-file {inputProt.referenceVolume.get()}" if inputProt.referenceVolume else "",
                f"--batch-size {inputProt.batchSizeD.get()}",
                f"--gpu-id {self.gpuList.get()}",
                f"--n-workers {inputProt.numWorkers.get()}",
                f"--dimensionality-reduction-method {self.dimRed.get()}",
                f"--cluster" if self.clusterize else ""
            ]
            self._insertFunctionStep(self.runTaskStep, params, needsGPU=True)

        elif self.doDeform:
            params = [
                "optimize-inverse-deformations",
                self._getExtraPath(),
                f"--checkpoint-file {checkpoint_file}",
                f"--batch-size {self.batchSizeI.get()}",
                f"--n-epochs {self.numEpochsI.get()}",
                f"--gpu-id {self.gpuList.get()}",
                "--preload-images" if self.allParticlesRam else "",
                f"--data-loader-threads {self.numberOfThreads.get()}",
                "--save-deformations" if self.storeDeforms else ""
            ]
            self._insertFunctionStep(self.runTaskStep, params, needsGPU=True)

            params = [
                "deformable-backprojection_correction",
                self._getExtraPath(),
                f"--mask-file {inputProt.referenceVolume.get()}" if inputProt.referenceVolume else "",
                f"--gpu-id {self.gpuList.get()}",
                f"--batch-size {self.batchSizeI.get()}",
                "--preload-images" if self.allParticlesRam else "",
                f"--data-loader-threads {self.numberOfThreads.get()}",
                f"--downsample {self.downFactor.get()}"
            ]
            self._insertFunctionStep(self.runTaskStep, params, needsGPU=True)

    def runTaskStep(self, params: List[str]):
        """ Run the actual job. """
        self.runProgram(params)

    def createOutputStep(self):
        parts = self._getInputParticles()
        if self.doVisualize:
            output = self._getExtraPath("maps", "map???_half?.mrc")
        elif self.doDeform:
            output = self._getExtraPath("backprojection", "map_half?.mrc")

        files = sorted(glob(output))
        if files:
            volumes = self._createSetOfVolumes()
            volumes.setSamplingRate(parts.getSamplingRate())

            for i, volFn in enumerate(files):
                vol = Volume()
                vol.setFileName(volFn)
                if self.doDeform:
                    vol.setObjLabel(f"half-map {i+1}")
                volumes.append(vol)

            self._defineOutputs(Volumes=volumes)
            self._defineSourceRelation(parts, volumes)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.isFinished() and hasattr(self, "Volumes"):
            if self.doVisualize:
                summary.append("Movies/maps generated along the deformation trajectory")
            else:
                summary.append("Half-maps generated by inverse deformation")

        return summary

    def _validate(self):
        errors = []
        tasks = [self.doVisualize.get(), self.doDeform.get()]

        if tasks.count(True) > 1:
            errors.append("You cannot select multiple tasks")

        if self.doContinue:
            inputProt = self.continueRun.get()
            inputProt._createFilenameTemplates()
            checkpoint_file = inputProt._getFileName('checkpoint_final')
            if not os.path.exists(checkpoint_file):
                errors.append(f"Cannot continue, {checkpoint_file} not found")

        return errors

    def _citations(self):
        return ['Schwab2024']

    # -------------------------- UTILS functions ------------------------------
    def runProgram(self, params: List[str]) -> None:
        program = "relion_python_dynamight"
        self.runJob(f"{Plugin.getActivationCmd()} && {program}",
                    " ".join(params))

    def _getEnviron(self):
        env = Plugin.getEnviron()
        if 'LD_LIBRARY_PATH' in env:
            # this is required to avoid conflict btw DynaMight Qt5 libs
            # and system Qt libs
            del env['LD_LIBRARY_PATH']

        return env
