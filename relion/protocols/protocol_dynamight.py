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
from typing import List

import pyworkflow.protocol.params as params
from pyworkflow.constants import NEW
from pwem.protocols import ProtAnalysis3D
from pwem.constants import ALIGN_PROJ

from relion import Plugin
import relion.convert as convert
from relion.protocols.protocol_base import ProtRelionBase


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
    _possibleOutputs = {}
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
                      pointerClass='Volume', allowsNull=True,
                      label="Input consensus map",
                      condition='not doContinue')

        form.addParam('numberOfGaussians', params.IntParam, default=10000,
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
        form.addParam('threshold', params.FloatParam, default=-1,
                      condition='not doContinue',
                      label="Initial map threshold (optional)",
                      help="If provided, this threshold will be used to position "
                           "initial Gaussians in the consensus map. If left "
                           "default (-1), an automated procedure will be used to "
                           "estimate the appropriate threshold.")
        form.addParam('regularizeFactor', params.IntParam, default=1,
                      condition='not doContinue',
                      label="Regularization factor",
                      help="This regularization factor defines the relative "
                           "weights between the data over the restraints. "
                           "Values higher than one will put more weights on the "
                           "restraints.")
        form.addParam('allParticlesRam', params.BooleanParam, default=False,
                      label='Pre-read all particles into RAM?',
                      help="If set to Yes, dynamight will preload images into "
                           "memory for learning the forward or inverse deformations "
                           "and for deformed backprojection. This will speed up "
                           "the calculations, but you need to make sure you have "
                           "enough RAM to do so.")

        form.addSection(label='Tasks', condition='doContinue')

        group = form.addGroup('Visualize')
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

        group = form.addGroup('Deformations')
        group.addParam('doDeformEstimate', params.BooleanParam, default=False,
                       condition='doContinue',
                       label="Do inverse-deformation estimation?",
                       help="If set to Yes, dynamight will be run to estimate "
                            "inverse-deformations. These are necessary if one "
                            "want to perform deformed backprojection to calculate "
                            "an improved consensus model.")
        group.addParam('numEpochs', params.IntParam, default=200,
                       condition='doContinue and doDeformEstimate',
                       label="Number of epochs to perform",
                       help="Number of epochs to perform inverse deformations. "
                            "You can monitor the convergence of the loss "
                            "function to assess how many are necessary. "
                            "Often 200 are enough")

        group = form.addGroup('Backprojection')
        group.addParam('doDeformBackProj', params.BooleanParam, default=False,
                       condition='doContinue',
                       label="Do deformed backprojection?",
                       help="If set to Yes, dynamight will be run to perform "
                            "a deformed backprojection, using "
                            "inverse-deformations from a previous task, to "
                            "get an improved consensus reconstruction.")
        group.addParam('batchSize', params.IntParam, default=10,
                       condition='doContinue and doDeformBackProj',
                       label="Backprojection batchsize",
                       help="Number of images to process in parallel. "
                            "This will speed up the calculation, but will "
                            "cost GPU memory. Try how high you can go on "
                            "your GPU, given your box size and size of the "
                            "neural network.")

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()

        if not self.doContinue:
            self._insertFunctionStep(self.convertInputStep)
            self._insertFunctionStep(self.runDynamightStep)
        else:
            self._insertFunctionStep(self.runTasksStep)

        #self._insertFunctionStep(self.createOutputStep)

    # -------------------------- STEPS functions ------------------------------
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        deform_path = "forward_deformations/checkpoints"
        myDict = {
            'input_particles': self._getExtraPath('input_particles.star'),
            'checkpoint_iter': self._getExtraPath(deform_path,
                                                  '%(iter)03d.pth'),
            'checkpoint_final': self._getExtraPath(deform_path,
                                                   'checkpoint_final.pth'),
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
            f"--n-gaussians {self.numberOfGaussians}",
            f"--initial-threshold {self.threshold}",
            f"--regularization-factor {self.regularizeFactor}",
            f"--n-threads {self.numberOfThreads}",
            f"--gpu-id {self.gpuList.get()}",
            "--preload-images" if self.allParticlesRam else ""
        ]

        self.runProgram(params)

    def runTasksStep(self):
        inputProt = self.continueRun.get()
        inputProt._createFilenameTemplates()
        checkpoint_file = inputProt._getFileName('checkpoint_final')

        if self.doVisualize:
            params = [
                "explore-latent-space",
                self._getExtraPath(),
                f"--half-set {self.halfSet.get()}",
                f"--checkpoint-file {checkpoint_file}",
                f"--gpu-id {self.gpuList.get()}"
            ]
        elif self.doDeformEstimate:
            params = [
                "optimize-inverse-deformations",
                self._getExtraPath(),
                f"--n-epochs {self.numEpochs.get()}",
                f"--checkpoint-file {checkpoint_file}",
                f"--gpu-id {self.gpuList.get()}",
                "--preload-images" if self.allParticlesRam else ""
            ]

        elif self.doDeformBackProj:
            params = [
                "deformable-backprojection",
                self._getExtraPath(),
                f"--batch-size {self.batchSize.get()}",
                f"--checkpoint-file {checkpoint_file}",
                f"--gpu-id {self.gpuList.get()}",
                "--preload-images"
            ]
        else:
            raise ValueError("Unrecognized task")

        self.runProgram(params)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        return []

    def _validate(self):
        errors = []
        tasks = [self.doVisualize.get(),
                 self.doDeformEstimate.get(),
                 self.doDeformBackProj.get()]

        if tasks.count(True) > 1:
            errors.append("You cannot select multiple tasks")

        return errors

    def _citations(self):
        return ['Schwab2023']

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
