# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)           [1]
# *              Eduardo García Delgado (eduardo.garcia@cnb.csic.es)  [2]
# *              David Herreros (dherreros@cnb.csic.es)               [2]
# *              Mikel Iceta (miceta@cnb.csic.es)                     [2]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
# * [2] Unidad de  Biocomputacion, Centro Nacional de Biotecnologia, CSIC (CNB-CSIC)
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
import shutil
from glob import glob
from typing import List
import numpy as np

import pyworkflow.protocol.params as params
from joblib.testing import param
from pyworkflow.constants import NEW
import pyworkflow.utils as pwutils
from pwem.protocols import ProtAnalysis3D, ProtFlexBase
from pwem.constants import ALIGN_PROJ
from pwem.objects import SetOfVolumes, Volume, ParticleFlex

import relion
from relion import Plugin
import relion.convert as convert
from relion.protocols.protocol_base import ProtRelionBase
from relion.constants import DYNAMIGHT

TSNE = 0
UMAP = 1
PCA = 2
ICA = 3

class ProtRelionDynaMight(ProtAnalysis3D, ProtRelionBase, ProtFlexBase):
    """
    Relion protocol for continuous flexibility analysis.

    AI Generated:

    DynaMight Flexibility Analysis (ProtRelionDynaMight) - User Manual
        Overview

        The DynaMight protocol performs continuous flexibility analysis
        of cryo-EM particle datasets using machine learning methods
        integrated within Relion. Its primary goal is to describe
        structural heterogeneity as continuous molecular motions rather
        than as a limited number of discrete conformational states.
        This approach is especially valuable for flexible proteins,
        molecular machines, membrane complexes, and assemblies that
        undergo coordinated structural rearrangements.

        The protocol models conformational variability by learning
        deformation fields directly from experimental particle images.
        These deformations are represented in a low-dimensional latent
        space that captures the principal modes of motion present in
        the dataset. From a biological perspective, this enables the
        exploration of continuous transitions between structural states,
        helping researchers visualize flexibility pathways and identify
        functionally relevant motions.

        In addition to deformation analysis, the protocol also supports
        inverse deformation estimation and deformable backprojection.
        These operations can improve consensus reconstructions by
        compensating for conformational variability during map
        reconstruction, potentially increasing structural quality in
        flexible regions.

        Inputs and Biological Context

        The protocol requires a set of aligned particles together with
        a consensus reference volume. The reference volume defines the
        structural framework from which deformations are learned. In
        practice, this reference is usually a consensus reconstruction
        obtained after standard cryo-EM refinement workflows.

        Optionally, a consensus mask may be provided. From a biological
        perspective, masking is highly important because it determines
        which regions contribute most strongly to the flexibility
        analysis. Applying an appropriate mask is particularly useful
        for excluding solvent regions, detergent micelles, disordered
        density, or highly noisy peripheral domains that could obscure
        meaningful motions.

        The quality of the consensus volume strongly influences the
        interpretability of the results. Poorly resolved references or
        highly heterogeneous reconstructions may lead to unstable
        deformation models or biologically ambiguous latent spaces.
        Therefore, it is generally advisable to begin with the best
        available refinement and a carefully inspected consensus map.

        Learning Molecular Motions

        DynaMight represents the structure using a collection of
        Gaussian components distributed throughout the molecular volume.
        The number of Gaussians determines the complexity and spatial
        detail of the deformation model. Smaller and simpler complexes
        generally require fewer components, whereas large ribosomes,
        spliceosomes, viral assemblies, or highly detailed maps may
        require substantially more.

        From a biological perspective, increasing the number of
        Gaussians improves the ability to capture localized motions and
        subtle conformational changes. However, larger models also
        increase computational cost and GPU memory requirements.
        Selecting an excessively large number of components may lead to
        unstable training or impractical execution times on limited
        hardware.

        The protocol also includes regularization controls that balance
        fidelity to experimental data against smoothness and physical
        plausibility of the deformations. Stronger regularization tends
        to suppress noisy or unrealistic motions, while weaker
        regularization allows more flexibility but may overfit noise.
        In practice, moderate regularization is usually preferred for
        biological interpretation.

        Latent Space Representation

        One of the central outputs of the protocol is the latent space,
        which provides a compact numerical representation of molecular
        conformations. Each particle is associated with coordinates in
        this latent space, allowing the visualization and exploration
        of continuous structural variability.

        The latent dimensionality controls how many independent modes
        of motion can be represented. Lower-dimensional spaces are
        easier to interpret and often sufficient for systems dominated
        by a few large motions. More complex systems may require higher
        dimensionality to capture coupled or independent rearrangements.

        For biological users, latent space exploration can reveal
        trajectories between conformations, identify clusters of
        structural states, and provide insight into functional dynamics.
        This is particularly useful when studying ligand binding,
        allosteric regulation, domain breathing, gating transitions, or
        assembly maturation processes.

        The protocol supports dimensionality reduction methods such as
        PCA, ICA, UMAP, and t-SNE for visualization. These approaches
        provide different perspectives on the organization of the
        conformational landscape. PCA is often a good starting point
        for general interpretation, while nonlinear methods such as
        UMAP or t-SNE may better separate complex conformational
        relationships.

        Visualization and Deformation Analysis

        The visualization stage allows users to inspect maps generated
        along deformation trajectories. These trajectories provide an
        intuitive view of how the structure changes across the latent
        space and may reveal biologically meaningful motions that are
        difficult to detect using discrete classification approaches.

        From a biological perspective, interpreting these trajectories
        requires caution. Continuous transitions inferred from latent
        representations are mathematical approximations of the observed
        variability and should be evaluated alongside biochemical,
        structural, or functional evidence. Flexible regions with weak
        signal may produce motions that are visually plausible but not
        necessarily biologically relevant.

        Visualization is particularly powerful for identifying domain
        movements, hinge motions, flexible loops, opening and closing
        transitions, or coordinated rearrangements within large
        assemblies. Generating movies from deformation trajectories can
        significantly improve qualitative interpretation and
        communication of structural dynamics.

        Inverse Deformations and Backprojection

        The protocol can also estimate inverse deformation fields and
        perform deformable backprojection. This advanced workflow aims
        to reconstruct improved consensus maps by accounting for
        conformational variability during reconstruction.

        Biologically, this is especially valuable when flexibility
        limits local resolution in standard refinements. By compensating
        for particle-specific deformations, the protocol may recover
        structural detail that would otherwise be blurred in highly
        dynamic regions.

        The inverse deformation process is computationally demanding
        and typically benefits from modern GPUs with sufficient memory.
        Batch sizes, image preloading, and memory optimization settings
        can strongly influence performance. Users working with very
        large datasets should carefully balance computational efficiency
        against available hardware resources.

        Practical Recommendations

        For most cryo-EM studies, it is advisable to begin with a
        moderate latent dimensionality and a conservative number of
        Gaussian components. Initial exploratory analyses can then be
        refined based on the observed flexibility patterns and training
        stability.

        Applying an appropriate mask often produces the largest
        improvement in biological interpretability, particularly for
        complexes with flexible appendages or poorly resolved solvent
        regions. Careful inspection of deformation trajectories is
        essential to distinguish meaningful conformational changes from
        noise-driven artifacts.

        GPU memory considerations are especially important for large
        structures and high-resolution analyses. Increasing model
        complexity, batch size, or deformation storage may accelerate
        calculations but can rapidly exceed available hardware
        capacity.

        When interpreting latent spaces, users should avoid assigning
        direct biological meaning to every observed axis or trajectory.
        Instead, latent coordinates should be considered as simplified
        representations of conformational variability that require
        validation through complementary structural or biochemical
        evidence.

        Final Perspective

        Continuous flexibility analysis represents an important advance
        in cryo-EM structural biology because many biological systems
        do not exist as a small number of rigid states. DynaMight
        provides a framework for studying molecular dynamics directly
        from experimental particle images, enabling a richer
        interpretation of conformational landscapes and functional
        mechanisms.

        For biological users, the protocol is most effective when
        combined with careful consensus refinement, thoughtful masking,
        realistic computational settings, and cautious interpretation
        of deformation trajectories. When used appropriately, it can
        reveal dynamic structural relationships that remain hidden in
        conventional discrete classification workflows.
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
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      condition='not doContinue',
                      label="Input particles",
                      help='Input particles to run with.')

        form.addParam('referenceVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Input consensus volume",
                      condition='not doContinue')

        form.addParam('referenceMask', params.PointerParam,
                      pointerClass='VolumeMask', allowsNull=True,
                      label="Input consensus mask",
                      expertLevel=params.LEVEL_ADVANCED,
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

        group.addParam('numEpochsD', params.IntParam, default=200,
                       condition='not doContinue',
                       label="Number of epochs",
                       help="Number of epochs for training network.")

        group.addParam('numWorkers', params.IntParam, default=8,
                       condition='not doContinue',
                       label="Number of workers in CPU",
                       help="Number of workers for multiple processes and loading in CPU.")

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

        group = form.addGroup('Inverse Deformations', condition='doContinue')
        doInverse = 'doContinue and doDeform'
        group.addParam('doDeform', params.BooleanParam, default=False,
                       condition='doContinue',
                       label="Estimate inverse deformation and backproject?",
                       help="If set to Yes, dynamight will be run to estimate "
                            "inverse-deformations first. These are necessary "
                            "to perform deformed backprojection to calculate "
                            "an improved consensus model.")

        group.addParam('numEpochsI', params.IntParam, default=200,
                       condition=doInverse,
                       label="Number of epochs to perform",
                       help="Number of epochs to perform inverse deformations. "
                            "You can monitor the convergence of the loss "
                            "function to assess how many are necessary. "
                            "Often 200 are enough.")

        group.addParam('storeDeforms', params.BooleanParam, default=False,
                       condition=doInverse,
                       label="Store deformations in RAM?",
                       help="If set to Yes, dynamight will store deformations "
                            "in the GPU memory, which will speed up the "
                            "calculations, but you need to have enough GPU "
                            "memory to do this.")

        group.addParam('batchSizeI', params.IntParam, default=10,
                       condition=doInverse,
                       label="Backprojection batch size",
                       help="Number of images to process in parallel. "
                            "This will speed up the calculation, but will "
                            "cost GPU memory. Try how high you can go on "
                            "your GPU, given your box size and size of the "
                            "neural network.")

        group.addParam('downFactor', params.IntParam, default=2,
                       condition=doInverse,
                       label='Downsampling factor for IT',
                       help='Downsampling factor to decrease IT computation to a smaller box. It is then upsampled'
                            ' to its original size.')

        form.addParam('allParticlesRam', params.BooleanParam, default=False,
                       label='Pre-read all particles into RAM?',
                       help="If set to Yes, dynamight will preload images into "
                            "memory for learning the forward or inverse deformations "
                            "and for deformed backprojection. This will speed up "
                            "the calculations, but you need to make sure you have "
                            "enough RAM to do so.")

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()

        if not self.doContinue:
            self._insertFunctionStep(self.convertInputStep, needsGPU=False)
            self._insertFunctionStep(self.runDynamightStep, needsGPU=True)
            self._insertFunctionStep(self.createOutputTrainingStep, needsGPU=False)
        else:
            self.runTasks()
            self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        deform_path = "forward_deformations/checkpoints"
        myDict = {
            'input_particles': self._getExtraPath('input_particles.star'),
            'input_mask': self._getExtraPath('input_mask.mrc'),
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

        if self.referenceMask.get() is not None:
            print("Mask was provided, loading mask...")
            maskFilename = self._getFileName('input_mask')
            inMask = self.referenceMask.get().getFileName()
            shutil.copy(inMask, maskFilename)
        else:
            print("No mask provided, continuing without mask...")

    def runDynamightStep(self):
        hasMask = self.referenceMask.get() is not None
        params = [
            "optimize-deformations",
            f"--refinement-star-file {self._getFileName('input_particles')}",
            f"--output-directory {self._getExtraPath()}",
            f"--initial-model {self._getRefArg()}",
            f"--initial-threshold {self.threshold.get()}",
            f"--mask-file {self._getFileName('input_mask')}" if hasMask else "",
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

        # Predict latent space
        script_dynamight_encode = os.path.join(os.path.dirname(relion.__file__), "dynamight", "dynamight_encode_latent_vectors.py")
        params = [
            f"--output_directory {self._getExtraPath()}",
            f"--checkpoint_file {self._getFileName('checkpoint_final')}",
            f"--gpu_id {self.gpuList.get()}",
        ]
        self.runPythonScript(script_dynamight_encode, params)

    def runTasks(self):
        inputProt = self.continueRun.get()
        pwutils.createLink(inputProt._getExtraPath("forward_deformations"),
                           self._getExtraPath("forward_deformations"))
        checkpoint_file = self._getFileName('checkpoint_final')
        hasMask = inputProt.referenceMask.get() is not None

        if hasMask:
            maskFilename = self._getFileName('input_mask')
            inMask = inputProt.referenceMask.get().getFileName()
            shutil.copy(inMask, maskFilename)

        if self.doVisualize:
            dimRed = ['TSNE', 'UMAP', 'PCA', 'ICA']
            dimRed = dimRed[self.dimRed.get()]

            params = [
                "explore-latent-space",
                self._getExtraPath(),
                f"--checkpoint-file {checkpoint_file}",
                f"--half-set {self.halfSet.get()}",
                f"--mask-file {self._getFileName('input_mask')}" if hasMask else "",
                f"--batch-size {inputProt.batchSizeD.get()}",
                f"--gpu-id {self.gpuList.get()}",
                f"--n-workers {inputProt.numWorkers.get()}",
                f"--dimensionality-reduction-method {dimRed}"
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
                "deformable-backprojection",
                self._getExtraPath(),
                f"--mask-file {self._getFileName('input_mask')}" if hasMask else "",
                f"--gpu-id {self.gpuList.get()}",
                f"--backprojection-batch-size {self.batchSizeI.get()}",
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

    def createOutputTrainingStep(self):
        parts = self._getInputParticles()

        partSet = self._createSetOfParticlesFlex(progName=DYNAMIGHT)

        partSet.copyInfo(parts)
        partSet.setHasCTF(parts.hasCTF())
        partSet.setAlignmentProj()
        partSet.getFlexInfo().setAttr("checkpoint_file", self._getFileName('checkpoint_final'))

        # Load encoded latent vectors
        latent_vectors = np.load(self._getExtraPath("latent_vectors.npy"))

        idx = 0
        for particle in parts.iterItems():
            outParticle = ParticleFlex(progName=DYNAMIGHT)
            outParticle.copyInfo(particle)

            outParticle.setZFlex(latent_vectors[idx])
            outParticle.getFlexInfo().setAttr("checkpoint_file", self._getFileName('checkpoint_final'))

            partSet.append(outParticle)

            idx += 1

        self._defineOutputs(Particles=partSet)
        self._defineSourceRelation(parts, partSet)

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

    def runPythonScript(self, script: str, params: List[str]) -> None:
        program = "python"
        self.runJob(f"{Plugin.getActivationCmd()} && {program} {script}",
                    " ".join(params))

    def _getEnviron(self):
        env = Plugin.getEnviron()
        if 'LD_LIBRARY_PATH' in env:
            # this is required to avoid conflict btw DynaMight Qt5 libs
            # and system Qt libs
            del env['LD_LIBRARY_PATH']

        return env
