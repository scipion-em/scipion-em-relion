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

from enum import Enum

from pyworkflow.protocol.params import (PointerParam, FloatParam,  
                                        StringParam, BooleanParam,
                                        EnumParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.constants import PROD
from pwem.objects import Volume
from pwem.protocols import ProtReconstruct3D
from pwem.constants import ALIGN_PROJ

import relion.convert as convert
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputVolume = Volume


class ProtRelionReconstruct(ProtReconstruct3D, ProtRelionBase):
    """
    Reconstructs a 3D cryo-EM volume from aligned particle images using Relion.
    The protocol combines experimental particle projections into a volumetric
    reconstruction while preserving the orientation information associated with
    each particle image.

    AI Generated:

    Reconstruct Volume (ProtRelionReconstruct) — User Manual
        Overview

        The Reconstruct Volume protocol generates a 3D density map from a set
        of aligned cryo-EM particle images using the Relion reconstruction
        framework. Its primary purpose is to transform two-dimensional particle
        projections into a biologically interpretable three-dimensional map that
        represents the underlying molecular structure. This protocol is commonly
        used after particle alignment, classification, or refinement procedures
        in which particle orientations have already been estimated.

        In practical cryo-EM workflows, reconstruction is one of the central
        stages connecting particle-based image analysis with structural
        interpretation. Biological users typically employ this protocol to
        generate consensus maps, reconstruct specific particle subsets, or
        evaluate structural heterogeneity after classification. The resulting
        maps may subsequently be used for visualization, model fitting,
        refinement, variability analysis, or downstream structural validation.

        Inputs and General Workflow

        The protocol requires a set of particles with projection alignment
        information already assigned. These angular assignments define how each
        particle contributes to the final three-dimensional reconstruction.
        Accurate alignments are therefore critical for obtaining biologically
        meaningful results. Poor angular assignments, excessive heterogeneity,
        or strong particle contamination may lead to blurred or distorted maps.

        The workflow reconstructs the volume directly from the experimental
        particle images while preserving the imaging geometry estimated during
        previous refinement or classification stages. Users may reconstruct the
        entire dataset or focus on specific subsets, such as individual classes
        or half-sets. This flexibility is especially important when exploring
        conformational variability or validating reconstruction consistency.

        Symmetry Considerations

        Symmetry selection is one of the most biologically significant decisions
        during reconstruction. Applying symmetry can substantially improve the
        signal-to-noise ratio by averaging equivalent structural views, often
        leading to higher apparent resolution and cleaner density maps.
        Appropriate symmetry application is therefore highly beneficial for
        highly symmetric assemblies such as viral capsids, molecular cages, or
        oligomeric protein complexes.

        However, imposing incorrect symmetry may artificially distort the
        biological structure and hide genuine asymmetry or flexibility.
        Biological users should therefore apply only symmetry that is strongly
        supported by experimental evidence or prior structural knowledge. When
        uncertainty exists, reconstructing initially with no symmetry is often
        the safest strategy.

        Resolution and Fourier Space Control

        The protocol allows the reconstruction to be limited to a specified
        maximum resolution. This parameter determines the highest spatial
        frequencies included during reconstruction and can influence both map
        sharpness and stability. Restricting reconstruction to lower resolution
        ranges may improve robustness during exploratory analyses or when
        particle quality is limited.

        Padding options control Fourier-space interpolation behavior during
        reconstruction. Larger padding values may improve numerical accuracy at
        the expense of additional computational cost. In most biological
        workflows, default values are generally appropriate unless specialized
        optimization is required.

        Subset and Class Reconstruction

        The protocol supports reconstruction from selected subsets of particles.
        Users may reconstruct all particles together or separately reconstruct
        the independent half-sets commonly used in gold-standard validation
        workflows. Half-set reconstructions are particularly important for
        assessing resolution reproducibility and avoiding overfitting during
        refinement.

        In addition, the protocol can reconstruct only particles belonging to a
        selected class. This functionality is biologically valuable when
        studying compositional or conformational heterogeneity. Separate class
        reconstructions often reveal distinct structural states that may be
        obscured in a consensus reconstruction.

        CTF Correction

        Contrast Transfer Function correction is an essential component of
        accurate cryo-EM reconstruction. The protocol supports reconstruction
        with CTF correction enabled, allowing recovery of structural signal
        affected by microscope imaging physics. Proper CTF handling generally
        improves map interpretability and resolution.

        In some cases, users may choose to preserve the low-frequency region up
        to the first CTF peak. This approach can help when low-resolution CTF
        modeling is uncertain or unstable. If the input particles were already
        phase-flipped during preprocessing, the protocol can account for this
        condition to maintain consistency throughout reconstruction.

        From a biological perspective, reliable CTF estimation is critical for
        obtaining accurate structural features, especially at intermediate and
        high resolution where secondary structure elements and side-chain
        densities become visible.

        Ewald Sphere Correction

        For very high-resolution cryo-EM datasets or particularly large
        molecular assemblies, Ewald sphere curvature effects may become
        significant. The protocol provides optional Ewald sphere correction to
        improve reconstruction accuracy under these demanding conditions.

        Biological users typically consider this option only when pursuing
        near-atomic or atomic resolution structures. The correction may improve
        fine structural detail, particularly in large particles where curvature
        effects are stronger. However, the increased computational complexity
        means that this option is generally unnecessary for moderate-resolution
        studies.

        Additional controls allow users to define masking behavior, weighting
        strategies, reconstruction box size, and curvature orientation. These
        advanced parameters are mainly relevant for expert users performing
        specialized high-resolution optimization.

        Outputs and Their Interpretation

        The protocol produces a reconstructed three-dimensional volume together
        with its associated sampling information. The resulting map represents
        the combined structural signal extracted from the aligned particle set.
        Depending on the selected reconstruction strategy, the output may
        correspond to a consensus structure, a class-specific state, or a
        validation half-map.

        Biologically, the interpretation of the reconstruction depends strongly
        on particle quality, structural homogeneity, alignment accuracy, and
        symmetry assumptions. Well-resolved reconstructions may reveal domain
        organization, secondary structure elements, ligand binding regions, or
        conformational differences between states.

        Practical Recommendations

        For most standard cryo-EM workflows, reconstructing all aligned
        particles with appropriate symmetry and CTF correction provides a
        reliable starting point. Users studying structural variability should
        reconstruct individual classes separately to better isolate distinct
        conformations or compositional states.

        When pursuing high-resolution reconstructions, careful attention should
        be paid to symmetry selection, half-set consistency, and CTF quality.
        Ewald sphere correction should generally be reserved for advanced
        refinement stages where subtle high-resolution improvements become
        biologically meaningful.

        It is often advisable to inspect reconstructions visually after
        completion to confirm that the resulting density agrees with known
        structural expectations and does not display artifacts caused by
        incorrect symmetry, alignment instability, or excessive heterogeneity.

        Final Perspective

        Three-dimensional reconstruction is one of the defining steps in
        cryo-EM structural analysis because it transforms collections of noisy
        particle images into interpretable molecular structures. Reliable
        reconstruction depends not only on computational accuracy but also on
        biologically informed decisions regarding symmetry, particle selection,
        and structural heterogeneity. Careful reconstruction strategies provide
        the foundation for meaningful structural interpretation and subsequent
        biological discovery.
    """
    _label = 'reconstruct'
    _devStatus = PROD
    _possibleOutputs = outputs

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group",
                      help='See [[https://relion.readthedocs.io/'
                           'en/latest/Reference/Conventions.html#symmetry]'
                           '[Relion Symmetry]] page for a description '
                           'of the symmetry format accepted by Relion')
        form.addParam('maxRes', FloatParam, default=-1,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).')
        form.addParam('pad', FloatParam, default=2,
                      label="Padding factor")
        form.addParam('subset', EnumParam, default=0,
                      choices=['all', 'half1', 'half2'],
                      display=EnumParam.DISPLAY_HLIST,
                      label='Subset to reconstruct',
                      help='Subset of images to consider.')
        form.addParam('classNum', IntParam, default=-1,
                      label='Use only this class',
                      help='Consider only this class (-1: use all classes)')
        
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters: ', 
                      help='Extra parameters to *relion_reconstruct* program. '
                           'Address to Relion to see full list of options.')
        form.addSection('CTF')
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Apply CTF correction?')
        form.addParam('ctfIntactFirstPeak', BooleanParam, default=False,
                      condition='doCTF',
                      label='Leave CTFs intact until first peak?')

        form.addSection('Ewald sphere')
        form.addParam('doEwald', BooleanParam, default=False,
                      label="Correct for Ewald-sphere curvature?")
        form.addParam('skipMask', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Skip masking?",
                      help="Do not apply real space mask during Ewald "
                           "sphere correction.")
        form.addParam('maskDiameterA', IntParam, default=-1,
                      condition='not skipMask',
                      label='Mask diameter (A)',
                      help='Diameter (in A) of mask for Ewald-sphere '
                           'curvature correction')
        form.addParam('edge', IntParam, default=3,
                      condition='not skipMask',
                      label='Add a soft-edge (px)',
                      help='Width (in pixels) of the soft edge on the mask.')
        form.addParam('reverseCurvature', BooleanParam, default=False,
                      label="Reverse curvature?")
        form.addParam('newBoxSize', IntParam, default=-1,
                      expertLevel=LEVEL_ADVANCED,
                      label="New box size (px)",
                      help="Box size of reconstruction after Ewald "
                           "sphere correction.")
        form.addParam('numSectors', IntParam, default=2,
                      expertLevel=LEVEL_ADVANCED,
                      label="Number of sectors",
                      help="Number of sectors for Ewald sphere correction.")
        form.addParam('skipWeight', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Skip weighting?",
                      help="Do not apply weighting during during Ewald "
                           "sphere correction.")
        
        form.addParallelSection(threads=0, mpi=1)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertReconstructStep()
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    def _insertReconstructStep(self):
        imgSet = self.inputParticles.get()

        params = ' --i %s' % self._getFileName('input_particles')
        params += ' --o %s' % self._getFileName('output_volume')
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --angpix %0.5f' % imgSet.getSamplingRate()
        params += ' --maxres %0.3f' % self.maxRes.get()
        params += ' --pad %0.3f' % self.pad.get()

        subset = -1 if self.subset.get() == 0 else self.subset
        params += ' --subset %d' % subset
        params += ' --class %d' % self.classNum.get()

        if self.doCTF:
            params += ' --ctf'
            if self.ctfIntactFirstPeak:
                params += ' --ctf_intact_first_peak'

            if imgSet.isPhaseFlipped():
                params += ' --ctf_phase_flipped'

        if self.extraParams.hasValue():
            params += " " + self.extraParams.get()

        if self.doEwald:
            params += " --ewald --sectors %d --newbox %d" % (self.numSectors,
                                                             self.newBoxSize)
            if self.skipMask:
                params += " --skip_mask"
            else:
                params += " --mask_diameter %d --width_mask_edge %d" % (
                    self.maskDiameterA, self.edge
                )
            if self.skipWeight:
                params += " --skip_weighting"
            if self.reverseCurvature:
                params += " --reverse_curvature"

        self._insertFunctionStep(self.reconstructStep, params, needsGPU=False)

    # -------------------------- STEPS functions ------------------------------
    def reconstructStep(self, params):
        self._runProgram('relion_reconstruct', params)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'output_volume': self._getExtraPath('output_volume.mrc')
            }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self.inputParticles.get()
        imgStar = self._getFileName('input_particles')

        # Pass stack file as None to avoid write the images files
        convert.writeSetOfParticles(imgSet, imgStar,
                                    outputDir=self._getTmpPath(),
                                    alignType=ALIGN_PROJ)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        volume = Volume()
        volume.setFileName(self._getFileName('output_volume'))
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(**{outputs.outputVolume.name: volume})
        self._defineSourceRelation(self.inputParticles, volume)
    
    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume not ready yet.")
        else:
            summary.append("Output volume has been reconstructed.")

        return summary
