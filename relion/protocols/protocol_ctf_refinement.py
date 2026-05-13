# ******************************************************************************
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
# ******************************************************************************

from enum import Enum

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.constants import ALIGN_PROJ
from pwem.objects import SetOfParticles
from pwem.protocols import ProtParticles

import relion
import relion.convert as convert
from relion.convert.convert31 import Reader, OpticsGroups
from .protocol_base import ProtRelionBase

from ..objects import CtfRefineGlobalInfo


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtRelionCtfRefinement(ProtParticles, ProtRelionBase):
    """
    Refines contrast transfer function parameters and optical aberrations
    for cryo-EM particle datasets using Relion CTF refinement procedures.
    The protocol improves the accuracy of particle imaging parameters in
    order to enhance high-resolution reconstruction quality and reduce
    systematic optical errors during single-particle analysis.

    AI Generated:

    CTF Refinement (ProtRelionCtfRefinement) — User Manual
        Overview

        The CTF Refinement protocol performs advanced refinement of imaging
        parameters for cryo-EM particle datasets. Its primary goal is to
        improve the accuracy of the contrast transfer function description
        after an initial reconstruction has already been obtained. In modern
        high-resolution cryo-EM workflows, this refinement step is often
        essential for reaching the highest possible map quality because small
        optical inaccuracies can significantly affect fine structural details.

        The protocol operates by combining previously aligned particles with
        information from a postprocessed reconstruction. Using the reconstructed
        signal as a reference, the method refines imaging parameters that may
        vary across particles, micrographs, or optics groups. Biological users
        typically apply this protocol during late stages of refinement once a
        reasonably accurate consensus reconstruction has already been achieved.

        Inputs and General Workflow

        The protocol requires a set of aligned particles together with a
        postprocessing result that provides the refined maps and solvent mask.
        The solvent mask plays an important biological role because it defines
        which regions contribute signal during refinement. In practice, the mask
        should encompass the entire molecular complex while excluding excessive
        solvent noise.

        For elongated or filamentous assemblies, users sometimes employ smaller
        masks during postprocessing to improve FSC estimation. However, for CTF
        refinement it is generally preferable to use broader masks that retain
        sufficient structural signal across the full particle region. Poor mask
        design may weaken parameter estimation and reduce refinement stability.

        CTF Parameter Refinement

        The protocol allows refinement of several imaging parameters, including
        defocus, astigmatism, phase shift, and B-factor terms. These parameters
        may be optimized globally at the micrograph level or individually for
        each particle, depending on dataset quality and resolution goals.

        Per-particle refinement is especially valuable for high-resolution
        datasets because local variations in ice thickness, particle height,
        charging effects, or beam-induced motion can introduce differences
        between particles collected within the same micrograph. Refining these
        effects individually often improves the consistency of particle
        alignment and reconstruction quality.

        Defocus refinement is commonly the most impactful option in routine
        biological workflows. Astigmatism refinement becomes more important
        when optical imperfections are present, while phase-shift refinement is
        particularly relevant for datasets collected using phase plates.
        B-factor refinement can help model signal attenuation and improve the
        description of high-frequency information.

        Beam Tilt and Higher-Order Aberrations

        The protocol can estimate beam tilt and additional higher-order optical
        aberrations. These corrections are mainly relevant for datasets aiming
        at near-atomic or atomic resolution, where subtle optical distortions
        become measurable and biologically significant.

        Beam tilt estimation corrects systematic phase shifts introduced by
        imperfect microscope alignment. For datasets extending beyond moderate
        resolution ranges, correcting beam tilt can noticeably improve map
        sharpness and interpretability.

        Trefoil and fourth-order aberration estimation provide even more
        detailed optical corrections. These options are generally recommended
        only for very high-resolution datasets because lower-resolution data
        rarely contain enough signal to support stable estimation. Biological
        users should interpret these advanced corrections cautiously and avoid
        unnecessary over-parameterization in weaker datasets.

        Anisotropic Magnification Refinement

        The protocol also supports estimation of anisotropic magnification
        distortions. These distortions arise when magnification differs slightly
        along different detector axes, producing subtle geometric stretching or
        compression effects in reconstructed maps.

        Correcting anisotropic magnification is particularly important when
        combining datasets from multiple optics groups or when pursuing the
        highest possible spatial accuracy. In many workflows, users alternate
        between anisotropic magnification correction and higher-order aberration
        refinement until optical parameters become stable.

        Resolution Limits and Stability

        The protocol allows users to define the minimum resolution used during
        fitting procedures. Restricting refinement to appropriate frequency
        ranges helps avoid overfitting noise or unstable signal regions. In
        routine biological analyses, conservative resolution thresholds often
        provide more reliable results than overly aggressive refinement.

        High-resolution refinements should only be attempted when the dataset
        quality, particle number, and reconstruction resolution justify the
        additional complexity. Attempting advanced aberration estimation on
        insufficient data may produce unstable corrections without meaningful
        biological improvement.

        Outputs and Their Interpretation

        After execution, the protocol produces a refined particle set with
        updated imaging parameters. These refined particles are typically used
        in subsequent rounds of three-dimensional refinement, polishing, or
        reconstruction.

        Additional diagnostic outputs describing optical distortions and
        aberration fits may also be generated. These outputs help users inspect
        systematic microscope behavior and assess whether the estimated
        corrections are physically meaningful.

        From a biological perspective, successful CTF refinement often leads to
        sharper density maps, improved side-chain visibility, cleaner secondary
        structure features, and better interpretability of flexible or weakly
        resolved regions.

        Practical Recommendations

        In standard cryo-EM workflows, it is often advisable to begin with
        conservative refinement options such as defocus correction before moving
        toward more advanced aberration estimation. Per-particle refinement is
        generally beneficial for high-quality datasets, but advanced optical
        corrections should only be enabled when the reconstruction resolution
        supports them.

        Beam tilt and higher-order aberration refinement are most valuable for
        datasets approaching high or near-atomic resolution. Users should
        visually inspect reconstruction improvements after each refinement cycle
        rather than relying exclusively on numerical resolution estimates.

        When anisotropic magnification and higher-order aberrations are both
        suspected, iterative refinement strategies frequently provide the most
        stable results. Applying one correction may improve the accuracy of the
        next refinement stage.

        Final Perspective

        CTF refinement represents one of the key late-stage optimization steps
        in modern single-particle cryo-EM. By improving the accuracy of optical
        parameter estimation, the protocol helps maximize the structural signal
        contained within experimental particles and enables more reliable
        biological interpretation of high-resolution reconstructions.
    """
    _label = 'ctf refinement'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _initialize(self):
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'output_star': self._getExtraPath("particles_ctf_refine.star"),
            'ctf_sqlite': self._getExtraPath("ctf_analyze.sqlite"),
            'mag_obs_x': self._getExtraPath("mag_disp_x_optics-group_%(og)d.mrc"),
            'mag_obs_y': self._getExtraPath("mag_disp_y_optics-group_%(og)d.mrc"),
            'mag_fit_x': self._getExtraPath("mag_disp_x_fit_optics-group_%(og)d.mrc"),
            'mag_fit_y': self._getExtraPath("mag_disp_y_fit_optics-group_%(og)d.mrc"),
            'tetrafoil_it_fit': self._getExtraPath("aberr_delta-phase_iter-fit_optics-group_%(og)d_N-4.mrc"),
            'tetrafoil_fit': self._getExtraPath("aberr_delta-phase_lin-fit_optics-group_%(og)d_N-4.mrc"),
            'tetrafoil_residual_fit': self._getExtraPath("aberr_delta-phase_lin-fit_optics-group_%(og)d_N-4_residual.mrc"),
            'tetrafoil_obs': self._getExtraPath("aberr_delta-phase_per-pixel_optics-group_%(og)d.mrc"),
            'beamtilt_it_fit': self._getExtraPath("beamtilt_delta-phase_iter-fit_optics-group_%(og)d.mrc"),
            'beamtilt_fit': self._getExtraPath("beamtilt_delta-phase_lin-fit_optics-group_%(og)d.mrc"),
            'trefoil_it_fit': self._getExtraPath("beamtilt_delta-phase_iter-fit_optics-group_%(og)d_N-3.mrc"),
            'trefoil_fit': self._getExtraPath("beamtilt_delta-phase_lin-fit_optics-group_%(og)d_N-3.mrc"),
            'trefoil_residual_fit': self._getExtraPath("beamtilt_delta-phase_lin-fit_optics-group_%(og)d_N-3_residual.mrc"),
            'beamtilt_obs': self._getExtraPath("beamtilt_delta-phase_per-pixel_optics-group_%(og)d.mrc")
        }

        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles for local CTF '
                           'refinement.')
        form.addParam('inputPostprocess', params.PointerParam,
                      important=True,
                      label='Input Postprocess',
                      pointerClass='ProtRelionPostprocess',
                      help='Select a PostProcess job. The mask used for this '
                           'postprocessing will be applied to the unfiltered '
                           'half-maps and should encompass the entire '
                           'complex. '
                           'The resulting FSC curve will be used for weighting'
                           ' the different frequencies.\n\n'
                           'Note that for helices it is common practice to '
                           'use '
                           'a mask only encompassing the central 30% or so of '
                           'the box. This gives higher resolution estimates, '
                           'as it disregards ill-defined regions near the box '
                           'edges. However, for ctf_refine it is better to '
                           'use a mask encompassing (almost) the entire box, '
                           'as otherwise there may not be enough signal. ')

        form.addSection(label='Fit')
        form.addParam('estimateAnisoMag', params.BooleanParam,
                      default=False,
                      label='Estimate (anisotropic) magnification?',
                      help="If set to Yes, then relion_ctf_refine will also "
                           "estimate the (anisotropic) magnification per optics"
                           " group. This option cannot be done simultaneously "
                           "with higher-order aberration estimation. It's "
                           "probably best to estimate the one that is most off "
                           "first, and the other one second. It might be worth "
                           "repeating the estimation if both are off.")

        group = form.addGroup('CTF', condition='not estimateAnisoMag')

        group.addParam('doCtfFitting', params.BooleanParam, default=False,
                       condition='not estimateAnisoMag',
                       label='Perform CTF parameter fitting?',
                       help='If set to Yes, then relion_ctf_refine will be '
                            'used to estimate the selected parameters below.')

        group.addParam('fitDefocus', params.EnumParam, default=relion.FIT_NO,
                       condition='doCtfFitting and not estimateAnisoMag',
                       choices=['no', 'per-micrograph', 'per-particle'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Fit defocus?',
                       help='If set to per-particle or per-micrograph, then '
                            'relion_ctf_refine will estimate a defocus values.')

        group.addParam('fitAstig', params.EnumParam, default=relion.FIT_NO,
                       condition='doCtfFitting and not estimateAnisoMag',
                       choices=['no', 'per-micrograph', 'per-particle'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Fit astigmatism? ',
                       help="If set to per-particle or per-micrograph, then "
                            "relion_ctf_refine will estimate astigmatism.")

        group.addParam('fitBfactor', params.EnumParam, default=relion.FIT_NO,
                       condition='doCtfFitting and not estimateAnisoMag',
                       choices=['no', 'per-micrograph', 'per-particle'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Fit B-factor?',
                       help='If set to per-particle or per-micrograph, then '
                            'relion_ctf_refine will estimate B-factors that '
                            'describe the signal falloff.')

        group.addParam('fitPhaseShift', params.EnumParam, default=relion.FIT_NO,
                       condition='doCtfFitting and not estimateAnisoMag',
                       choices=['no', 'per-micrograph', 'per-particle'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Fit phase-shift? ',
                       help="If set to per-particle or per-micrograph, then "
                            "relion_ctf_refine will estimate astigmatism.")

        form.addParam('doBeamtiltEstimation', params.BooleanParam, default=False,
                      label='Estimate beamtilt?',
                      condition='not estimateAnisoMag',
                      help='If set to Yes, then relion_ctf_refine will '
                           'also estimate the beamtilt per optics group. '
                           'This option is only recommended for data sets '
                           'that extend beyond 4.5 Angstrom resolution.')
        form.addParam('doEstimateTrefoil', params.BooleanParam, default=False,
                      condition='doBeamtiltEstimation and not estimateAnisoMag',
                      label='Also estimate trefoil?',
                      help='If set to Yes, then relion_ctf_refine will also '
                           'estimate the trefoil (3-fold astigmatism) per '
                           'optics group. This option is only recommended for '
                           'data sets that extend beyond 3.5 Angstrom '
                           'resolution.')

        form.addParam('doEstimate4thOrder', params.BooleanParam, default=False,
                      label='Estimate 4th order aberrations?',
                      condition='not estimateAnisoMag',
                      help='If set to Yes, then relion_ctf_refine will also '
                           'estimate the Cs and the tetrafoil (4-fold '
                           'astigmatism) per optics group. This option is only '
                           'recommended for data sets that extend beyond 3 '
                           'Angstrom resolution.')

        form.addParam('minResolution', params.FloatParam, default=30,
                      label='Minimum resolution for fits (A)',
                      help="The minimum spatial frequency (in Angstrom) used "
                           "in the beam tilt fit.")

        form.addParam('extraParams', params.StringParam,
                      default='',
                      label='Additional arguments',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options "
                           "and/or expert use of the program")

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.refineCtfStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

        if self.doCtfFitting:
            self._insertFunctionStep(self.createGlobalInfoStep, needsGPU=False)

    def convertInputStep(self):
        inputParts = self.inputParticles.get()
        imgStar = self._getPath('input_particles.star')

        self.info("Converting set from '%s' into '%s'" %
                  (inputParts.getFileName(), imgStar))

        convert.writeSetOfParticles(inputParts, imgStar,
                                    outputDir=self._getExtraPath(),
                                    alignType=ALIGN_PROJ,
                                    fillMagnification=True)

    def refineCtfStep(self):
        args = "--i %s " % self._getPath('input_particles.star')
        args += "--o %s " % self._getExtraPath()
        inputProt = self.inputPostprocess.get()
        postStar = inputProt._getExtraPath('postprocess.star')
        args += "--f %s " % postStar
        args += "--angpix_ref %0.5f " % inputProt.solventMask.get().getSamplingRate()
        minRes = '%0.3f' % self.minResolution

        if self.estimateAnisoMag:
            args += " --fit_aniso --kmin_mag %s " % minRes
        else:
            if self.doCtfFitting:
                def _letter(option):
                    options = ['f', 'm', 'p']
                    return options[self.getAttributeValue(option)]

                args += "--fit_defocus --kmin_defocus %s " % minRes
                args += "--fit_mode %s%s%sf%s " % (_letter('fitPhaseShift'),
                                                   _letter('fitDefocus'),
                                                   _letter('fitAstig'),
                                                   _letter('fitBfactor'))

            if self.doBeamtiltEstimation:
                args += "--fit_beamtilt --kmin_tilt %s " % minRes
                if self.doEstimateTrefoil:
                    args += " --odd_aberr_max_n 3 "

            if self.doEstimate4thOrder:
                args += '--fit_aberr '

        args += "--j %d " % self.numberOfThreads

        if self.extraParams.hasValue():
            args += ' ' + self.extraParams.get()

        self._runProgram("relion_ctf_refine", args)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        outImgsFn = self._getFileName("output_star")
        imgSet.setAlignmentProj()

        mdIter = convert.Table.iterRows('particles@' + outImgsFn,
                                        key='rlnImageId', types=convert.LABELS_DICT)
        outImgSet.copyItems(imgSet,
                            updateItemCallback=self._updateItem,
                            itemDataIterator=mdIter,
                            doClone=False)
        og = OpticsGroups.fromStar(outImgsFn)
        og.toImages(outImgSet)

        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(self.inputParticles, outImgSet)

    def createGlobalInfo(self, filename):
        pwutils.cleanPath(filename)
        ctfInfo = CtfRefineGlobalInfo(filename=filename)
        ctfInfo.loadFromParticles(self.inputParticles.get(),
                                  self.outputParticles)
        return ctfInfo

    def createGlobalInfoStep(self):
        self.createGlobalInfo(self._getFileName("ctf_sqlite"))

    def _updateItem(self, particle, row):
        Reader.rowToCtf(row, particle.getCTF())
        # Reader.rowToAcquisition(self._optics[row.rlnOpticsGroup],
        #                         particle.getAcquisition())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.estimateAnisoMag:
            summary.append("Estimate anisotropic magnification: *Yes*")
        else:
            if self.doCtfFitting:
                summary.append("CTF parameter fitting: *Yes*")
                for p in ['fitPhaseShift', 'fitDefocus', 'fitAstig', 'fitBfactor']:
                    summary.append("   - %s: *%s*" % (self.getParam(p).getLabel(),
                                                      self.getEnumText(p)))

            if self.doBeamtiltEstimation:
                trefoil = '*Yes*' if self.doEstimateTrefoil else 'No'
                summary.append("Estimate beamtilt: *Yes*, trefoil: " + trefoil)

            if self.doEstimate4thOrder:
                summary.append("Estimate 4th order aberrations: *Yes*")

        return summary

    def _validate(self):
        errors = []
        return errors
