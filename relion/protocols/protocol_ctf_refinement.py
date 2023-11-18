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
    """ Wrapper protocol for the Relion's CTF refinement. """
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
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('refineCtfStep')
        self._insertFunctionStep('createOutputStep')

        if self.doCtfFitting:
            self._insertFunctionStep('createGlobalInfoStep')

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
