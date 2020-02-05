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

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pwem.emlib.metadata as md
from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtParticles
from pyworkflow.object import Float

import relion
import relion.convert as convert
from ..objects import CtfRefineGlobalInfo
from ..convert.metadata import Table


class ProtRelionCtfRefinement(ProtParticles):
    """ Wrapper protocol for the Relion's per-particle CTF refinement. """
    _label = 'ctf refinement'

    def _defineParams(self, form):
        form.addSection(label='Input')
        # TODO: conditions on particles?
        form.addParam('inputParticles', params.PointerParam,
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

        form.addParam('doCtfFitting', params.BooleanParam, default=True,
                      condition='not estimateAnisoMag',
                      label='Perform CTF parameter fitting?',
                      help='If set to Yes, then relion_ctf_refine will be '
                           'used to estimate the selected parameters below.')

        form.addParam('fitDefocus', params.EnumParam, default=relion.FIT_NO,
                      condition='doCtfFitting',
                      choices=['no', 'per-micrograph', 'per-particle'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Fit defocus?',
                      help='If set to per-particle or per-micrograph, then '
                           'relion_ctf_refine will estimate a defocus values.')

        form.addParam('fitAstig', params.EnumParam, default=relion.FIT_NO,
                      condition='doCtfFitting',
                      choices=['no', 'per-micrograph', 'per-particle'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Fit astigmatism? ',
                      help="If set to per-particle or per-micrograph, then "
                           "relion_ctf_refine will estimate astigmatism.")

        form.addParam('fitBfactor', params.EnumParam, default=relion.FIT_NO,
                      condition='doCtfFitting',
                      choices=['no', 'per-micrograph', 'per-particle'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Fit B-factor?',
                      help='If set to per-particle or per-micrograph, then '
                           'relion_ctf_refine will estimate B-factors that '
                           'describe the signal falloff.')

        form.addParam('fitPhaseShift', params.EnumParam, default=relion.FIT_NO,
                      condition='doCtfFitting',
                      choices=['no', 'per-micrograph', 'per-particle'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Fit phase-shift? ',
                      help="If set to per-particle or per-micrograph, then "
                           "relion_ctf_refine will estimate astigmatism.")

        form.addParam('doBeamtiltEstimation', params.BooleanParam, default=True,
                      label='Estimate beamtilt?',
                      help='If set to Yes, then relion_ctf_refine will '
                           'also estimate the beamtilt per optics group. '
                           'This option is only recommended for data sets '
                           'that extend beyond 4.5 Angstrom resolution.')
        form.addParam('doEstimateTrefoil', params.BooleanParam, default=False,
                      condition='doBeamtiltEstimation',
                      label='Also estimate trefoil?',
                      help='If set to Yes, then relion_ctf_refine will also '
                           'estimate the trefoil (3-fold astigmatism) per '
                           'optics group. This option is only recommended for '
                           'data sets that extend beyond 3.5 Angstrom '
                           'resolution.')

        form.addParam('doEstimate4thOrder', params.BooleanParam, default=False,
                      label='Estimate 4th order aberrations?',
                      help='If set to Yes, then relion_ctf_refine will also '
                           'estimate the Cs and the tetrafoil (4-fold '
                           'astigmatism) per optics group. This option is only '
                           'recommended for data sets that extend beyond 3 '
                           'Angstrom resolution.')

        form.addParam('minResolution', params.FloatParam, default=30,
                      label='Minimum resolution for fits (A)',
                      help="The minimum spatial frequency (in Angstrom) used "
                           "in the beamtilt fit.")

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('refineCtfStep')
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep('createGlobalInfoStep')

    def convertInputStep(self):
        inputParts = self.inputParticles.get()
        imgStar = self._getPath('input_particles.star')

        self.info("Converting set from '%s' into '%s'" %
                  (inputParts.getFileName(), imgStar))

        convert.writeSetOfParticles(inputParts, imgStar,
                                    self._getExtraPath(),
                                    alignType=ALIGN_PROJ,
                                    fillMagnification=True,
                                    fillRandomSubset=True)

    def _getInputVolumes(self, postStar):
        """ Parse the input volumes: halves and mask
        from the postprocess.star file. """
        table = Table(fileName=postStar, tableName='general')
        row = table[0]
        return (row.rlnUnfilteredMapHalf1,
                row.rlnUnfilteredMapHalf2,
                row.rlnMaskName)

    def refineCtfStep(self):
        """
        `which relion_ctf_refine` --i particles.star --f postprocess.star --o CtfRefine/job001/

        --fit_beamtilt --kmin_tilt 30 --j 1
        --fit_beamtilt --kmin_tilt 30 --odd_aberr_max_n 3 --j 1
        --fit_beamtilt --kmin_tilt 30 --odd_aberr_max_n 3 --fit_aberr
        --fit_beamtilt --kmin_tilt 3

        --fit_aniso --kmin_mag 30 --j 1

        --fit_defocus --kmin_defocus 30 --fit_mode fmfff --j 4
        --fit_defocus --kmin_defocus 30 --fit_mode fpfff --j 4
        --fit_defocus --kmin_defocus 30 --fit_mode ffmff --j 4
        --fit_defocus --kmin_defocus 30 --fit_mode ffpff --j 4

        --fit_defocus --kmin_defocus 30 --fit_mode pppfp --j 4  --pipeline_control CtfRefine/job001/
        --fit_defocus --kmin_defocus 30 --fit_mode mppfp --j 4  --pipeline_control CtfRefine/job001/

        	// Always either do anisotropic magnification, or CTF,tilt-odd,even
        if (joboptions["do_aniso_mag"].getBoolean())
        {
            command += " --fit_aniso";
            command += " --kmin_mag " + joboptions["minres"].getString();
        }
        else
        {
            if (joboptions["do_ctf"].getBoolean())
            {
                command += " --fit_defocus --kmin_defocus " + joboptions["minres"].getString();
                std::string fit_options = "";

                fit_options += getStringFitOption(joboptions["do_phase"].getString());
                fit_options += getStringFitOption(joboptions["do_defocus"].getString());
                fit_options += getStringFitOption(joboptions["do_astig"].getString());
                fit_options += "f"; // always have Cs refinement switched off
                fit_options += getStringFitOption(joboptions["do_bfactor"].getString());

                command += " --fit_mode " + fit_options;
            }

            // do not allow anisotropic magnification to be done simultaneously with higher-order aberrations
            if (joboptions["do_tilt"].getBoolean())
            {
                command += " --fit_beamtilt";
                command += " --kmin_tilt " + joboptions["minres"].getString();

                if (joboptions["do_trefoil"].getBoolean())
                {
                    command += " --odd_aberr_max_n 3";
                }
            }

            if (joboptions["do_4thorder"].getBoolean())
            {
                command += " --fit_aberr";
            }
        }
        """

        args = "--i %s " % self._getPath('input_particles.star')
        args += "--o %s " % self._getExtraPath()
        inputProt = self.inputPostprocess.get()
        postStar = inputProt._getExtraPath('postprocess.star')
        args += "--f %s " % postStar
        args += "--angpix_ref %0.3f " % inputProt.solventMask.get().getSamplingRate()

        # TODO: Check why we were using -m1 and -m2 options here
        # Maybe not needed in R3.1???
        #postVols = convert.getVolumesFromPostprocess(postStar)
        #args += "--m1 %s --m2 %s --mask %s " % postVols


        minRes = '%0.3f' % self.minResolution

        # New command line string taken from here:
        # https://github.com/3dem/relion/blob/a5d691e8a9507c1efcce2810ba74f5fac6d2a098/src/pipeline_jobs.cpp#L5050
        if self.estimateAnisoMag:
            args += " --fit_aniso --kmin_mag %s" % minRes
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
        prog = "relion_ctf_refine" + ("_mpi" if self.numberOfMpi > 1 else "")
        self.runJob(prog, args)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        outImgsFn = self.fileWithRefinedCTFName()
        imgSet.setAlignmentProj()
        rowIterator = md.iterRows('particles@' + outImgsFn,
                                  sortByLabel=md.RLN_IMAGE_ID)
        outImgSet.copyItems(imgSet,
                            updateItemCallback=self._updateItemCtfBeamTilt,
                            itemDataIterator=rowIterator)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

    def createGlobalInfo(self, filename):
        pwutils.cleanPath(filename)
        ctfInfo = CtfRefineGlobalInfo(filename=filename)
        ctfInfo.loadFromParticles(self.inputParticles.get(),
                                  self.outputParticles)
        return ctfInfo

    def createGlobalInfoStep(self):
        self.createGlobalInfo(self.fileWithAnalyzeInfo())

    def _updateItemCtfBeamTilt(self, particle, row):
        particle.setCTF(convert.rowToCtfModel(row))
        # TODO: Add other field from the .star file when other options?
        # check if beamtilt is available and save it
        if row.hasLabel('rlnBeamTiltX'):
            particle._rlnBeamTiltX = Float(row.getValue('rlnBeamTiltX', 0))
            particle._rlnBeamTiltY = Float(row.getValue('rlnBeamTiltY', 0))

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        return errors

    def fileWithRefinedCTFName(self):
        return self._getExtraPath('particles_ctf_refine.star')

    def fileWithPhaseDifferenceName(self):
        return self._getExtraPath(
            'beamtilt_delta-phase_per-pixel_class_0.mrc:mrc')

    def fileWithModelFitterName(self):
        return self._getExtraPath(
            'beamtilt_delta-phase_lin-fit_class_0.mrc:mrc')

    def fileWithAnalyzeInfo(self):
        return self._getExtraPath('ctf_analyze.sqlite')
