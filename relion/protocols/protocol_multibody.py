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
from emtable import Table

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.objects import Volume, Float
from pwem.protocols import ProtAnalysis3D

import relion.convert as convert
from ..constants import ANGULAR_SAMPLING_LIST
from .protocol_base import ProtRelionBase


class ProtRelionMultiBody(ProtAnalysis3D, ProtRelionBase):
    """
    Relion protocol for multi-body refinement.

    This approach models flexible complexes as a user-defined number of rigid
    bodies that move independently from each other.
    Using separate focused refinements with iteratively improved partial
    signal subtraction, improved reconstructions are generated for
    each of the defined bodies.

    Moreover, using PCA on the relative orientations of the bodies
    over all particle images in the data set, we generate movies that describe
    the most important motions in the data.
    """
    _label = '3D multi-body'
    _devStatus = PROD

    def _getInputPath(self, *paths):
        return self._getPath('input', *paths)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'finalModel': self._getExtraPath("run_model.star"),
            'input': self._getExtraPath("input_body.star"),
            'body': self._getExtraPath("run_body%(body)03d.mrc"),
            'half1': self._getExtraPath("run_half1_body%(body)03d_unfil.mrc"),
            'half2': self._getExtraPath("run_half2_body%(body)03d_unfil.mrc")
        }
        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineConstants()

        form.addSection(label='Input')
        form.addParam('protRefine', params.PointerParam,
                      pointerClass="ProtRefine3D",
                      label='Consensus refinement protocol',
                      help='Select any previous refinement protocol from '
                           'where to run the multi-body refinement. '
                           'The output volume will be used and some '
                           'parameters from the optimiser.star file. ')
        # FIXME: Find an easy way to avoid input a file here
        form.addParam('bodyStarFile', params.FileParam,
                      label='Body STAR file',
                      help='Provide the STAR file with all information '
                           'about the bodies to be used in multi-body '
                           'refinement. An example for a three-body '
                           'refinement would look like this:\n'
                           'data_\n'
                           'loop_\n'
                           '_rlnBodyMaskName\n'
                           '_rlnBodyRotateRelativeTo\n'
                           '_rlnBodySigmaAngles\n'
                           '_rlnBodySigmaOffset\n'
                           'large_body_mask.mrc 2 10 2\n'
                           'small_body_mask.mrc 1 10 2 \n'
                           'head_body_mask.mrc 2 10 2 \n')

        """
 Where each data line represents a different body, and:
  - rlnBodyMaskName contains the name of a soft-edged mask with values in [0,1] that define the body;
 - rlnBodyRotateRelativeTo defines relative to which other body this body rotates (first body is number 1);
 - rlnBodySigmaAngles and _rlnBodySigmaOffset are the standard deviations (widths) of Gaussian priors on the consensus rotations and translations;

 Optionally, there can be a fifth column with _rlnBodyReferenceName. Entries can be 'None' (without the ''s) or the name of a MRC map with an initial reference for that body. In case the entry is None, the reference will be taken from the density in the consensus refinement.

Also note that larger bodies should be above smaller bodies in the STAR file. For more information, see the multi-body paper.')
        """

        form.addParam('recSubtractedBodies', params.BooleanParam, default=True,
                      label='Reconstruct subtracted bodies?',
                      help='If set to Yes, then the reconstruction of each of '
                           'the bodies will use the subtracted images. This '
                           'may give useful insights about how well the '
                           'subtraction worked. If set to No, the original '
                           'particles are used for reconstruction (while the '
                           'subtracted ones are still used for alignment). '
                           'This will result in fuzzy densities for bodies '
                           'outside the one used for refinement.')

        group = form.addGroup('Auto-Sampling')
        group.addParam('initialAngularSampling', params.EnumParam, default=4,
                       choices=ANGULAR_SAMPLING_LIST,
                       label='Initial angular sampling (deg)',
                       help='There are only a few discrete angular samplings'
                            ' possible because we use the HealPix library to'
                            ' generate the sampling of the first two Euler '
                            'angles on the sphere. The samplings are '
                            'approximate numbers and vary slightly over '
                            'the sphere. \n\n'
                            'Note that this will only be the value for the '
                            'first few iteration(s): the sampling rate will '
                            'be increased automatically after that.')
        group.addParam('initialOffsetRange', params.FloatParam, default=3,
                       label='Initial offset range (pix)',
                       help='Probabilities will be calculated only for '
                            'translations in a circle with this radius (in '
                            'pixels). The center of this circle changes at '
                            'every iteration and is placed at the optimal '
                            'translation for each image in the previous '
                            'iteration. \n\n'
                            'Note that this will only be the value for the '
                            'first few iteration(s): the sampling rate will '
                            'be increased automatically after that.')
        group.addParam('initialOffsetStep', params.FloatParam, default=0.75,
                       label='Initial offset step (pix)',
                       help='Translations will be sampled with this step-size '
                            '(in pixels). Translational sampling is also done '
                            'using the adaptive approach. Therefore, if '
                            'adaptive=1, the translations will first be '
                            'evaluated on a 2x coarser grid. \n\n'
                            'Note that this will only be the value for the '
                            'first few iteration(s): the sampling rate will '
                            'be increased automatically after that.')

        form.addSection(label='Analyse')

        form.addParam('runFlexAnalysis', params.BooleanParam, default=True,
                      label='Run flexibility analysis?',
                      help='If set to Yes, after the multi-body refinement has '
                           'completed, a PCA analysis will be run on the '
                           'orientations all all bodies in the data set. This '
                           'can be set to No initially, and then the job can '
                           'be continued afterwards to only perform this '
                           'analysis.')
        form.addParam('numberOfEigenvectors', params.IntParam, default=3,
                      condition='runFlexAnalysis',
                      label='Number of eigenvector movies:',
                      help='Series of ten output maps will be generated along '
                           'this many eigenvectors. These maps can be opened '
                           'as a "Volume Series" in UCSF Chimera, and then '
                           'displayed as a movie. They represent the principal '
                           'motions in the particles.')
        form.addParam('selectByEigenvalues', params.BooleanParam, default=False,
                      condition='runFlexAnalysis',
                      label='Select particles based on eigenvalues?',
                      help='If set to Yes, a particles.star file is written '
                           'out with all particles that have the below '
                           'indicated eigenvalue in the selected range.')
        form.addParam('selectEigenvalueNumber', params.IntParam, default=1,
                      condition='runFlexAnalysis and selectByEigenvalues',
                      label='Select on eigenvalue:',
                      help='This is the number of the eigenvalue to be used '
                           'in the particle subset selection '
                           '(start counting at 1).')
        line = form.addLine('Eigenvalue',
                            condition='runFlexAnalysis and selectByEigenvalues',
                            help='Minimum and maximum values for the selected '
                                 'eigenvalue; only particles with the selected '
                                 'eigenvalue within that range (min, max) will '
                                 'be included in the output particles.star file.')
        line.addParam('minEigenvalue', params.IntParam, default=-999, label='min')
        line.addParam('maxEigenvalue', params.IntParam, default=999, label='max')

        form.addSection('Compute')
        self._defineComputeParams(form)
        form.addParam('extraParams', params.StringParam,
                      default='',
                      label='Additional parameters',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. ")

        form.addParallelSection(threads=1, mpi=3)
    
    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        objId = self.protRefine.get().getObjId()
        self._insertFunctionStep('convertInputStep', objId)
        self._insertFunctionStep('multibodyRefineStep',
                                 self._getRefineArgs())
        if self.runFlexAnalysis:
            self._insertFunctionStep('flexAnalysisStep',
                                     self._getAnalyseArgs())
        self._insertFunctionStep('createOutputStep')
    
    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self, protId):
        self.info("Relion version:")
        self.runJob("relion_refine --version", "", numberOfMpi=1)

        pwutils.copyFile(self.bodyStarFile.get(),
                         self._getExtraPath('input_body.star'))

    def multibodyRefineStep(self, args):
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.items()])
        self._runProgram('relion_refine', params)

    def flexAnalysisStep(self, args):
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.items()])
        # use runJob since MPI is not allowed
        self.runJob('relion_flex_analyse', params)

    def createOutputStep(self):
        protRefine = self.protRefine.get()
        sampling = protRefine._getInputParticles().getSamplingRate()
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(sampling)

        self._loadVolsInfo()

        for item in range(1, self._getNumberOfBodies() + 1):
            vol = Volume()
            self._updateVolume(item, vol)
            volumes.append(vol)

        self._defineOutputs(outputVolumes=volumes)
        vol = self.protRefine.get().outputVolume
        self._defineSourceRelation(vol, volumes)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        bodyFn = self.bodyStarFile.get()
        if not os.path.exists(bodyFn):
            errors.append("Input body star file %s does not exist." % bodyFn)
        else:
            table = Table(fileName=bodyFn)
            missing = []
            for row in table:
                if not os.path.exists(row.rlnBodyMaskName):
                    missing.append(row.rlnBodyMaskName)
                ref = getattr(row, 'rlnBodyReferenceName', 'None')
                if ref != 'None' and not os.path.exists(ref):
                    missing.append(ref)
            if missing:
                errors.append("Missing files from input star file: ")
                for f in missing:
                    errors.append(" - %s" % f)
        return errors

    def _citations(self):
        return ['Nakane2018']

    def _summary(self):
        summary = []
        return summary

    # -------------------------- UTILS functions ------------------------------
    def _loadVolsInfo(self):
        """ Read some information about the produced Relion bodies
        from the *model.star file.
        """
        self._volsInfo = {}
        mdTable = Table(fileName=self._getFileName('finalModel'),
                        tableName='model_bodies')

        for body, row in enumerate(mdTable):
            self._volsInfo[body + 1] = row

    def _getNumberOfBodies(self):
        table = Table(fileName=self._getFileName('input'))
        return int(table.size())

    def _updateVolume(self, bodyNum, item):
        item.setFileName(self._getFileName('body', body=bodyNum))
        half1 = self._getFileName('half1', body=bodyNum)
        half2 = self._getFileName('half2', body=bodyNum)
        item.setHalfMaps([half1, half2])

        row = self._volsInfo[bodyNum]
        # TODO: check if the following makes sense for Volume
        item._rlnAccuracyRotations = Float(row.rlnAccuracyRotations)
        item._rlnAccuracyTranslationsAngst = Float(row.rlnAccuracyTranslationsAngst)

    def _getRefineArgs(self):
        """ Define all parameters to run relion_refine.
        """
        protRefine = self.protRefine.get()
        protRefine._initialize()
        fnOptimiser = protRefine._getFileName('optimiser',
                                              iter=protRefine._lastIter())
        healpix = self.initialAngularSampling.get()

        args = {
            '--continue': fnOptimiser,
            '--multibody_masks': self._getExtraPath('input_body.star'),
            '--o': self._getExtraPath('run'),
            '--solvent_correct_fsc': '',
            '--oversampling': 1,
            '--pad': 1 if self.skipPadding else 2,
            '--healpix_order': healpix,
            '--auto_local_healpix_order': healpix,
            '--offset_range': self.initialOffsetRange.get(),
            '--offset_step': self.initialOffsetStep.get()
        }

        if self.skipGridding:
            args['--skip_gridding'] = ''

        # restore mask for previous refinement protocol
        # it's not used by multibody but required by relion cmd
        if protRefine.referenceMask.hasValue():
            tmp = protRefine._getTmpPath()
            newDim = protRefine._getInputParticles().getXDim()
            newPix = protRefine._getInputParticles().getSamplingRate()
            convert.convertMask(protRefine.referenceMask.get(), tmp,
                                newPix, newDim)

        self._setComputeArgs(args)

        if self.recSubtractedBodies:
            args['--reconstruct_subtracted_bodies'] = ''

        return args

    def _getAnalyseArgs(self):
        args = {
            '--PCA_orient': '',
            '--model': self._getExtraPath('run_model.star'),
            '--data': self._getExtraPath('run_data.star'),
            '--bodies': self._getExtraPath('input_body.star'),
            '--o': self._getExtraPath('analyse'),
            '--do_maps': '',
            '--k': self.numberOfEigenvectors.get()
        }

        if self.selectByEigenvalues:
            args.update({
                '--select_eigenvalue': self.selectEigenvalueNumber.get(),
                '--select_eigenvalue_min': self.minEigenvalue.get(),
                '--select_eigenvalue_max': self.maxEigenvalue.get()
            })

        return args
