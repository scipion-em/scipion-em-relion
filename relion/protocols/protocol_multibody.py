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
from enum import Enum
from emtable import Table

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.objects import Volume, Float, SetOfVolumes, SetOfParticles
from pwem.protocols import ProtAnalysis3D
from pwem.constants import ALIGN_PROJ

from relion import Plugin
import relion.convert as convert
from ..constants import ANGULAR_SAMPLING_LIST, LABELS_DICT, PARTICLE_EXTRA_LABELS
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputVolumes = SetOfVolumes
    outputParticles = SetOfParticles


class ProtRelionMultiBody(ProtAnalysis3D, ProtRelionBase):
    """
    Relion protocol for multi-body refinement.

    This approach models flexible complexes as a user-defined number of rigid
    bodies that move independently of each other.
    Using separate focused refinements with iteratively improved partial
    signal subtraction, improved reconstructions are generated for
    each of the defined bodies.

    Moreover, using PCA on the relative orientations of the bodies
    over all particle images in the data set, we generate movies that describe
    the most important motions in the data.
    """
    _label = '3D multi-body'
    _devStatus = PROD
    _possibleOutputs = outputs
    IS_CLASSIFY = False
    IS_3D_MB = True
    PREFIXES = ['half1_', 'half2_']

    def _initialize(self):
        """ This function is meant to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

    def _getInputPath(self, *paths):
        return self._getPath('input', *paths)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineConstants()

        form.addSection(label='Input')
        form.addParam('doContinue', params.BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous '
                           'MultiBody protocol and most of the input parameters '
                           'will be taken from it.')
        form.addParam('protRefine', params.PointerParam,
                      condition='not doContinue',
                      pointerClass="ProtRefine3D",
                      label='Consensus refinement protocol',
                      help='Select any previous refinement protocol from '
                           'where to run the multi-body refinement. '
                           'The output volume will be used and some '
                           'parameters from the optimiser.star file. ')
        # FIXME: Find an easy way to avoid input a file here
        form.addParam('bodyStarFile', params.FileParam,
                      condition='not doContinue',
                      label='Body STAR file',
                      help='Provide the STAR file with all information '
                           'about the bodies to be used in multi-body '
                           'refinement. An example for a three-body '
                           'refinement would look like this:\n\n'
                           'data_\n'
                           'loop_\n'
                           '_rlnBodyMaskName\n'
                           '_rlnBodyRotateRelativeTo\n'
                           '_rlnBodySigmaAngles\n'
                           '_rlnBodySigmaOffset\n'
                           'large_body_mask.mrc 2 10 2\n'
                           'small_body_mask.mrc 1 10 2\n'
                           'head_body_mask.mrc 2 10 2\n\n'
                           ''

                           """
 Where each data line represents a different body, and:
 - rlnBodyMaskName contains the name of a soft-edged mask with values in [0,1] that define the body; the mask name should be relative to the project folder;
 - rlnBodyRotateRelativeTo defines relative to which other body this body rotates (first body is number 1);
 - rlnBodySigmaAngles and _rlnBodySigmaOffset are the standard deviations (widths) of Gaussian priors on the consensus rotations and translations;

 Optionally, there can be a fifth column with _rlnBodyReferenceName. Entries can be 'None' (without the ''s) or the name of a MRC map with an initial reference for that body. In case the entry is None, the reference will be taken from the density in the consensus refinement.

Also note that larger bodies should be above smaller bodies in the STAR file. For more information, see the multi-body paper.
                           """)

        form.addParam('recSubtractedBodies', params.BooleanParam, default=True,
                      condition='not doContinue',
                      label='Reconstruct subtracted bodies?',
                      help='If set to Yes, then the reconstruction of each of '
                           'the bodies will use the subtracted images. This '
                           'may give useful insights about how well the '
                           'subtraction worked. If set to No, the original '
                           'particles are used for reconstruction (while the '
                           'subtracted ones are still used for alignment). '
                           'This will result in fuzzy densities for bodies '
                           'outside the one used for refinement.')
        if Plugin.IS_GT50():
            form.addParam('useBlush', params.BooleanParam, default=False,
                          label='Use Blush regularisation?',
                          help='If set to Yes, relion_refine will use a neural '
                               'network to perform regularisation by denoising '
                               'at every iteration, instead of the standard '
                               'smoothness regularisation.')
        form.addParam('continueRun', params.PointerParam,
                      pointerClass='ProtRelionMultiBody',
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('continueIter', params.StringParam, default='last',
                      condition='doContinue',
                      label='Continue from iteration',
                      help='Select from which iteration do you want to '
                           'continue. If you use *last*, then the last '
                           'iteration will be used. Otherwise, a valid '
                           'iteration number should be provided.')

        group = form.addGroup('Auto-Sampling',
                              condition='not doContinue')
        group.addParam('initialAngularSampling', params.EnumParam, default=4,
                       condition='not doContinue',
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
                       condition='not doContinue',
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
                       condition='not doContinue',
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
                      label='Additional arguments',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. They will "
                           "be appended to the relion_refine command.")

        form.addParallelSection(threads=1, mpi=3)
    
    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()
        if self.doContinue:
            objId = self.continueRun.get().getObjId()
        else:
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
        if self.doContinue:
            bodyFn = self.continueRun.get().bodyStarFile.get()
        else:
            bodyFn = self.bodyStarFile.get()

        pwutils.copyFile(bodyFn, self._getExtraPath('input_body.star'))

    def multibodyRefineStep(self, args):
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.items()])
        if self.extraParams.hasValue():
            params += ' ' + self.extraParams.get()

        self._runProgram('relion_refine', params)

    def flexAnalysisStep(self, args):
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.items()])
        # use runJob since MPI is not allowed
        self.runJob('relion_flex_analyse', params, numberOfMpi=1,
                    numberOfThreads=1)

    def createOutputStep(self):
        protRefine = self._getProtRefine()
        if self.doContinue:
            # get original 3D refine protocol
            protRefine = protRefine.protRefine.get()

        # get refine 3d output parts pointer
        inputPartsSet = protRefine.outputParticles
        sampling = inputPartsSet.getSamplingRate()
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(sampling)

        self._loadVolsInfo()

        for item in range(1, self._getNumberOfBodies() + 1):
            vol = Volume()
            self._updateVolume(item, vol)
            volumes.append(vol)

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(inputPartsSet)
        self._fillDataFromIter(inputPartsSet, outImgSet, self._lastIter())

        self._defineOutputs(**{outputs.outputVolumes.name: volumes})
        self._defineSourceRelation(protRefine.outputVolume, volumes)
        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(inputPartsSet, outImgSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if self.doContinue:
            continueProtocol = self.continueRun.get()
            if (continueProtocol is not None and
                    continueProtocol.getObjId() == self.getObjId()):
                errors.append('In Scipion you must create a new Relion run '
                              'and select the continue option rather than '
                              'select continue from the same run.\n')

            continueProtocol._initialize()
            lastIter = continueProtocol._lastIter()

            if self.continueIter.get() == 'last':
                continueIter = lastIter
            else:
                continueIter = int(self.continueIter.get())

            if continueIter > lastIter:
                errors.append("You can continue only from the iteration %01d or less" % lastIter)
        else:
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
        self._initialize()
        lastIter = self._lastIter()

        if lastIter is not None:
            iterMsg = 'Iteration %d' % lastIter
        else:
            iterMsg = 'No iteration finished yet.'

        summary = [iterMsg]

        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
            it = self._lastIter() or -1
            if it >= 1:
                table = Table(fileName=self._getFileName('half1_model', iter=it),
                              tableName='model_general')
                row = table[0]
                resol = float(row.rlnCurrentResolution)
                summary.append("Current resolution: *%0.2f A*" % resol)
        else:
            table = Table(fileName=self._getFileName('modelFinal'),
                          tableName='model_general')
            row = table[0]
            resol = float(row.rlnCurrentResolution)
            summary.append("Final resolution: *%0.2f A*" % resol)

        return summary

    # -------------------------- UTILS functions ------------------------------
    def _loadVolsInfo(self):
        """ Read some information about the produced Relion bodies
        from the *model.star file.
        """
        self._volsInfo = {}
        mdTable = Table(fileName=self._getFileName('modelFinal'),
                        tableName='model_bodies')

        for body, row in enumerate(mdTable):
            self._volsInfo[body + 1] = row

    def _getNumberOfBodies(self):
        table = Table(fileName=self._getExtraPath("input_body.star"))
        return int(table.size())

    def _updateVolume(self, bodyNum, item):
        item.setFileName(self._getFileName('finalvolume_mbody', ref3d=bodyNum))
        half1 = self._getFileName('final_half1_volume_mbody', ref3d=bodyNum)
        half2 = self._getFileName('final_half2_volume_mbody', ref3d=bodyNum)
        item.setHalfMaps([half1, half2])

        row = self._volsInfo[bodyNum]
        item._rlnAccuracyRotations = Float(row.rlnAccuracyRotations)
        item._rlnAccuracyTranslationsAngst = Float(row.rlnAccuracyTranslationsAngst)

    def _getRefineArgs(self):
        """ Define all parameters to run relion_refine. """
        args = {'--o': self._getExtraPath('relion')}
        protRefine = self._getProtRefine()
        protRefine._initialize()

        if self.doContinue:
            continueIter = self._getContinueIter()
            fnOptimiser = protRefine._getFileName('optimiser',
                                                  iter=continueIter)

            if protRefine.recSubtractedBodies:
                args['--reconstruct_subtracted_bodies'] = ''

        else:
            fnOptimiser = protRefine._getOptimiserFile()
            args.update({
                '--multibody_masks': self._getExtraPath('input_body.star'),
                '--solvent_correct_fsc': '',
                '--oversampling': 1,
                '--pad': 1 if self.skipPadding else 2,
                '--healpix_order': self.initialAngularSampling.get(),
                '--auto_local_healpix_order': self.initialAngularSampling.get(),
                '--offset_range': self.initialOffsetRange.get(),
                '--offset_step': self.initialOffsetStep.get()
            })

            if self.recSubtractedBodies:
                args['--reconstruct_subtracted_bodies'] = ''

            if Plugin.IS_GT50() and self.useBlush:
                args['--blush'] = ''

            # Due to Relion bug we create a fake mask from previous refinement protocol
            # it's not used by multi-body
            if protRefine.referenceMask.hasValue():
                table = Table(fileName=fnOptimiser,
                              tableName='optimiser_general',
                              types=LABELS_DICT)
                maskFn = table[0].rlnSolventMaskName
                bodyFn = self.bodyStarFile.get()
                maskBody1 = Table(fileName=bodyFn)[0].rlnBodyMaskName
                os.makedirs(os.path.dirname(maskFn), exist_ok=True)
                pwutils.createAbsLink(os.path.abspath(maskBody1), maskFn)

        args['--continue'] = fnOptimiser

        self._setComputeArgs(args)

        return args

    def _getAnalyseArgs(self):
        args = {
            '--PCA_orient': '',
            '--model': self._getFileName('modelFinal'),
            '--data': self._getFileName('dataFinal'),
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

    def _getProtRefine(self):
        return self.continueRun.get() if self.doContinue else self.protRefine.get()

    def _fillDataFromIter(self, inputSet, outSet, iteration):
        outImgsFn = self._getFileName('data', iter=iteration)
        outSet.setAlignmentProj()
        self.reader = convert.createReader(alignType=ALIGN_PROJ,
                                           pixelSize=outSet.getSamplingRate())

        mdIter = Table.iterRows('particles@' + outImgsFn, key='rlnImageId',
                                types=convert.LABELS_DICT)
        outSet.copyItems(inputSet, doClone=False,
                         updateItemCallback=self._updateParticle,
                         itemDataIterator=mdIter)

    def _updateParticle(self, particle, row):
        self.reader.setParticleTransform(particle, row)

        if getattr(self, '__updatingFirst', True):
            self.reader.createExtraLabels(particle, row, PARTICLE_EXTRA_LABELS)
            self.__updatingFirst = False
        else:
            self.reader.setExtraLabels(particle, row)
