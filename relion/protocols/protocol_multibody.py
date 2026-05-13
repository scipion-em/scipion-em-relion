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

    AI Generated:

    Multi-Body Refinement (ProtRelionMultiBody) — User Manual
        Overview

        The Multi-Body Refinement protocol is designed to analyze structural
        flexibility in cryo-EM reconstructions by dividing a macromolecular
        complex into several rigid regions, referred to as bodies, that can
        move independently relative to one another. This strategy is especially
        useful for large assemblies, molecular machines, or complexes that
        contain flexible domains whose motion cannot be adequately represented
        by a single rigid reconstruction.

        In biological practice, many complexes exhibit continuous conformational
        variability rather than existing in only a few discrete states. Standard
        consensus refinement often averages these motions together, leading to
        blurred densities and reduced local resolution. Multi-body refinement
        addresses this limitation by refining the orientations and positions of
        independently defined bodies while preserving their internal structure.
        The result is improved local reconstructions together with a quantitative
        description of relative body motions.

        Inputs and Biological Context

        The protocol starts from a previously refined cryo-EM reconstruction and
        its associated particle set. This consensus refinement provides the
        global orientation framework from which flexible motion is analyzed.
        The user must additionally define a set of body masks describing the
        regions that should move independently.

        Each body mask should correspond to a structurally meaningful region of
        the complex. Typical examples include mobile domains, rotating subunits,
        flexible heads, peripheral arms, or independently moving membrane and
        cytosolic regions. The masks should overlap smoothly and avoid abrupt
        edges because poorly designed masks may introduce instability or
        artificial discontinuities in the refinement.

        From a biological perspective, the definition of bodies is the most
        important conceptual step in the workflow. Bodies should represent
        coherent structural units that are expected to behave approximately as
        rigid objects. Excessively small bodies may become unstable during
        refinement, whereas overly large bodies may fail to capture the relevant
        flexibility.

        Body Definition and STAR File Organization

        The protocol uses a STAR file describing all bodies involved in the
        refinement. Each entry specifies the mask associated with a body together
        with information about relative motion priors and optional initial
        references.

        The ordering of bodies is biologically meaningful because larger and
        structurally dominant regions are generally expected to appear before
        smaller flexible regions. Relative rotation definitions establish how
        motions are interpreted with respect to neighboring bodies. Gaussian
        priors on rotational and translational variability help stabilize the
        refinement by restricting unrealistic motions.

        In practice, conservative priors are often appropriate for relatively
        rigid complexes, while broader priors may be beneficial when large-scale
        conformational changes are expected. However, excessively permissive
        priors may lead to unstable or biologically implausible motions.

        Focused Refinement and Signal Subtraction

        During refinement, the protocol iteratively improves the reconstruction
        of each body using focused alignment strategies and partial signal
        subtraction. This allows the refinement of a selected body while reducing
        interference from the remaining regions of the complex.

        Biologically, this is particularly valuable when small or flexible
        regions are poorly resolved in the consensus map. By concentrating the
        alignment on the body of interest, the protocol often recovers secondary
        structure elements or conformational details that were previously hidden
        by averaging.

        The protocol also provides the option to reconstruct bodies either from
        signal-subtracted particles or from the original particle images.
        Reconstructions based on subtracted particles may better isolate the
        body of interest and help evaluate subtraction quality, whereas using
        the original particles preserves the global molecular context but may
        introduce blurred densities outside the refined region.

        Sampling and Refinement Strategy

        Multi-body refinement relies on adaptive angular and translational
        searches to optimize body orientations throughout the refinement process.
        Initial angular sampling determines the coarseness of the orientation
        search during the early stages, while offset ranges and translational
        steps control the exploration of positional variability.

        In most biological workflows, moderate initial sampling values provide
        a good balance between robustness and computational efficiency. Complexes
        exhibiting large conformational variability may require broader searches,
        whereas systems already close to convergence can benefit from finer
        refinement parameters.

        The protocol progressively refines the sampling strategy during
        iterations, improving precision as convergence is approached. This
        adaptive behavior is especially useful for heterogeneous systems where
        the initial relative orientations between bodies may not be accurately
        known.

        Flexibility Analysis and Principal Component Motions

        One of the major strengths of the protocol is its ability to analyze
        continuous conformational variability after refinement. The protocol
        performs principal component analysis on the relative orientations of
        all bodies across the particle dataset in order to identify dominant
        collective motions.

        The resulting eigenvectors describe the principal modes of structural
        variability within the sample. For each selected eigenvector, the
        protocol generates a series of reconstructed maps that can be visualized
        sequentially as movies. These animations provide an intuitive
        representation of biologically relevant motions such as domain opening,
        rotational rearrangements, hinge bending, or coordinated subunit
        movements.

        From a biological interpretation standpoint, the principal motions often
        correspond to functional transitions associated with ligand binding,
        catalytic cycles, transport mechanisms, or allosteric regulation.
        However, users should remember that principal components represent the
        dominant variance in the dataset and may combine multiple underlying
        physical processes.

        Particle Selection Based on Flexibility

        The protocol also supports selection of particle subsets according to
        eigenvalue ranges associated with specific principal motions. This
        capability allows users to isolate particles corresponding to particular
        conformational regions along a motion trajectory.

        In practical cryo-EM analysis, this feature can be useful for separating
        extreme conformations, studying transition intermediates, or generating
        focused reconstructions representing distinct states along a continuous
        motion landscape.

        Outputs and Interpretation

        The protocol produces one refined volume for each defined body together
        with updated particle information and motion-related metadata. Each body
        reconstruction reflects the improved local alignment achieved through
        focused refinement.

        Additionally, the flexibility analysis generates principal component
        trajectories and eigenvector-associated map series suitable for
        visualization in molecular graphics software. These outputs provide both
        structural and dynamical insight into the conformational organization of
        the complex.

        Resolution estimates associated with the refined bodies help assess the
        quality of local refinement. In many cases, flexible regions that were
        poorly resolved in the consensus map become significantly clearer after
        multi-body analysis.

        Practical Recommendations

        Successful multi-body refinement depends strongly on biologically
        meaningful body definitions. It is generally advisable to begin with a
        small number of large, clearly identifiable bodies before attempting
        more detailed decompositions. Over-partitioning the complex often leads
        to unstable refinement and difficult interpretation.

        Smooth masks with limited overlap usually provide the most robust
        behavior. Visual inspection of the masks before refinement is highly
        recommended to ensure that flexible regions are represented properly
        without introducing disconnected fragments.

        The protocol is particularly effective for ribosomes, spliceosomes,
        membrane transporters, chaperones, viral assemblies, and other large
        complexes with coordinated domain movements. For relatively rigid
        particles with limited conformational variability, standard consensus
        refinement may remain sufficient.

        Final Perspective

        Multi-body refinement transforms cryo-EM analysis from a purely static
        reconstruction problem into a framework for studying structural dynamics.
        By combining focused refinement with quantitative motion analysis, the
        protocol enables biological interpretation of flexibility directly from
        experimental particle images.

        For many modern cryo-EM studies, understanding conformational landscapes
        is as important as achieving high nominal resolution. Careful body
        definition, thoughtful interpretation of principal motions, and
        validation against known biochemical behavior are essential for obtaining
        biologically meaningful insights from multi-body refinement.
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
        self._insertFunctionStep(self.convertInputStep, objId, needsGPU=False)
        self._insertFunctionStep(self.multibodyRefineStep,
                                 self._getRefineArgs(),
                                 needsGPU=self.usesGpu())
        if self.runFlexAnalysis:
            self._insertFunctionStep(self.flexAnalysisStep,
                                     self._getAnalyseArgs(),
                                     needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)
    
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
