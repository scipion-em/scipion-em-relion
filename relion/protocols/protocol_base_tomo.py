# *
# * Authors:     Scipion Team
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# *  e-mail address 'scipion-users@lists.sourceforge.net'
# *
# **************************************************************************
"""
This module contains the protocol base class for Relion protocols
"""

import re
from glob import glob
from os.path import exists, abspath

import pyworkflow.protocol.params as params
from pwem.objects import SetOfClasses3D
from pwem.protocols import EMProtocol
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import cleanPath

from pwem.convert.transformations import euler_from_matrix, translation_from_matrix

import relion
from relion.convert import convertBinaryVol, MASK_FILL_ZERO, Table, writeSetOfSubtomograms
from relion import ANGULAR_SAMPLING_LIST


class ProtRelionBaseTomo(EMProtocol):
    """ This class cointains the common functionalities for all Relion protocols.
    In subclasses there should be little changes about how to create the command line
    and the files produced.

    Most of the Relion protocols, have two modes: NORMAL or CONTINUE. That's why
    some of the function have a template pattern approach to define the behaivour
    depending on the case.
    """
    _label = '3d classify'
    OUTPUT_TYPE = SetOfClasses3D
    FILE_KEYS = ['data', 'optimiser', 'sampling']
    PREFIXES = ['']
    modelTable = Table()
    claasesTable = Table()
    dataTable = Table()
    nClasses = 0

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()

        self.ClassFnTemplate = '%(rootDir)s/relion_it%(iter)03d_class%(ref)03d.mrc:mrc'
        if not self.doContinue:
            self.continueRun.set(None)
        else:
            self.referenceVolume.set(None)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.extraIter = self._getExtraPath('relion_it%(iter)03d_')
        myDict = {
            'input_star': self._getPath('input_subtomograms.star'),
            'data_scipion': self.extraIter + 'data_scipion.sqlite',
            'projections': self.extraIter + '%(half)sclass%(ref3d)03d_projections.sqlite',
            'classes_scipion': self.extraIter + 'classes_scipion.sqlite',
            'model': self.extraIter + 'model.star',
            'shiny': self._getExtraPath('shiny.star'),
            'optimiser': self.extraIter + 'optimiser.star',
            'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
            'all_avgPmax_xmipp': self._getTmpPath('iterations_avgPmax_xmipp.xmd'),
            'all_changes_xmipp': self._getTmpPath('iterations_changes_xmipp.xmd'),
            'selected_volumes': self._getTmpPath('selected_volumes_xmipp.xmd'),
            'dataFinal': self._getExtraPath("relion_data.star"),
            'modelFinal': self._getExtraPath("relion_model.star"),
            'finalvolume': self._getExtraPath("relion_class%(ref3d)03d.mrc:mrc")

        }
        # add to keys, data.star, optimiser.star and sampling.star
        for key in self.FILE_KEYS:
            myDict[key] = self.extraIter + '%s.star' % key
            key_xmipp = key + '_xmipp'
            myDict[key_xmipp] = self.extraIter + '%s.xmd' % key
        # add other keys that depends on prefixes
        for p in self.PREFIXES:
            myDict['%smodel' % p] = self.extraIter + '%smodel.star' % p
            myDict['%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc:mrc'

        self._updateFilenamesDict(myDict)

    def _createIterTemplates(self):
        """ Setup the regex on how to find iterations. """
        self._iterTemplate = self._getFileName('data', iter=0).replace('000', '???')
        # Iterations will be identify by _itXXX_ where XXX is the iteration number
        # and is restricted to only 3 digits.
        self._iterRegex = re.compile('_it(\d{3,3})_')

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Some hidden variables to be used for conditions
        form.addHidden('isClassify', params.BooleanParam, default=self.IS_CLASSIFY)

        form.addParam('doContinue', params.BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                           'run of type *%s* class and most of the input parameters'
                           'will be taken from it.' % self.getClassName())
        form.addParam('inputSubtomograms', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      condition='not doContinue',
                      important=True,
                      label="Input subtomograms",
                      help='Select the input subtomograms from the project.')
        form.addParam('maskDiameterA', params.IntParam,
                      default=-1,
                      condition='not doContinue',
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a soft circular mask '
                           'with this <diameter>. '
                           'Make sure this diameter is not set too small because that may mask '
                           'away part of the signal! If set to a value larger than the image '
                           'size no masking will be performed.\n\n'
                           'The same diameter will also be used for a spherical mask of the '
                           'reference structures if no user-provided mask is specified.')
        form.addParam('continueRun', params.PointerParam,
                      pointerClass=self.getClassName(),
                      condition='doContinue',
                      allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('continueIter', params.StringParam, default='last',
                      condition='doContinue',
                      label='Continue from iteration',
                      help='Select from which iteration do you want to continue.'
                           'if you use *last*, then the last iteration will be used.'
                           'otherwise, a valid iteration number should be provided.')

        form.addParam('numberOfClasses', params.IntParam, default=3,
                      condition='not doContinue and isClassify',
                      label='Number of classes:',
                      help='The number of classes (K) for a multi-reference refinement.'
                           'These classes will be made in an unsupervised manner from a single'
                           'reference by division of the data into random subsets during the'
                           'first iteration.')
        group = form.addGroup('Reference 3D map',
                              condition='not doContinue')
        group.addParam('referenceVolume', params.PointerParam,
                       pointerClass='Volume',
                       allowsNull=True,
                       label="Input volume",
                       condition='not doContinue',
                       help='Initial reference 3D map, it should have the same '
                            'dimensions and the same pixel size as your input particles.')
        group.addParam('isMapAbsoluteGreyScale', params.BooleanParam, default=False,
                       label="Is initial 3D map on absolute greyscale?",
                       help='The probabilities are based on squared differences, '
                            'so that the absolute grey scale is important. \n'
                            'Probabilities are calculated based on a Gaussian noise model,'
                            'which contains a squared difference term between the reference and the experimental '
                            'image. This has a consequence that the reference needs to be on the same absolute '
                            'intensity greyscale as the experimental images. RELION and XMIPP reconstruct maps at '
                            'their absolute intensity greyscale. Other packages may perform internal normalisations of '
                            'the reference density, which will result in incorrect grey-scales. Therefore: if the map '
                            'was reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No. '
                            'If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the '
                            'first iteration, and prior to the second iteration the map will be filtered again using '
                            'the initial low-pass filter. This procedure is relatively quick and typically does not '
                            'negatively affect the outcome of the subsequent MAP refinement. Therefore, if in doubt it '
                            'is recommended to set this option to No.')

        self.addSymmetry(group)
        group.addParam('initialLowPassFilterA', params.FloatParam, default=60,
                       label='Initial low-pass filter (A)',
                       help='It is recommended to strongly low-pass filter your initial reference map. '
                            'If it has not yet been low-pass filtered, it may be done internally using this option. '
                            'If set to 0, no low-pass filter will be applied to the initial reference(s).')

        form.addSection(label='CTF')
        form.addParam('contuinueMsg', params.LabelParam, default=True,
                      label='CTF parameters are not available in continue mode', condition='doContinue', )
        form.addParam('doCTF', params.BooleanParam, default=True,
                      label='Do CTF-correction?', condition='not doContinue',
                      help='If set to Yes, CTFs will be corrected inside the MAP refinement. '
                           'The resulting algorithm intrinsically implements the optimal linear, '
                           'or Wiener filter. Note that input particles should contains CTF parameters.')
        form.addParam('ignoreCTFUntilFirstPeak', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore CTFs until first peak?', condition='doCTF and not doContinue',
                      help='If set to Yes, then CTF-amplitude correction will only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF model is inadequate at the lowest '
                           'resolution. Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) '
                           'often yields better results. Therefore, this option is not generally recommended.')
        form.addParam('ctfFlipPhases', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Perform only CTF phase-flipping?', condition='doCTF and not doContinue')
        form.addParam('hasReferenceCTFCorrected', params.BooleanParam, default=False,
                      condition='not doContinue',
                      label='Has reference been CTF-corrected?',
                      help='Set this option to Yes if the reference map represents CTF-unaffected density, '
                           'e.g. it was created using Wiener filtering inside RELION or from a PDB. If set to No, '
                           'then in the first iteration, the Fourier transforms of the reference projections '
                           'are not multiplied by the CTFs.')
        form.addParam('ctfMultiplied', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Has the data been pre-multiplied by the CTF?', condition='not doContinue')
        form.addParam('ctfPhaseFlipped', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Has the data been CTF phase-flipped?', condition='not doContinue')

        form.addSection(label='Optimisation')
        if self.IS_CLASSIFY:
            form.addParam('numberOfIterations', params.IntParam, default=25,
                          label='Number of iterations',
                          help='Number of iterations to be performed. Note that the current implementation does NOT '
                               'comprise a convergence criterium. Therefore, the calculations will need to be stopped '
                               'by the user if further iterations do not yield improvements in resolution or classes. '
                               'If continue option is True, you going to do this number of new iterations (e.g. if '
                               '*Continue from iteration* is set 3 and this param is set 25, the final iteration of the '
                               'protocol will be the 28th.')
            form.addParam('regularisationParamT', params.IntParam, default=2,
                          label='Regularisation parameter T',
                          help='Bayes law strictly determines the relative weight between the contribution of the '
                               'experimental data and the prior. '
                               'However, in practice one may need to adjust this weight to put slightly more weight on '
                               'the experimental '
                               'data to allow optimal results. Values greater than 1 for this regularisation parameter '
                               '(T in the JMB2011 paper) put more weight on the experimental data. Values around 2-4 '
                               'have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. '
                               'Too small values yield too-low resolution structures; too high values result in '
                               'over-estimated resolutions and overfitting.')
        form.addParam('maskZero', params.EnumParam, default=0,
                      choices=['Yes, fill with zeros', 'No, fill with random noise'],
                      label='Mask particles with zeros?', condition='not doContinue',
                      help='If set to <Yes>, then in the individual particles, the area outside a circle with the '
                           'radius of the particle will be set to zeros prior to taking the Fourier transform. '
                           'This will remove noise and therefore increase sensitivity in the alignment and '
                           'classification. However, it will also introduce correlations between the Fourier '
                           'components that are not modelled. When set to <No>, then the solvent area is filled with '
                           'random noise, which prevents introducing correlations.High-resolution refinements (e.g. '
                           'in 3D auto-refine) tend to work better when filling the solvent area with random noise, '
                           'some classifications go better when using zeros.')
        form.addParam('referenceMask', params.PointerParam, pointerClass='VolumeMask',
                      label='Reference mask (optional)', allowsNull=True,
                      help='A volume mask containing a (soft) mask with the same dimensions '
                           'as the reference(s), and values between 0 and 1, with 1 being 100% protein '
                           'and 0 being 100% solvent. The reconstructed reference map will be multiplied '
                           'by this mask. If no mask is given, a soft spherical mask based on the <radius> '
                           'of the mask for the experimental images will be applied.\n\n'
                           'In some cases, for example for non-empty icosahedral viruses, it is also useful '
                           'to use a second mask. Check _Advaced_ level and select another volume mask')
        form.addParam('solventMask', params.PointerParam, pointerClass='VolumeMask',
                      expertLevel=LEVEL_ADVANCED, allowsNull=True,
                      label='Second reference mask (optional)',
                      help='For all white (value 1) pixels in this second mask the '
                           'corresponding pixels in the reconstructed map are set to the average value of '
                           'these pixels. Thereby, for example, the higher density inside the virion may be '
                           'set to a constant. Note that this second mask should have one-values inside the '
                           'virion and zero-values in the capsid and the solvent areas.')

        if self.IS_CLASSIFY:
            form.addParam('limitResolEStep', params.FloatParam, default=-1,
                          label='Limit resolution E-step to (A)', condition="not doContinue",
                          help='If set to a positive number, then the expectation step '
                               '(i.e. the alignment) will be done only including the Fourier '
                               'components up to this resolution (in Angstroms). This is useful '
                               'to prevent overfitting, as the classification runs in RELION are '
                               'not to be guaranteed to be 100% overfitting-free (unlike the '
                               '_3D auto-refine_ with its gold-standard FSC). In particular for very '
                               'difficult data sets, e.g. of very small or featureless particles, '
                               'this has been shown to give much better class averages. In such '
                               'cases, values in the range of 7-12 Angstroms have proven useful.')
            # Change the Sampling section name depending if classify or refine 3D
            form.addSection('Sampling')
        else:
            form.addSection('Auto-Sampling')
            form.addParam('noteAutoSampling', params.LabelParam,
                          label='Note that initial sampling rates will be auto-incremented!')

        form.addParam('doImageAlignment', params.BooleanParam, default=False,
                      label='Perform image alignment?', condition="isClassify",
                      help='If set to No, then rather than performing both alignment '
                           'and classification, only classification will be performed. '
                           'This allows the use of very focused masks.This requires '
                           'that the optimal orientations of all particles are already '
                           'calculated.')
        form.addParam('angularSamplingDeg', params.EnumParam, default=2,
                      choices=ANGULAR_SAMPLING_LIST,
                      label='Angular sampling interval (deg)', condition='not isClassify or doImageAlignment',
                      help='There are only a few discrete angular samplings possible because '
                           'we use the HealPix library to generate the sampling of the first '
                           'two Euler angles on the sphere. The samplings are approximate numbers '
                           'and vary slightly over the sphere.')
        form.addParam('offsetSearchRangePix', params.FloatParam, default=5,
                      condition='not isClassify or doImageAlignment',
                      label='Offset search range (pix)',
                      help='Probabilities will be calculated only for translations in a circle '
                           'with this radius (in pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal translation for each '
                           'image in the previous iteration.')
        form.addParam('offsetSearchStepPix', params.FloatParam, default=1.0,
                      condition='not isClassify or doImageAlignment',
                      label='Offset search step (pix)',
                      help='Translations will be sampled with this step-size (in pixels). '
                           'Translational sampling is also done using the adaptive approach. '
                           'Therefore, if adaptive=1, the translations will first be evaluated'
                           'on a 2x coarser grid.')
        if self.IS_CLASSIFY:
            form.addParam('localAngularSearch', params.BooleanParam, default=False,
                          condition='doImageAlignment',
                          label='Perform local angular search?',
                          help='If set to Yes, then rather than performing exhaustive angular searches, '
                               'local searches within the range given below will be performed. A prior '
                               'Gaussian distribution centered at the optimal orientation in the previous '
                               'iteration and with a stddev of 1/3 of the range given below will be enforced. '
                               'NOTE: activate this option will enable the gpu calculating in Compute tab.')
            form.addParam('localAngularSearchRange', params.FloatParam, default=5.0,
                          condition='localAngularSearch',
                          label='Local angular search range',
                          help='Local angular searches will be performed within +/- '
                               'the given amount (in degrees) from the optimal orientation '
                               'in the previous iteration. A Gaussian prior (also see '
                               'previous option) will be applied, so that orientations '
                               'closer to the optimal orientation in the previous iteration '
                               'will get higher weights than those further away.')
        else:
            form.addParam('localSearchAutoSamplingDeg', params.EnumParam, default=4,
                          choices=ANGULAR_SAMPLING_LIST,
                          label='Local search from auto-sampling (deg)',
                          help='In the automated procedure to increase the angular samplings,\n'
                               'local angular searches of -6/+6 times the sampling rate will\n'
                               'be used from this angular sampling rate onwards.')

        form.addSection('Compute')
        self._defineComputeParams(form)

        form.addSection('Additional')
        joinHalves = ("--low_resol_join_halves 40 (only not continue mode)"
                      if not self.IS_CLASSIFY else "")

        form.addParam('keepOnlyLastIterFiles', params.BooleanParam,
                      default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Keep only files from last iteration?",
                      help="If Yes is chosen, only the files which correspond to the last "
                           "iteration will be saved in the protocol's extra directory. Otherwise, "
                           "files corresponding to each iteration will be kept.")
        form.addParam('oversampling', params.IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label="Over-sampling",
                      help="Adaptive oversampling order to speed-up "
                           "calculations (0=no oversampling, 1=2x, 2=4x, etc)")

        form.addParam('extraParams', params.StringParam,
                      default='',
                      label='Additional arguments',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options "
                           "and/or expert use of the program, e.g: \n"
                           "--verb 1\n"
                           "--pad 2\n" + joinHalves)

        form.addParallelSection(threads=1, mpi=3)

    def addSymmetry(self, container):
        container.addParam('symmetryGroup', params.StringParam, default='c1',
                           label="Symmetry",
                           help='If the molecule is asymmetric, set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n'
                                '* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196-207)\n'
                                '* I2: Crowther 222                                                                     \n'
                                '* I3: 52-setting (as used in SPIDER?)                                                  \n'
                                '* I4: A different 52 setting                                                           \n'
                                'The command *relion_refine --sym D2 --print_symmetry_ops* prints a list of all symmetry operators for symmetry group D2. RELION uses XMIPP\'s libraries for symmetry operations. Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry')

    @ staticmethod
    def _defineComputeParams(form):
        form.addParam('useParallelDisk', params.BooleanParam, default=True,
                      label='Use parallel disc I/O?',
                      help='If set to Yes, all MPI slaves will read '
                           'their own images from disc. Otherwise, only '
                           'the master will read images and send them '
                           'through the network to the slaves. Parallel '
                           'file systems like gluster of fhgfs are good '
                           'at parallel disc I/O. NFS may break with many '
                           'slaves reading in parallel.')
        form.addParam('pooledSubtomos', params.IntParam, default=3,
                      label='Number of pooled particles:',
                      help='Particles are processed in individual batches '
                           'by MPI slaves. During each batch, a stack of '
                           'particle images is only opened and closed '
                           'once to improve disk access times. All '
                           'particle images of a single batch are read '
                           'into memory together. The size of these '
                           'batches is at least one particle per thread '
                           'used. The nr_pooled_particles parameter '
                           'controls how many particles are read together '
                           'for each thread. If it is set to 3 and one '
                           'uses 8 threads, batches of 3x8=24 particles '
                           'will be read together. This may improve '
                           'performance on systems where disk access, and '
                           'particularly metadata handling of disk '
                           'access, is a problem. It has a modest cost of '
                           'increased RAM usage.')
        form.addParam('allParticlesRam', params.BooleanParam, default=False,
                      label='Pre-read all particles into RAM?',
                      help='If set to Yes, all particle images will be '
                           'read into computer memory, which will greatly '
                           'speed up calculations on systems with slow '
                           'disk access. However, one should of course be '
                           'careful with the amount of RAM available. '
                           'Because particles are read in '
                           'float-precision, it will take \n'
                           '( N * (box_size)^2 * 4 / (1024 * 1024 '
                           '* 1024) ) Giga-bytes to read N particles into '
                           'RAM. For 100 thousand 200x200 images, that '
                           'becomes 15Gb, or 60 Gb for the same number of '
                           '400x400 particles. Remember that running a '
                           'single MPI slave on each node that runs as '
                           'many threads as available cores will have '
                           'access to all available RAM.\n\n'
                           'If parallel disc I/O is set to No, then only '
                           'the master reads all particles into RAM and '
                           'sends those particles through the network to '
                           'the MPI slaves during the refinement '
                           'iterations.')
        form.addParam('combineItersDisc', params.BooleanParam, default=False,
                      label='Combine iterations through disc?',
                      help='If set to Yes, at the end of every iteration '
                           'all MPI slaves will write out a large file '
                           'with their accumulated results. The MPI '
                           'master will read in all these files, combine '
                           'them all, and write out a new file with the '
                           'combined results. All MPI slaves will then '
                           'read in the combined results. This reduces '
                           'heavy load on the network, but increases load '
                           'on the disc I/O. This will affect the time it '
                           'takes between the progress-bar in the '
                           'expectation step reaching its end (the mouse '
                           'gets to the cheese) and the start of the '
                           'ensuing maximisation step. It will depend on '
                           'your system setup which is most efficient.')
        form.addParam('doGpu', params.BooleanParam, default=False,
                      condition='doImageAlignment or not isClassify',
                      label='Use GPU acceleration?',
                      help='If set to Yes, the job will try to use GPU '
                           'acceleration. This will be ignored if skipping alignment '
                           '(From relion_refine_mpi, classification case: '
                           'you cannot use accelerators when '
                           'skipping alignments).')
        form.addParam('gpusToUse', params.StringParam,
                      default='', condition='doGpu and (doImageAlignment or not isClassify)',
                      label='GPUs to use:',
                      help='It can be used to provide a list of which GPUs '
                           '(e. g. "0:1:2:3") to use. MPI-processes are '
                           'separated by ":", threads by ",". '
                           'For example: "0,0:1,1:0,0:1,1"')
        form.addParam('scratchDir', params.PathParam,
                      condition='not allParticlesRam',
                      label='Copy particles to scratch directory (optional): ',
                      help='If a directory is provided here, then the job '
                           'will create a sub-directory in it called '
                           'relion_volatile. If that relion_volatile '
                           'directory already exists, it will be wiped. '
                           'Then, the program will copy all input '
                           'particles into a large stack inside the '
                           'relion_volatile subdirectory. Provided this '
                           'directory is on a fast local drive (e.g. an '
                           'SSD drive), processing in all the iterations '
                           'will be faster. If the job finishes '
                           'correctly, the relion_volatile directory will '
                           'be wiped. If the job crashes, you may want to '
                           'remove it yourself.')

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('convertInputStep')
        self._insertRelionStep()
        self._insertFunctionStep('createOutputStep')

    def _insertRelionStep(self):
        """ Prepare the command line arguments before calling Relion. """

        # Join in a single line all key, value pairs of the args dict
        args = {}

        if self.doContinue:
            self._setContinueArgs(args)
        else:
            self._setNormalArgs(args)
        params = ' '.join(['%s %s' % (k, str(v)) for k, v in args.items()])

        if self.extraParams.hasValue():
            params += ' ' + self.extraParams.get()

        self._insertFunctionStep('runRelionStep', params)

    # --------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion."""

        subtomoSet = self._getInputParticles()
        subtomosStar = self._getFileName('input_star')
        self.info("Converting set from '%s' into '%s'" % (subtomoSet.getFileName(), subtomosStar))
        writeSetOfSubtomograms(subtomoSet, subtomosStar)

    def runRelionStep(self, params):
        """ Execute the relion steps with the give params. """

        params += ' --j %d' % self.numberOfThreads.get()
        self.runJob(self._getProgram(), params)

    def createOutputStep(self):
        pass  # should be implemented in subclasses

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        if self.doContinue:
            continueProtocol = self.continueRun.get()
            if (continueProtocol is not None and
                    continueProtocol.getObjId() == self.getObjId()):
                errors.append('In Scipion you must create a new Relion run')
                errors.append('and select the continue option rather than')
                errors.append('select continue from the same run.')
                errors.append('')  # add a new line
            errors += self._validateContinue()
        else:
            if self._getInputParticles().isOddX():
                errors.append("Relion only works with even values for the image dimensions!")

        if self.ctfPhaseFlipped and self.ctfMultiplied:
            errors.append("CTF options\n"
                          "*Has the data been pre-multiplied by the CTF?* \n"
                          "and \n"
                          "*Has the data been CTF phase-flipped?* \n"
                          "are *exclusive*. Only one of them can be enabled at the same time.")

            errors += self._validateNormal()
        return errors

    def _validateNormal(self):
        """ Should be overriden in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []

    def _validateContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []

    def _citations(self):
        cites = []
        return cites

    def _summary(self):
        self._initialize()
        iter = self._lastIter()
        if iter >= 0:
            iterMsg = 'Calculating iteration %d' % iter
            if self.hasAttribute('numberOfIterations'):
                iterMsg += '/%d' % self._getnumberOfIters()
            summary = [iterMsg]
        else:
            summary = ['Initializing...']
        if self.doContinue:
            summary += self._summaryContinue()
        summary += self._summaryNormal()
        return summary

    def _summaryNormal(self):
        """ Should be overriden in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        return []

    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []

    def _methods(self):
        """ Should be overriden in each protocol.
        """
        return []

    # --------------------------- UTILS functions --------------------------------------------
    def _getProgram(self, program='relion_refine'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def _getInputParticles(self):
        if self.doContinue:
            self.inputSubtomograms.set(self.continueRun.get().inputSubtomograms.get())
        return self.inputSubtomograms.get()

    def _getIterNumber(self, index):
        """ Return the list of iteration files, give the iterTemplate. """
        result = -1
        files = sorted(glob(self._iterTemplate))
        if files:
            f = files[index]
            s = self._iterRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is 3 digits iteration number
        return result

    def _lastIter(self):
        return self._getIterNumber(-1)

    def _firstIter(self):
        return self._getIterNumber(0) or 1

    def _getIterClasses(self, it, clean=False):
        """ Return a classes .sqlite file for this iteration.
        If the file doesn't exists, it will be created by
        converting from this iteration data.star file.
        """
        data_classes = self._getFileName('classes_scipion', iter=it)

        if clean:
            cleanPath(data_classes)

        if not exists(data_classes):
            clsSet = self.OUTPUT_TYPE(filename=data_classes)
            clsSet.setImages(self.inputParticles.get())
            self._fillClassesFromIter(clsSet, it)
            clsSet.write()
            clsSet.close()

        return data_classes

    def _getContinueIter(self):
        continueRun = self.continueRun.get()

        if continueRun is not None:
            continueRun._initialize()

        if self.doContinue:
            if self.continueIter.get() == 'last':
                continueIter = continueRun._lastIter()
            else:
                continueIter = int(self.continueIter.get())
        else:
            continueIter = 0

        return continueIter

    def _getnumberOfIters(self):
        return self._getContinueIter() + self.numberOfIterations.get()

    def _setNormalArgs(self, args):
        # Params only appear in the command line if their value is different than the default
        maskDiameter = self.maskDiameterA.get()
        if maskDiameter <= 0:
            x, _, _ = self._getInputParticles().getDim()
            maskDiameter = self._getInputParticles().getSamplingRate() * x

        args.update({'--i': abspath(self._getFileName('input_star')),
                     '--particle_diameter': maskDiameter,
                     '--angpix': self._getInputParticles().getSamplingRate(),
                     })
        self._setMaskArgs(args)
        self._setCTFArgs(args)

        if self.maskZero.get() == MASK_FILL_ZERO:
            # --zero_mask (false) : Mask surrounding background in particles to zero
            # (by default the solvent area is filled with random noise)
            args['--zero_mask'] = ''

        if self.IS_CLASSIFY:
            # --K (1) : Number of references to be refined
            args['--K'] = self.numberOfClasses.get()
            if self.limitResolEStep.get() > 0:
                # --strict_highres_exp (-1) : Resolution limit (in Angstrom) to
                # restrict probability calculations in the expectation step
                args['--strict_highres_exp'] = self.limitResolEStep.get()
        refVol = self.referenceVolume.get()
        if refVol is not None:
            args['--ref'] = abspath(convertBinaryVol(refVol, self._getTmpPath()))
        if self.isMapAbsoluteGreyScale.get():
            # --firstiter_cc (false) : Perform CC-calculation in the first iteration (use
            # this if references are not on the absolute intensity scale)
            args['--firstiter_cc'] = ''
        # --ini_high (-1) : Resolution (in Angstroms) to which to limit refinement in
        # the first iteration
        args['--ini_high'] = self.initialLowPassFilterA.get()
        # Symmetry group
        args['--sym'] = self.symmetryGroup.get()

        # args['--memory_per_thread'] = self.memoryPreThreads.get()
        self._setBasicArgs(args)
        self._setComputeArgs(args)

    def _setContinueArgs(self, args):
        continueRun = self.continueRun.get()
        continueRun._initialize()

        if self.IS_CLASSIFY:
            self.copyAttributes(continueRun, 'regularisationParamT')
        self._setBasicArgs(args)

        args['--continue'] = continueRun._getFileName('optimiser', iter=self._getContinueIter())

    def _setComputeArgs(self, args):
        # Params only appear in the command line if their value is different than the default
        if self.useParallelDisk.get():
            # --no_parallel_disc_io (false) : Do NOT let parallel (MPI) processes access the
            # disc simultaneously (use this option with NFS)
            args['--no_parallel_disc_io'] = ''
        # --pool (1) : Number of images to pool for each thread task
        args['--pool'] = self.pooledSubtomos.get()
        if self.allParticlesRam.get():
            # --preread_images (false) : Use this to let the master process read all particles
            # into memory. Be careful you have enough RAM for large data sets!
            args['--preread_images'] = ''
        scracthDir = self.scratchDir.get()
        if scracthDir:
            # --scratch_dir () : If provided, particle stacks will be copied to this local
            # scratch disk prior to refinement.
            args['--scratch_dir'] = ''
        if not self.combineItersDisc.get():
            # --dont_combine_weights_via_disc (false) : Send the large arrays of summed
            # weights through the MPI network, instead of writing large files to disc
            args['--dont_combine_weights_via_disc'] = ''
        if self.doGpu.get():
            if (self.IS_CLASSIFY and self.doImageAlignment.get()) or not self.IS_CLASSIFY:
                # Captured from relion_refine_mpi:
                # ERROR: you cannot use accelerators when skipping alignments
                args['--gpu'] = self.gpusToUse.get()

    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args.update({'--flatten_solvent': '',
                     '--norm': '',
                     '--scale': '',
                     '--o': abspath(self._getExtraPath('relion')),
                     '--oversampling': self.oversampling.get()
                     })

        if self.IS_CLASSIFY:
            args['--tau2_fudge'] = self.regularisationParamT.get()
            args['--iter'] = self._getnumberOfIters()

        self._setSamplingArgs(args)

    def _setCTFArgs(self, args):
        # Params only appear in the command line if their value is different than the default
        if self.doCTF.get():
            # --ctf (false) : Perform CTF correction?
            args['--ctf'] = ''

        if self.hasReferenceCTFCorrected.get():
            # --ctf_corrected_ref (false) : Have the input references been CTF-amplitude corrected?
            args['--ctf_corrected_ref'] = ''

        if self.ignoreCTFUntilFirstPeak.get():
            # --ctf_intact_first_peak (false) : Ignore CTFs until their first peak?
            args['--ctf_intact_first_peak'] = ''

        if self.ctfMultiplied.get():
            # --ctf_multiplied (false) : Have the data been premultiplied with their CTF?
            args['--ctf_multiplied'] = ''
        if self.ctfPhaseFlipped.get():
            # --ctf_phase_flipped (false) : Have the data been CTF phase-flipped?
            args['--ctf_multiplied'] = ''
        if self.ctfFlipPhases.get():
            # --only_flip_phases (false) : Only perform CTF phase-flipping?
            # (default is full amplitude-correction)
            args['--only_flip_phases'] = ''

    def _setMaskArgs(self, args):
        if self.referenceMask.hasValue():
            args[
                '--solvent_mask'] = abspath(self.referenceMask.get().getFileName())  # FIXE: CHANGE BY LOCATION, convert if necessary

        if self.solventMask.hasValue():
            args[
                '--solvent_mask2'] = abspath(self.solventMask.get().getFileName())  # FIXME: CHANGE BY LOCATION, convert if necessary

    def _cleanUndesiredFiles(self):
        """Remove all files generated by relion_refine excepting the ones which
        correspond to the last iteration."""
        pass

    def _getRlnCurrentResolution(self, starFile):
        with open(starFile) as fid:
            self.modelTable.readStar(fid, 'model_general')
            # Model table has only one row, while classes table has the same number of rows as classes found
            return float(self.modelTable._rows[0].rlnCurrentResolution)

    @ staticmethod
    def _getCTFFileFromSubtomo(subtomo):
        return subtomo.getCoordinate3D()._3dcftMrcFile.get()

    def IS_GT30(self):
        return relion.Plugin.IS_GT30()
