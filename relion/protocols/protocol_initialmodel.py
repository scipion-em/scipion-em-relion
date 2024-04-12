# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
# * [2] SciLifeLab, Stockholm University
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
from emtable import Table

from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtInitialVolume
from pwem.objects import Volume, SetOfVolumes, SetOfParticles
from pyworkflow.constants import PROD
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, IntParam,
                                        StringParam, BooleanParam)

import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionInitialModel(ProtInitialVolume, ProtRelionBase):
    """ This protocols creates a 3D initial model using Relion.

    Generate a 3D initial model _de novo_ from 2D particles using
    Relion Stochastic Gradient Descent (SGD) algorithm.
    """
    _label = '3D initial model'
    _devStatus = PROD
    _possibleOutputs = {
        "outputVolume": Volume,
        "outputVolumeSymmetrized": Volume,
        "outputVolumes": SetOfVolumes,
        "outputParticles": SetOfParticles
    }
    IS_CLASSIFY = False
    IS_3D_INIT = True
    IS_2D = False
    CHANGE_LABELS = ['rlnChangesOptimalOrientations',
                     'rlnChangesOptimalOffsets']

    def __init__(self, **args):
        ProtRelionBase.__init__(self, **args)

    def _initialize(self):
        """ This function is meant to be called after the
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        self._createFilenameTemplates()
        self._createIterTemplates()
        if not self.doContinue:
            self.continueRun.set(None)
        self.maskZero = False
        self.copyAlignment = False
        self.hasReferenceCTFCorrected = False

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineConstants()

        self.IS_3D = not self.IS_2D
        form.addSection(label='Input')
        form.addParam('doContinue', BooleanParam, default=False,
                      label='Continue from a previous run?',
                      help='If you set to *Yes*, you should select a previous'
                           'run of type *%s* class and most of the input '
                           'parameters will be taken from it.'
                           % self.getClassName())
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not doContinue',
                      important=True,
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('maskDiameterA', IntParam, default=-1,
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a '
                           'soft circular mask with this <diameter>. '
                           'Make sure this diameter is not set too small '
                           'because that may mask away part of the signal! If '
                           'set to a value larger than the image size no '
                           'masking will be performed.\n\n'
                           'The same diameter will also be used for a '
                           'spherical mask of the reference structures if no '
                           'user-provided mask is specified.')
        form.addParam('continueRun', PointerParam,
                      pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('continueIter', StringParam, default='last',
                      condition='doContinue',
                      label='Continue from iteration',
                      help='Select from which iteration do you want to '
                           'continue. If you use *last*, then the last '
                           'iteration will be used. Otherwise, a valid '
                           'iteration number should be provided.')

        self.addSymmetry(form)

        form.addSection(label='CTF')
        form.addParam('continueMsg', LabelParam, default=True,
                      condition='doContinue',
                      label='CTF parameters are not available in continue mode')
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-correction?', condition='not doContinue',
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('haveDataBeenPhaseFlipped', LabelParam,
                      condition='not doContinue',
                      label='Have data been phase-flipped?      '
                            '(Don\'t answer, see help)',
                      help='The phase-flip status is recorded and managed by '
                           'Scipion. \n In other words, when you import or '
                           'extract particles, \nScipion will record whether '
                           'or not phase flipping has been done.\n\n'
                           'Note that CTF-phase flipping is NOT a necessary '
                           'pre-processing step \nfor MAP-refinement in '
                           'RELION, as this can be done inside the internal\n'
                           'CTF-correction. However, if the phases have been '
                           'flipped, the program will handle it.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      label='Ignore CTFs until first peak?',
                      condition='not doContinue',
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')
        form.addParam('doCtfManualGroups', BooleanParam, default=False,
                      label='Do manual grouping ctfs?',
                      condition='not doContinue',
                      help='Set this to Yes the CTFs will grouping manually.')
        form.addParam('defocusRange', FloatParam, default=1000,
                      label='Defocus range for group creation (in Angstroms)',
                      condition='doCtfManualGroups and not doContinue',
                      help='Particles will be grouped by defocus.'
                           'This parameter is the bin for a histogram.'
                           'All particles assigned to a bin form a group')
        form.addParam('numParticles', FloatParam, default=10,
                      label='minimum size for defocus group',
                      condition='doCtfManualGroups and not doContinue',
                      help='If defocus group is smaller than this value, '
                           'it will be expanded until number of particles '
                           'per defocus group is reached')

        form.addSection('Optimisation')
        form.addParam('numberOfIter', IntParam, default=200,
                      label="Number of VDAM mini-batches",
                      help="How many iterations (i.e. mini-batches) "
                           "to perform with the VDAM algorithm?")
        form.addParam('regularisationParamT', FloatParam,
                      default=4,
                      label='Regularisation parameter T',
                      help='Bayes law strictly determines the relative '
                           'weight between the contribution of the '
                           'experimental data and the prior. '
                           'However, in practice one may need to adjust '
                           'this weight to put slightly more weight on the '
                           'experimental data to allow optimal results. '
                           'Values greater than 1 for this regularisation '
                           'parameter (T in the JMB2011 paper) put more '
                           'weight on the experimental data. Values around '
                           '2-4 have been observed to be useful for 3D '
                           'initial model calculations.')
        form.addParam('numberOfClasses', IntParam, default=1,
                      condition='not doContinue',
                      label='Number of classes',
                      help='The number of classes (K) for a multi-reference '
                           'ab initio SGD refinement. These classes will be '
                           'made in an unsupervised manner, starting from a '
                           'single reference in the initial iterations of '
                           'the SGD, and the references will become '
                           'increasingly dissimilar during the in between '
                           'iterations.')

        form.addParam('doFlattenSolvent', BooleanParam, default=True,
                      condition='not doContinue',
                      label='Flatten and enforce non-negative solvent?',
                      help='If set to Yes, the job will apply a spherical '
                           'mask and enforce all values in the reference '
                           'to be non-negative.')

        form.addParam('symmetryGroup', StringParam, default='c1',
                      condition='not doContinue',
                      label="Symmetry",
                      help='The initial model is always generated in C1 '
                           'and then aligned to and symmetrized with the '
                           'specified point group. If the automatic alignment '
                           'fails, please manually rotate run_itNNN_class001.mrc '
                           '(NNN is the number of iterations) so that it '
                           'conforms the symmetry convention.')

        form.addParam('runInC1', BooleanParam, default=True,
                      condition='not doContinue',
                      label='Run in C1 and apply symmetry later?',
                      help='If set to Yes, the gradient-driven optimisation '
                           'is run in C1 and the symmetry orientation is searched '
                           'and applied later. If set to No, the entire '
                           'optimisation is run in the symmetry point group '
                           'indicated above.')

        form.addSection('Compute')
        self._defineComputeParams(form)

        form.addParam('extraParams', StringParam, default='',
                      label='Additional arguments',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options "
                           "and/or expert use of the program, e.g: \n"
                           "--dont_combine_weights_via_disc\n"
                           "--verb 1\n"
                           "--pad 2\n")

        form.addParallelSection(threads=1, mpi=0)

    def addSymmetry(self, container):
        pass

    # -------------------------- INSERT steps functions -----------------------

    # -------------------------- STEPS functions ------------------------------
    def _getVolumes(self):
        """ Return the list of volumes generated.
        The number of volumes in the list will be equal to
        the number of classes requested by the user in the protocol. """
        k = self.numberOfClasses.get()
        pixelSize = self._getInputParticles().getSamplingRate()
        lastIter = self._lastIter()
        volumes = []

        for i in range(1, k + 1):
            vol = Volume(self._getExtraPath('relion_it%03d_class%03d.mrc')
                         % (lastIter, i))
            vol.setSamplingRate(pixelSize)
            volumes.append(vol)

        return volumes

    def createOutputStep(self):
        imgSet = self._getInputParticles()
        volumes = self._getVolumes()

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet, self._lastIter())

        if len(volumes) > 1:
            output = self._createSetOfVolumes()
            output.setSamplingRate(imgSet.getSamplingRate())
            for vol in volumes:
                output.append(vol)
            self._defineOutputs(outputVolumes=output)
        else:
            output = volumes[0]
            self._defineOutputs(outputVolume=output)

        self._defineSourceRelation(self.inputParticles, output)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

        # Run symmetrization for the largest volume if necessary
        sym = self.symmetryGroup.get()
        if sym.lower() != "c1" and self.runInC1:
            inFn = self._getFileName('model', iter=self._lastIter())
            symFn = self._getPath('volume_sym%s.mrc' % sym)
            pixSize = imgSet.getSamplingRate()

            self.runJob("relion_align_symmetry",
                        "--apply_sym --select_largest_class "
                        "--i %s --o %s --sym %s --angpix %0.5f" % (
                            inFn, symFn, sym, pixSize))

            vol = Volume()
            vol.copyInfo(volumes[0])
            vol.setLocation(symFn)
            self._defineOutputs(outputVolumeSymmetrized=vol)

    # -------------------------- INFO functions -------------------------------
    def _validateNormal(self):
        errors = []

        return errors

    def _validateContinue(self):
        errors = []
        continueRun = self.continueRun.get()
        continueRun._initialize()
        lastIter = continueRun._lastIter()

        if self.continueIter.get() == 'last':
            continueIter = lastIter
        else:
            continueIter = int(self.continueIter.get())

        if continueIter > lastIter:
            errors += ["The iteration from you want to continue must be "
                       "%01d or less" % lastIter]

        return errors

    def _summaryNormal(self):
        summary = []
        it = self._lastIter() or -1
        if it >= 1:
            table = Table(fileName=self._getFileName('model', iter=it),
                          tableName='model_general')
            row = table[0]
            resol = float(row.rlnCurrentResolution)
            summary.append("Current resolution: *%0.2f*" % resol)
        return summary

    def _summaryContinue(self):
        summary = ["Continue from iteration %01d" % self._getContinueIter()]
        return summary

    # -------------------------- UTILS functions ------------------------------
    def _setBasicArgs(self, args):
        """ Return a dictionary with basic arguments. """
        args.update({'--o': self._getExtraPath('relion'),
                     '--oversampling': 1,
                     '--pad': 1,
                     '--tau2_fudge': self.regularisationParamT.get()
                     })

        if self.doFlattenSolvent:
            args['--flatten_solvent'] = ''
        if not self.doContinue:
            args['--sym'] = 'C1' if self.runInC1 else self.symmetryGroup.get()

        self._setGradArgs(args)
        self._setSamplingArgs(args)

    def _setGradArgs(self, args):
        args['--iter'] = self.numberOfIter.get()
        args['--grad'] = ''
        args['--K'] = self.numberOfClasses.get()

        if not self.doContinue:
            args['--denovo_3dref'] = ''

    def _setSamplingArgs(self, args):
        """ Set sampling related params"""
        if not self.doContinue:
            args['--healpix_order'] = 1
            args['--offset_range'] = 6
            args['--offset_step'] = 2
            args['--auto_sampling'] = ''

    def _fillDataFromIter(self, imgSet, iteration):
        outImgsFn = self._getFileName('data', iter=iteration)
        imgSet.setAlignmentProj()
        px = imgSet.getSamplingRate()
        self.reader = convert.createReader(alignType=ALIGN_PROJ,
                                           pixelSize=px)
        mdIter = convert.Table.iterRows('particles@' + outImgsFn, key='rlnImageId',
                                        types=convert.LABELS_DICT)
        imgSet.copyItems(self._getInputParticles(), doClone=False,
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=mdIter)

    def _createItemMatrix(self, item, row):
        self.reader.setParticleTransform(item, row)
