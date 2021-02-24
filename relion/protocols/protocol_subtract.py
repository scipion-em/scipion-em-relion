# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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

from pyworkflow.protocol.constants import LEVEL_ADVANCED
# import pwem.emlib.metadata as md
from pyworkflow.object import String, Integer
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, LabelParam
from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtOperateParticles

import relion.convert as convert
from relion import Plugin

import os
import sys

class ProtRelionSubtract(ProtOperateParticles):
    """ Signal subtraction protocol of Relion.

    Subtract volume projections from the experimental particles.
    The particles must have projection alignment in order to
    properly generate volume projections.
    """
    _label = 'subtract projection'

    def _initialize(self):
        self._createFilenameTemplates()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  'input_star': self._getExtraPath('input_particles.star'),
                  'output_star': self._getExtraPath('particles_subtracted.star')
                  }
        #if not self.isRelionInput:
        #    myDict['output'] = self._getExtraPath('particles_subtracted')
        self._updateFilenamesDict(myDict)
    
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('relionInput', BooleanParam,
                      default = True,
                      important=True,
                      label="Start from Relion 3D refinement?",
                      help="Set to yes if you wish to use as input"
                           " a Relion 3D refinement protocol "
                           "Otherwise set it to No")
        # input volume
        form.addParam('inputProtocol', PointerParam,
                      important=True,
                      pointerClass='ProtRelionRefine3D, ProtRelionClassify3D',
                      label="Input 3D protocol",
                      condition="relionInput",
                      help="Select the 3D refinement/classification run which "
                           "you want to use for subtraction. It will use the "
                           "maps from this run for the subtraction.")
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input map to be projected",
                      important=True,
                      condition="not relionInput",
                      help='Provide the input volume that will be used to '
                           'calculate projections, which will be subtracted '
                           'from the experimental particles. Make sure this '
                           'map was calculated by RELION from the same '
                           'particles as above, and preferably with those '
                           'orientations, as it is crucial that the absolute '
                           'greyscale is the same as in the experimental '
                           'particles.')
        # input particles
        form.addParam('useAll', BooleanParam, default=True,
                      label="Use all particles from input protocol?",
                      condition="relionInput",
                      help="If No, then you need to provide a subset of "
                           "particles below.")

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not useAll',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles subset",
                      help='Select the particles which are a SUBSET of the '
                           'input protocol provided above.')

        form.addParam('inputParticlesAll', PointerParam,
                      pointerClass='SetOfParticles',
                      condition="not relionInput",
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the experimental particles.')

        #reference mask
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask of the signal to keep',
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is black (0) and the density you "
                           "wish to keep is white (1).\n"
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to KEEP.*")

        form.addParam('invertMask', BooleanParam, default=False,
                      label='Invert Mask',
                      condition="not relionInput",
                      help="Invert the provided mask so 0 became 1 and viceversa"
                      )

        form.addSection('Centering')
        form.addParam('help1', LabelParam,
                      label="Empty section (see help) ",
                      condition="not relionInput",
                      help="This section is only used if "
                            "'Start from Relion 3D refinement' is set to True")
        form.addParam('centerOnMask', BooleanParam, default=True,
                      label="Do center subtracted images on mask?",
                      condition="relionInput",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the "
                           "center-of-mass of the input mask.")
        form.addParam('centerOnCoord', BooleanParam, default=False,
                      condition='(not centerOnMask) and relionInput',
                      label="Do center on my coordinates?",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the x,y,z "
                           "coordinates below. The unit is pixel, not "
                           "angstrom. The origin is at the center of the box, "
                           "not at the corner.")

        line = form.addLine('Center coordinate (px)',
                            condition='centerOnCoord and relionInput',
                            help='Coordinate of the 3D center (in pixels).')
        line.addParam('cX', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='X')
        line.addParam('cY', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='Y')
        line.addParam('cZ', IntParam, default=0, condition='centerOnCoord and relionInput',
                      label='Z')

        form.addParam('newBoxSize', IntParam, default=-1,
                      condition="relionInput",
                      label="New box size",
                      help="Provide a non-negative value to re-window the "
                           "subtracted particles in a smaller box size.")

        form.addSection(label='CTF')
        form.addParam('help', LabelParam,
                      label="Empty section (see help) ",
                      condition="relionInput",
                      help="This section is only used if "
                            "'Start from Relion 3D refinement' is set to False")
        form.addParam('doCTF', BooleanParam, default=True,
                      label='Do CTF-correction?',
                      condition="not relionInput",
                      help='If set to Yes, CTFs will be corrected inside the '
                           'MAP refinement. The resulting algorithm '
                           'intrinsically implements the optimal linear, or '
                           'Wiener filter. Note that input particles should '
                           'contains CTF parameters.')
        form.addParam('haveDataBeenPhaseFlipped', LabelParam,
                      label='Have data been phase-flipped?      '
                            '(Don\'t answer, see help)',
                      condition="not relionInput",
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
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore CTFs until first peak?',
                      condition="not relionInput",
                      help='If set to Yes, then CTF-amplitude correction will '
                           'only be performed from the first peak '
                           'of each CTF onward. This can be useful if the CTF '
                           'model is inadequate at the lowest resolution. '
                           'Still, in general using higher amplitude contrast '
                           'on the CTFs (e.g. 10-20%) often yields better '
                           'results. Therefore, this option is not generally '
                           'recommended.')

        form.addParallelSection(threads=0, mpi=1)
    
    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self.isRelionInput = self.relionInput.get()
        self._initialize()

        if not self.useAll or not self.isRelionInput:
            self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('subtractStep')
        self._insertFunctionStep('createOutputStep')
    
    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ Write the input images as a Relion star file. """
        if not self.isRelionInput:
            self.inputParticles = self.inputParticlesAll
        imgSet = self.inputParticles.get()

        convert.writeSetOfParticles(
            imgSet, self._getFileName('input_star'),
            outputDir=self._getExtraPath(), alignType=ALIGN_PROJ)

    def subtractStep(self):
        if self.isRelionInput:
            self.subtractStepRelion()
        else:
            self.subtractStepNoRelion()

    def subtractStepNoRelion(self):
        self.numberOfMpi = Integer(1)
        volume = self.inputVolume.get()
        volFn = convert.convertBinaryVol(volume,
                                         self._getExtraPath())
        params = ' --i %s --subtract_exp --angpix %0.3f' % (volFn,
                                                            volume.getSamplingRate())

        if self.refMask.hasValue():
            tmp = self._getTmpPath()
            newDim = self._getInputParticles().getXDim()
            newPix = self._getInputParticles().getSamplingRate()

            maskFn = convert.convertMaskThreshold(self.refMask.get(),
                                                  tmp,
                                                  newPix,
                                                  newDim,
                                                  self.invertMask)
            params += ' --mask %s' % maskFn

        if self._getInputParticles().isPhaseFlipped():
            params += ' --ctf_phase_flip'

        if self.doCTF:
            params += ' --ctf'
            if self.ignoreCTFUntilFirstPeak:
                params += ' --ctf_intact_first_peak'

        params += ' --ang %s  --o %s ' % (
            self._getFileName('input_star'),
            self._getFileName('output_star').replace(".star", ""))

        try:
            self.runJob('relion_project', params)
        except Exception as ex:
            fn = self._getFileName('output_star')
            if not os.path.exists(fn):
                sys.stderr.write('The file %s was not produced\n' % fn)
                raise ex
            else:
                sys.stderr.write('----Everything OK-----\n')

    def subtractStepRelion(self):
        inputProt = self.inputProtocol.get()
        inputProt._initialize()
        fnOptimiser = inputProt._getFileName('optimiser',
                                             iter=inputProt._lastIter())
        params = " --i %s --o %s --new_box %s" % (fnOptimiser,
                                                  self._getExtraPath(),
                                                  self.newBoxSize.get())

        if not self.useAll:
            params += " --data %s" % self._getFileName('input_star')
        if self.centerOnMask:
            params += " --recenter_on_mask"
        elif self.centerOnCoord:
            params += " --center_x %d --center_y %d --center_z %d" % (
                self.cX, self.cY, self.cZ)

        tmp = self._getTmpPath()
        newDim = self._getInputParticles().getXDim()
        newPix = self._getInputParticles().getSamplingRate()
        maskFn = convert.convertMask(self.refMask.get(), tmp, newPix, newDim)
        params += ' --mask %s' % maskFn

        prog = "relion_particle_subtract" + ("_mpi" if self.numberOfMpi > 1 else "")
        self.runJob(prog, params)

    def createOutputStep(self):
        imgSet = self._getInputParticles()
        outImgSet = self._createSetOfParticles()
        outImgsFn = self._getFileName('output_star')
        outImgSet.copyInfo(imgSet)
        outImgSet.setAlignmentProj()

        self.reader = convert.createReader(alignType=ALIGN_PROJ)
        tableName = 'particles@' if self.IS_GT30() else ''
        if self.isRelionInput:
            doClone = False
        else:
            doClone = True
            self.reader.setPixelSize(3.5)

        mdIter = convert.Table.iterRows(tableName + outImgsFn)

        outImgSet.copyItems(imgSet, doClone=doClone,
                            updateItemCallback=self._updateItem,
                            itemDataIterator=mdIter)

        self._defineOutputs(outputParticles=outImgSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        if not self.useAll:
            self._validateDim(self.inputParticles(),
                              self._getInputParticles().getXDim(),
                              errors, 'Input particles subset',
                              'Input particles from 3D protocol')

        return errors
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output is not ready yet.")
        else:
            summary.append('Projections of the masked input volume were '
                           'subtracted from original particles.')

        return summary
    
    # -------------------------- UTILS functions ------------------------------
    def _updateItem(self, particle, row):
        if self.isRelionInput:
            # FIXME: check if other attrs need saving
            particle._rlnRandomSubset = Integer(row.rlnRandomSubset)
            self.reader.setParticleTransform(particle, row)
        particle._rlnImageOriginalName = String(row.rlnImageOriginalName)
        newFn = row.rlnImageName
        newLoc = convert.relionToLocation(newFn)
        particle.setLocation(newLoc)

    def _getInputParticles(self):
        if self.isRelionInput:
            inputProt = self.inputProtocol.get()
            return inputProt.outputParticles
        else:
            return self.inputParticlesAll.get()

    def IS_GT30(self):
        return Plugin.IS_GT30()
