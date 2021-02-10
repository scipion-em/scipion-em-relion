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

from pyworkflow.object import String, Integer
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam
from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtOperateParticles

import relion.convert as convert
from relion import Plugin


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
        self._updateFilenamesDict(myDict)
    
    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputProtocol', PointerParam,
                      important=True,
                      pointerClass='ProtRelionRefine3D, ProtRelionClassify3D',
                      label="Input 3D protocol",
                      help="Select the 3D refinement/classification run which "
                           "you want to use for subtraction. It will use the "
                           "maps from this run for the subtraction.")

        form.addParam('useAll', BooleanParam, default=True,
                      label="Use all particles from input protocol?",
                      help="If No, then you need to provide a subset of "
                           "particles below.")

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not useAll',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles subset",
                      help='Select the particles which are a SUBSET of the '
                           'input protocol provided above.')

        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask of the signal to keep',
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is black (0) and the density you "
                           "wish to keep is white (1).\n"
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to KEEP.*")

        form.addSection('Centering')
        form.addParam('centerOnMask', BooleanParam, default=True,
                      label="Do center subtracted images on mask?",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the "
                           "center-of-mass of the input mask.")
        form.addParam('centerOnCoord', BooleanParam, default=False,
                      condition='not centerOnMask',
                      label="Do center on my coordinates?",
                      help="If set to Yes, the subtracted particles will "
                           "be centered on projections of the x,y,z "
                           "coordinates below. The unit is pixel, not "
                           "angstrom. The origin is at the center of the box, "
                           "not at the corner.")

        line = form.addLine('Center coordinate (px)',
                            condition='centerOnCoord',
                            help='Coordinate of the 3D center (in pixels).')
        line.addParam('cX', IntParam, default=0, condition='centerOnCoord',
                      label='X')
        line.addParam('cY', IntParam, default=0, condition='centerOnCoord',
                      label='Y')
        line.addParam('cZ', IntParam, default=0, condition='centerOnCoord',
                      label='Z')

        form.addParam('newBoxSize', IntParam, default=-1,
                      label="New box size",
                      help="Provide a non-negative value to re-window the "
                           "subtracted particles in a smaller box size.")

        form.addParallelSection(threads=0, mpi=1)
    
    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initialize()

        if not self.useAll:
            self._insertFunctionStep('convertInputStep')

        self._insertFunctionStep('subtractStep')
        self._insertFunctionStep('createOutputStep')
    
    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ Write the input images as a Relion star file. """
        imgSet = self.inputParticles.get()
        convert.writeSetOfParticles(
            imgSet, self._getFileName('input_star'),
            outputDir=self._getExtraPath(), alignType=ALIGN_PROJ)
    
    def subtractStep(self):
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
        mdIter = convert.Table.iterRows(tableName + outImgsFn)
        outImgSet.copyItems(imgSet, doClone=False,
                            updateItemCallback=self._updateItem,
                            itemDataIterator=mdIter)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)
    
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
        self.reader.setParticleTransform(particle, row)
        # FIXME: check if other attrs need saving
        particle._rlnImageOriginalName = String(row.rlnImageOriginalName)
        particle._rlnRandomSubset = Integer(row.rlnRandomSubset)

        newLoc = convert.relionToLocation(row.rlnImageName)
        particle.setLocation(newLoc)

    def _getInputParticles(self):
        inputProt = self.inputProtocol.get()
        return inputProt.outputParticles

    def IS_GT30(self):
        return Plugin.IS_GT30()
