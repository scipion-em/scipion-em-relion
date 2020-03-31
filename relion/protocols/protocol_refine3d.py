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

from pyworkflow.object import String

import pwem
import pwem.emlib.metadata as md
from pwem.objects import Volume, FSC
from pwem.protocols import ProtRefine3D

import relion
import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionRefine3D(ProtRefine3D, ProtRelionBase):
    """ Protocol to refine a 3D map using Relion.

Relion employs an empirical Bayesian approach to refinement
of (multiple) 3D reconstructions
or 2D class averages in electron cryo-microscopy (cryo-EM). Many
parameters of a statistical model are learned from the data,which
leads to objective and high-quality results.
    """    
    _label = '3D auto-refine'
    IS_CLASSIFY = False
    CHANGE_LABELS = [md.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, 
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS,
                     md.RLN_OPTIMISER_ACCURACY_ROT,
                     md.RLN_OPTIMISER_ACCURACY_TRANS]

    PREFIXES = ['half1_', 'half2_']
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'

    # -------------------------- INSERT steps functions -----------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params"""
        args['--auto_local_healpix_order'] = self.localSearchAutoSamplingDeg.get()
        
        if not self.doContinue:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            f = self._getSamplingFactor()
            args['--offset_step'] = self.offsetSearchStepPix.get() * f
            args['--auto_refine'] = ''
            args['--split_random_halves'] = ''
            
            joinHalves = "--low_resol_join_halves"
            if joinHalves not in self.extraParams.get():
                args['--low_resol_join_halves'] = 40

            if relion.Plugin.IS_GT30() and self.useFinerSamplingFaster:
                args['--auto_ignore_angles'] = ''
                args['--auto_resol_angles'] = ''

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        imgSet = self._getInputParticles()
        vol = Volume()
        vol.setFileName(self._getExtraPath('relion_class001.mrc'))
        vol.setSamplingRate(imgSet.getSamplingRate())
        half1 = self._getFileName("final_half1_volume", ref3d=1)
        half2 = self._getFileName("final_half2_volume", ref3d=1)
        vol.setHalfMaps([half1, half2])

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet, self._lastIter())

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

        fsc = FSC(objLabel=self.getRunName())
        blockName = 'model_class_%d@' % 1
        fn = blockName + self._getExtraPath("relion_model.star")
        mData = md.MetaData(fn)
        fsc.loadFromMd(mData,
                       md.RLN_RESOLUTION,
                       md.RLN_MLMODEL_FSC_HALVES_REF)
        self._defineOutputs(outputFSC=fsc)
        self._defineSourceRelation(vol, fsc)

    # -------------------------- INFO functions -------------------------------
    def _validateNormal(self):
        """ Should be overwritten in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        # We we scale the input volume to have the same size as the particles...
        # so no need to validate the following
        # self._validateDim(self._getInputParticles(), self.referenceVolume.get(),
        #                   errors, 'Input particles', 'Reference volume')

        if self.IS_3D and self.solventFscMask and not self.referenceMask.get():
            errors.append('When using solvent-corrected FSCs, '
                          'please provide a reference mask.')

        return errors
    
    def _validateContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        errors = []
        continueRun = self.continueRun.get()
        continueRun._initialize()
        lastIter = continueRun._lastIter()
        
        if self.continueIter.get() == 'last':
            continueIter = lastIter
        else:
            continueIter = int(self.continueIter.get())
        
        if continueIter > lastIter:
            errors += ["You can continue only from the iteration %01d or less" % lastIter]
        
        return errors
    
    def _summaryNormal(self):
        """ Should be overwritten in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume not ready yet.")
            it = self._lastIter()
            if it >= 1 and it > self._getContinueIter():
                row = md.getFirstRow('model_general@' +
                                     self._getFileName('half1_model',
                                                       iter=it))
                resol = row.getValue("rlnCurrentResolution")
                summary.append("Current resolution: *%0.2f A*" % resol)
        else:
            row = md.getFirstRow('model_general@' +
                                 self._getFileName('modelFinal'))
            resol = row.getValue("rlnCurrentResolution")
            summary.append("Final resolution: *%0.2f A*" % resol)

        return summary
    
    def _summaryContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return ["Continue from iteration %01d" % self._getContinueIter()]

    # -------------------------- UTILS functions ------------------------------
    def _fillDataFromIter(self, imgSet, iteration):
        tableName = '' if relion.Plugin.IS_30() else 'particles@'
        outImgsFn = self._getFileName('data', iter=iteration)
        imgSet.setAlignmentProj()
        self.reader = convert.Reader(alignType=pwem.ALIGN_PROJ)
        mdIter = md.iterRows(tableName + outImgsFn, sortByLabel=md.RLN_IMAGE_ID)
        imgSet.copyItems(self._getInputParticles(), doClone=False,
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=mdIter)

    def _createItemMatrix(self, particle, row):
        self.reader.setParticleTransform(particle, row)
        convert.setRelionAttributes(particle, row,
                                    md.RLN_PARTICLE_RANDOM_SUBSET)

    def _updateParticle(self, particle, row):
        particle._coordinate._micName = String(row.getValue('rlnMicrographName'))
