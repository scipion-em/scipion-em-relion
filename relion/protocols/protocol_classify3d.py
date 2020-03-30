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

from pyworkflow.object import String, Float
import pwem
from pwem.protocols import ProtClassify3D
import pwem.emlib.metadata as md

import relion
import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionClassify3D(ProtClassify3D, ProtRelionBase):
    """
    Protocol to classify 3D using Relion Bayesian approach.
    Relion employs an empirical Bayesian approach to refinement of (multiple)
    3D reconstructions or 2D class averages in electron cryo-EM. Many
    parameters of a statistical model are learned from the data, which
    leads to objective and high-quality results.
    """

    _label = '3D classification'
    CHANGE_LABELS = [md.RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, 
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS,
                     md.RLN_OPTIMISER_ACCURACY_ROT,
                     md.RLN_OPTIMISER_ACCURACY_TRANS,
                     md.RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES]
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
    
    # -------------------------- INSERT steps functions -----------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        if self.doImageAlignment:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * self._getSamplingFactor()

            if self.localAngularSearch:
                args['--sigma_ang'] = self.localAngularSearchRange.get() / 3.

            if relion.Plugin.IS_GT30() and self.allowCoarserSampling:
                args['--allow_coarser_sampling'] = ''

        else:
            args['--skip_align'] = ''
    
    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        partSet = self.inputParticles.get()
        classes3D = self._createSetOfClasses3D(partSet)
        self._fillClassesFromIter(classes3D, self._lastIter())
        
        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(self.inputParticles, classes3D)

        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(partSet.getSamplingRate())
        
        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)
        
        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputParticles, volumes)
        
        if not self.doContinue:
            self._defineSourceRelation(self.referenceVolume, classes3D)
            self._defineSourceRelation(self.referenceVolume, volumes)
    
    # -------------------------- INFO functions -------------------------------
    def _validateNormal(self):
        """ Should be overwritten in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
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
        it = self._lastIter()
        if it >= 1:
            row = md.getFirstRow('model_general@' + self._getFileName('model', iter=it))
            resol = row.getValue("rlnCurrentResolution")
            summary.append("Current resolution: *%0.2f A*" % resol)
        
        summary.append("Input Particles: *%d*\n"
                       "Classified into *%d* 3D classes\n"
                       % (self.inputParticles.get().getSize(),
                          self.numberOfClasses.get()))
        
        return summary
    
    def _summaryContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        summary = list()
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        return summary
    
    def _methods(self):
        strline = ''
        if hasattr(self, 'outputClasses'):
            strline += 'We classified %d particles into %d 3D classes using Relion Classify3d. ' %\
                           (self.inputParticles.get().getSize(), self.numberOfClasses.get())
        return [strline]
    
    # -------------------------- UTILS functions ------------------------------
    def _loadClassesInfo(self, iteration):
        """ Read some information about the produced Relion 3D classes
        from the *model.star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
         
        modelStar = md.MetaData('model_classes@' +
                                self._getFileName('model', iter=iteration))
        
        for classNumber, row in enumerate(md.iterRows(modelStar)):
            index, fn = convert.relionToLocation(row.getValue('rlnReferenceImage'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration            
            self._classesInfo[classNumber+1] = (index, fn, row.clone())
    
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses3D from a given iteration. """
        self._loadClassesInfo(iteration)
        tableName = '' if relion.Plugin.IS_30() else 'particles@'
        dataStar = self._getFileName('data', iter=iteration)
        self.reader = convert.Reader(alignType=pwem.ALIGN_PROJ)
        mdIter = md.iterRows(tableName + dataStar, sortByLabel=md.RLN_IMAGE_ID)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=mdIter,
                             doClone=False)
    
    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        self.reader.setParticleTransform(item, row)

        item._rlnLogLikeliContribution = Float(row.getValue('rlnLogLikeliContribution'))
        item._rlnMaxValueProbDistribution = Float(row.getValue('rlnMaxValueProbDistribution'))
        item._rlnGroupName = String(row.getValue('rlnGroupName'))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            fn += ":mrc"
            item.setAlignmentProj()
            item.getRepresentative().setLocation(index, fn)
            item._rlnClassDistribution = Float(row.getValue('rlnClassDistribution'))
            item._rlnAccuracyRotations = Float(row.getValue('rlnAccuracyRotations'))
            item._rlnAccuracyTranslations = Float(row.getValue('rlnAccuracyTranslations'))
