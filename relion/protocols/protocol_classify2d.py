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
import pwem.emlib.metadata as md
from pwem.objects import SetOfClasses2D
from pwem.protocols import ProtClassify2D

import relion
import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionClassify2D(ProtRelionBase, ProtClassify2D):
    """ This protocol runs Relion 2D classification."""

    _label = '2D classification'
    IS_2D = True
    OUTPUT_TYPE = SetOfClasses2D
    
    def __init__(self, **args):        
        ProtRelionBase.__init__(self, **args)
        
    def _initialize(self):
        """ This function is mean to be called after the 
        working dir for the protocol have been set.
        (maybe after recovery from mapper)
        """
        ProtRelionBase._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'

    # --------------------------- INSERT steps functions ----------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        if self.doImageAlignment:
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * self._getSamplingFactor()
            args['--psi_step'] = self.inplaneAngularSamplingDeg.get() * self._getSamplingFactor()

            if relion.Plugin.IS_GT30() and self.allowCoarserSampling:
                args['--allow_coarser_sampling'] = ''

        else:
            args['--skip_align'] = ''

    # --------------------------- STEPS functions -----------------------------
    def _loadClassesInfo(self, iteration):
        """ Read some information about the produced Relion 2D classes
        from the *model.star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
         
        modelStar = md.MetaData('model_classes@%s' %
                                self._getFileName('model', iter=iteration))
        
        for classNumber, row in enumerate(md.iterRows(modelStar)):
            index, fn = convert.relionToLocation(row.getValue('rlnReferenceImage'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration            
            self._classesInfo[classNumber+1] = (index, fn, row.clone())
    
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses2D from a given iteration. """
        self._loadClassesInfo(iteration)
        tableName = '' if relion.Plugin.IS_30() else 'particles@'
        dataStar = self._getFileName('data', iter=iteration)
        self.reader = convert.Reader(alignType=pwem.ALIGN_2D)
        mdIter = md.iterRows(tableName + dataStar, sortByLabel=md.RLN_IMAGE_ID)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=mdIter,
                             doClone=False)
        
    def createOutputStep(self):
        partSet = self.inputParticles.get()       
        
        classes2D = self._createSetOfClasses2D(partSet)
        self._fillClassesFromIter(classes2D, self._lastIter())
        
        self._defineOutputs(outputClasses=classes2D)
        self._defineSourceRelation(self.inputParticles, classes2D)
        
    # --------------------------- INFO functions ------------------------------
    def _validateNormal(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
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
        summary = list()
        summary.append("Input Particles: %s" % self.getObjectTag('inputParticles'))
        summary.append("Classified into *%d* classes." % self.numberOfClasses)
        summary.append("Output set: %s" % self.getObjectTag('outputClasses'))
        
        return summary
    
    def _summaryContinue(self):
        """ Should be overwritten in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        summary = list()
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        
        return summary
    
    def _methods(self):
        methods = ''
        if hasattr(self, 'outputClasses'):
            methods += "We classified input particles %s (%d items) " % (
                self.getObjectTag('inputParticles'),
                self.inputParticles.get().getSize())
            methods += "into %d classes using Relion Classify2d. " % self.numberOfClasses
            methods += 'Output classes: %s' % self.getObjectTag('outputClasses')
        return [methods]
    
    # --------------------------- UTILS functions -----------------------------
    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        self.reader.setParticleTransform(item, row)

        item._rlnNormCorrection = Float(row.getValue('rlnNormCorrection'))
        item._rlnLogLikeliContribution = Float(row.getValue('rlnLogLikeliContribution'))
        item._rlnMaxValueProbDistribution = Float(row.getValue('rlnMaxValueProbDistribution'))
        item._rlnGroupName = String(row.getValue('rlnGroupName'))
        
    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            item.setAlignment2D()
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = Float(row.getValue('rlnClassDistribution'))
            item._rlnAccuracyRotations = Float(row.getValue('rlnAccuracyRotations'))
            item._rlnAccuracyTranslations = Float(row.getValue('rlnAccuracyTranslations'))
