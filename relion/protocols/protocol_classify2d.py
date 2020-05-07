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
from pwem.objects import SetOfClasses2D
from pwem.protocols import ProtClassify2D

import relion.convert as convert
from relion.convert.metadata import Table
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

            if self.IS_GT30() and self.allowCoarserSampling:
                args['--allow_coarser_sampling'] = ''

        else:
            args['--skip_align'] = ''

    # --------------------------- STEPS functions -----------------------------
    def _loadClassesInfo(self, iteration):
        """ Read some information about the produced Relion 2D classes
        from the *model.star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
        modelFn = self._getFileName('model', iter=iteration)
        modelIter = Table.iterRows('model_classes@' + modelFn)
        
        for classNumber, row in enumerate(modelIter):
            index, fn = convert.relionToLocation(row.rlnReferenceImage)
            # Store info indexed by id
            self._classesInfo[classNumber+1] = (index, fn, row)
    
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses2D from a given iteration. """
        self._loadClassesInfo(iteration)
        tableName = 'particles@' if self.IS_GT30() else ''
        dataStar = self._getFileName('data', iter=iteration)
        self.reader = convert.Reader(alignType=pwem.ALIGN_2D)

        mdIter = Table.iterRows(tableName + dataStar, key='rlnImageId')
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
        return []
    
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
            errors += ["You can continue only from the iteration %01d or less" % lastIter]
        
        return errors
    
    def _summaryNormal(self):
        summary = list()
        summary.append("Input Particles: %s" % self.getObjectTag('inputParticles'))
        summary.append("Classified into *%d* classes." % self.numberOfClasses)
        summary.append("Output set: %s" % self.getObjectTag('outputClasses'))
        
        return summary
    
    def _summaryContinue(self):
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
        item.setClassId(row.rlnClassNumber)
        self.reader.setParticleTransform(item, row)

        item._rlnNormCorrection = Float(row.rlnNormCorrection)
        item._rlnLogLikeliContribution = Float(row.rlnLogLikeliContribution)
        item._rlnMaxValueProbDistribution = Float(row.rlnMaxValueProbDistribution)

        if hasattr(row, 'rlnGroupName'):
            item._rlnGroupName = String(row.rlnGroupName)
        
    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            item.setAlignment2D()
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = Float(row.rlnClassDistribution)
            item._rlnAccuracyRotations = Float(row.rlnAccuracyRotations)
            if self.IS_GT30():
                item._rlnAccuracyTranslationsAngst = Float(row.rlnAccuracyTranslationsAngst)
            else:
                item._rlnAccuracyTranslations = Float(row.rlnAccuracyTranslations)
