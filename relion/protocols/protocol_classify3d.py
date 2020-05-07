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

import relion
from relion import Plugin
import relion.convert as convert
from relion.convert.metadata import Table
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
    CHANGE_LABELS = ['rlnChangesOptimalOrientations',
                     'rlnChangesOptimalOffsets',
                     'rlnOverallAccuracyRotations',
                     'rlnOverallAccuracyTranslationsAngst' if Plugin.IS_GT30() else 'rlnOverallAccuracyTranslations',
                     'rlnChangesOptimalClasses']
    
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
            errors += ["You can continue only from the iteration %01d or less" % lastIter]
        
        return errors
    
    def _summaryNormal(self):
        summary = []
        it = self._lastIter()
        if it >= 1:
            table = Table(fileName=self._getFileName('model', iter=it),
                          tableName='model_general')
            row = table[0]
            resol = float(row.rlnCurrentResolution)
            summary.append("Current resolution: *%0.2f A*" % resol)
        
        summary.append("Input Particles: *%d*\n"
                       "Classified into *%d* 3D classes\n"
                       % (self.inputParticles.get().getSize(),
                          self.numberOfClasses.get()))
        
        return summary
    
    def _summaryContinue(self):
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

        modelFn = self._getFileName('model', iter=iteration)
        modelIter = Table.iterRows('model_classes@' + modelFn)

        for classNumber, row in enumerate(modelIter):
            index, fn = convert.relionToLocation(row.rlnReferenceImage)
            # Store info indexed by id
            self._classesInfo[classNumber+1] = (index, fn, row)
    
    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses3D from a given iteration. """
        self._loadClassesInfo(iteration)
        tableName = 'particles@' if self.IS_GT30() else ''
        dataStar = self._getFileName('data', iter=iteration)
        self.reader = convert.Reader(alignType=pwem.ALIGN_PROJ)
        mdIter = Table.iterRows(tableName + dataStar, key='rlnImageId')
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=mdIter,
                             doClone=False)
    
    def _updateParticle(self, item, row):
        item.setClassId(row.rlnClassNumber)
        self.reader.setParticleTransform(item, row)

        item._rlnLogLikeliContribution = Float(row.rlnLogLikeliContribution)
        item._rlnMaxValueProbDistribution = Float(row.rlnMaxValueProbDistribution)

        if hasattr(row, 'rlnGroupName'):
            item._rlnGroupName = String(row.rlnGroupName)

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            fn += ":mrc"
            item.setAlignmentProj()
            item.getRepresentative().setLocation(index, fn)
            item._rlnClassDistribution = Float(row.rlnClassDistribution)
            item._rlnAccuracyRotations = Float(row.rlnAccuracyRotations)
            if self.IS_GT30():
                item._rlnAccuracyTranslationsAngst = Float(row.rlnAccuracyTranslationsAngst)
            else:
                item._rlnAccuracyTranslations = Float(row.rlnAccuracyTranslations)
