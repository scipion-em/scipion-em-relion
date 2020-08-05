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
This module contains the protocol for 3d classification with relion.
"""
from os import remove
from os.path import abspath, exists

from pwem.protocols import ProtClassify3D, params
from .protocol_base_tomo import ProtRelionBaseTomo
from tomo.protocols import ProtTomoBase


class ProtRelionSubtomoClassif3D(ProtClassify3D, ProtRelionBaseTomo, ProtTomoBase):
    """
    Protocol to classify 3D using Relion. Relion employs an empirical
    Bayesian approach to refinement of (multiple) 3D reconstructions
    or 2D class averages in electron cryo-microscopy (cryo-EM). Many
    parameters of a statistical model are learned from the data,which
    leads to objective and high-quality results.
    """
    _label = '3D subtomogram classification'
    IS_CLASSIFY = True

    def __init__(self, **args):
        ProtRelionBaseTomo.__init__(self, **args)

    # @classmethod
    # def isDisabled(cls):
    #     return Plugin.IS_30()

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBaseTomo._initialize(self)

    # --------------------------- INSERT steps functions --------------------------------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params. """
        if self.doImageAlignment.get():
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * 2
            if self.localAngularSearch.get():
                args['--sigma_ang'] = self.localAngularSearchRange.get() / 3.
        else:
            args['--skip_align'] = ''

    # --------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):
        subtomoSet = self._getInputParticles()
        classes3D = self._createSetOfClassesSubTomograms(subtomoSet)
        self._fillClassesFromIter(classes3D, self._lastIter())

        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(subtomoSet, classes3D)

        # Create a SetOfVolumes and define its relations
        volumes = self._createSetOfAverageSubTomograms()
        volumes.setSamplingRate(subtomoSet.getSamplingRate())

        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(subtomoSet, volumes)

        if not self.doContinue:
            self._defineSourceRelation(self.referenceVolume, classes3D)
            self._defineSourceRelation(self.referenceVolume, volumes)

        if self.keepOnlyLastIterFiles:
            self._cleanUndesiredFiles()

    # --------------------------- INFO functions --------------------------------------------
    def _validateNormal(self):
        """ Should be overriden in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        errors = []
        if self.referenceVolume.get() is not None:
            partSizeX, _, _ = self._getInputParticles().getDim()
            volSizeX, _, _ = self.referenceVolume.get().getDim()
            if partSizeX != volSizeX:
                errors.append('Volume and particles dimensions must be equal!!!')

        return errors

    def _validateContinue(self):
        """ Should be overriden in subclasses to
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
            errors += ["The iteration from you want to continue must be %01d or less" % lastIter]

        return errors

    def _summaryNormal(self):
        """ Should be overriden in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        summary = []
        it = self._lastIter()
        if it:
            if it >= 1:
                modelStar = self._getFileName('model', iter=it)
                resol = self._getRlnCurrentResolution(modelStar)
                summary.append("Current resolution: *%0.2f* Å" % resol)
                summary.append("Input Particles: *%d*\nClassified into *%d* 3D classes\n" %
                               (self.inputSubtomograms.get().getSize(), self.numberOfClasses.get()))
        return summary

    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        summary = []
        summary.append("Continue from iteration %01d" % self._getContinueIter())
        return summary

    def _methods(self):
        strline = ''
        if hasattr(self, 'outputClasses'):
            strline += 'We classified %d particles into %d 3D classes using Relion Classify3d. ' % \
                       (self.inputParticles.get().getSize(), self.numberOfClasses.get())
        return [strline]

    # --------------------------- UTILS functions --------------------------------------------
    def _loadClassifyInfo(self, iteration):
        """ Read some information about the produced Relion 3D classes
        from the *model.star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id
        modelStar = self._getFileName('model', iter=iteration)
        with open(modelStar) as fid:
            self.modelTable.readStar(fid, 'model_general')
            self.claasesTable.readStar(fid, 'model_classes')
        dataStar = self._getFileName('data', iter=iteration)
        with open(dataStar) as fid:
            self.dataTable.readStar(fid)

        # Model table has only one row, while classes table has the same number of rows as classes found
        self.nClasses = int(self.modelTable._rows[0].rlnNrClasses)
        # Adapt data to the format expected by classifyItems and its callbacks
        for i, row in zip(range(self.nClasses), self.claasesTable.__iter__()):
            self._classesInfo[i + 1] = (i + 1, row)

    def _fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses3D from a given iteration. """
        self._loadClassifyInfo(iteration)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=self.dataTable.__iter__())

    def _updateParticle(self, item, row):
        item.setClassId(int(row.rlnClassNumber))#rlnGroupNumber))
        item._rlnLogLikeliContribution = params.Float(row.rlnLogLikeliContribution)
        item._rlnMaxValueProbDistribution = params.Float(row.rlnMaxValueProbDistribution)

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, row = self._classesInfo[classId]
            fn = row.rlnReferenceImage + ":mrc"
            item.setAlignment3D()
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = params.Float(row.rlnClassDistribution)
            item._rlnAccuracyRotations = params.Float(row.rlnAccuracyRotations)
            item._rlnAccuracyTranslations = params.Float(row.rlnAccuracyTranslations)

    def _cleanUndesiredFiles(self):
        """Remove all files generated by relion_classify 3d excepting the ones which
        correspond to the last iteration. Example for iteration 25:
        relion_it025_class002.mrc
        relion_it025_class001.mrc
        relion_it025_model.star
        relion_it025_sampling.star
        relion_it025_optimiser.star
        relion_it025_data.star
        """
        itPref = 'relion_it'
        clPref = 'class'
        starExt = '.star'
        mrcExt = '.mrc'
        # Classify calculations related files
        calcFiles = ['data', 'model', 'optimiser', 'sampling']
        for i in range(self._lastIter()):
            for calcFile in calcFiles:
                fn = abspath(self._getExtraPath('{}{:03d}_{}{}'.format(
                    itPref, i, calcFile, starExt)))
                if exists(fn):
                    remove(fn)
            # Classes related files
            for itr in range(1, self.nClasses + 1):
                fn = abspath(self._getExtraPath('{}{:03d}_{}{:03d}{}'.format(
                    itPref, i, clPref, itr, mrcExt)))
                if exists(fn):
                    remove(fn)

