# **************************************************************************
# *
# * Author:     Scipion Team
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
# *  e-mail address scipion-users@lists.sourceforge.net
# *
# **************************************************************************
"""
This module contains the protocol for 3d refinement with Relion.
"""
from os import remove
from os.path import abspath, exists

from pwem.objects import Volume, ALIGN_PROJ
from pwem.protocols import ProtRefine3D
from relion.convert import createItemMatrix
from relion.protocols.protocol_base_tomo import ProtRelionBaseTomo
from relion import Plugin


class ProtRelionSubtomoRefine3D(ProtRefine3D, ProtRelionBaseTomo):
    """Protocol to refine a 3D map using Relion. Relion employs an empirical
Bayesian approach to refinement of (multiple) 3D reconstructions
or 2D class averages in electron cryo-microscopy (cryo-EM). Many
parameters of a statistical model are learned from the data,which
leads to objective and high-quality results.
    """
    _label = '3D subtomogram auto-refine'
    IS_CLASSIFY = False

    PREFIXES = ['half1_', 'half2_']

    def __init__(self, **args):
        ProtRelionBaseTomo.__init__(self, **args)

    @classmethod
    def isDisabled(cls):
        return Plugin.IS_30()

    def _initialize(self):
        """ This function is mean to be called after the
        working dir for the protocol have been set. (maybe after recovery from mapper)
        """
        ProtRelionBaseTomo._initialize(self)
        self.ClassFnTemplate = '%(ref)03d@%(rootDir)s/relion_it%(iter)03d_classes.mrcs'
        self.numberOfClasses.set(1)  # For refinement we only need just one "class"

    # --------------------------- INSERT steps functions --------------------------------------------
    def _setSamplingArgs(self, args):
        """ Set sampling related params"""
        # Sampling stuff
        args['--auto_local_healpix_order'] = self.localSearchAutoSamplingDeg.get()

        if not self.doContinue:
            args['--healpix_order'] = self.angularSamplingDeg.get()
            args['--offset_range'] = self.offsetSearchRangePix.get()
            args['--offset_step'] = self.offsetSearchStepPix.get() * 2
            args['--auto_refine'] = ''
            args['--split_random_halves'] = ''

            joinHalves = "--low_resol_join_halves"
            if not joinHalves in self.extraParams.get():
                args['--low_resol_join_halves'] = 40

    # --------------------------- STEPS functions --------------------------------------------
    def createOutputStep(self):

        imgSet = self._getInputParticles()
        vol = Volume()
        vol.setFileName(self._getExtraPath('relion_refined_volume.mrc'))
        vol.setSamplingRate(imgSet.getSamplingRate())

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(imgSet, vol)

        if self.keepOnlyLastIterFiles:
            self._cleanUndesiredFiles()

    # --------------------------- INFO functions --------------------------------------------
    def _validateNormal(self):
        """ Should be overriden in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        errors = []
        if self.referenceVolume.get() is not None:
            particlesDim = self._getInputParticles().getDim()
            volumeDim = self.referenceVolume.get().getDim()

            if particlesDim is None:
                errors.append('Can not get dimensions from input particles!!!')

            elif volumeDim is None:
                errors.append('Can not get dimensions from reference volume!!!')

            elif particlesDim[0] != volumeDim[0]:
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
                modelStar = self._getFileName('half1_model', iter=it)
                summary.append("Current resolution %s: *%0.2f* Å" %
                               self._getRlnCurrentResolution(modelStar))
        return summary

    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        summary = ["Continue from iteration %01d" % self._getContinueIter()]
        return summary

    # --------------------------- UTILS functions --------------------------------------------
    @ staticmethod
    def _createItemMatrix(item, row):
        createItemMatrix(item, row, align=ALIGN_PROJ)

    def _cleanUndesiredFiles(self):
        """Remove all files generated by relion_refine 3d excepting the ones which
        correspond to the last iteration. Example for iteration 12:
        relion_it012_half1_class001.mrc
        relion_it012_half2_class001_angdist.bild
        relion_it012_half2_class001.mrc
        relion_it012_half2_model.star
        relion_it012_optimiser.star
        relion_it012_half1_class001_angdist.bild
        relion_it012_data.star
        relion_it012_sampling.star
        relion_it012_half1_model.star
        """
        itPref = 'relion_it'
        clPref = 'class'
        starExt = '.star'
        mrcExt = '.mrc'
        bildExt = '_angdist.bild'
        halfPref = 'half'
        # Refine calculations related files
        calcFiles = ['data', '%s1_model' % halfPref, '%s2_model' % halfPref, 'optimiser', 'sampling']
        for i in range(self._lastIter()):
            for calcFile in calcFiles:
                fn = abspath(self._getExtraPath('{}{:03d}_{}{}'.format(
                    itPref, i, calcFile, starExt)))
                if exists(fn):
                    remove(fn)
            # Classes related files
            for itr in range(1, self.nClasses + 1):
                for h in range(1, 3):  # Two halves
                    fn = abspath(self._getExtraPath('{}{:03d}_{}{}_{}{:03d}'.format(
                        itPref, i, halfPref, h, clPref, itr)))
                    mrcFile = fn + mrcExt
                    bildFile = fn + bildExt
                    if exists(mrcFile):
                        remove(mrcFile)
                    if exists(bildFile):
                        remove(bildFile)
