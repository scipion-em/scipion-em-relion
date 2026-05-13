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
from enum import Enum
from emtable import Table

from pyworkflow.constants import PROD
from pwem.constants import ALIGN_PROJ
from pwem.objects import Volume, FSC, SetOfParticles
from pwem.protocols import ProtRefine3D

import relion.convert as convert
from ..constants import PARTICLE_EXTRA_LABELS
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputVolume = Volume
    outputParticles = SetOfParticles
    outputFSC = FSC


class ProtRelionRefine3D(ProtRefine3D, ProtRelionBase):
    """
    Refines a 3D cryo-EM reconstruction using the Bayesian optimization
    strategy implemented in Relion. The protocol iteratively improves
    particle alignments and map quality in order to generate high-resolution
    density maps together with refined particle orientations and
    reconstruction statistics.

    AI Generated:

    3D Auto-Refine (ProtRelionRefine3D) — User Manual
        Overview

        The 3D Auto-Refine protocol performs high-resolution refinement of
        cryo-EM particle datasets using the empirical Bayesian framework
        implemented in Relion. Its primary goal is to improve both particle
        alignment parameters and the reconstructed density map through an
        iterative optimization process that balances signal recovery with
        noise suppression.

        In practical cryo-EM workflows, this protocol is commonly used after
        obtaining an initial 3D reference volume from ab initio
        reconstruction, stochastic initialization, or low-resolution
        refinement. The refinement progressively improves angular assignments,
        translational alignment, and reconstruction consistency, ultimately
        producing a biologically interpretable map suitable for downstream
        structural analysis.

        Biological Purpose and Context

        For structural biology users, the protocol represents one of the
        central stages of single-particle cryo-EM analysis. It transforms an
        approximate initial reconstruction into a refined consensus structure
        by integrating information from all aligned particles. During this
        process, the refinement attempts to maximize structural consistency
        while minimizing overfitting and noise amplification.

        The resulting map may reveal secondary structure elements, ligand
        densities, conformational features, or interaction interfaces,
        depending on the quality of the data and the achieved resolution.
        Because the refinement is driven directly by experimental particle
        images, the biological reliability of the final structure depends
        strongly on the quality of the input dataset and the correctness of
        the initial reference.

        Inputs and Initial Requirements

        The protocol requires a set of particles with associated alignment
        metadata together with an initial 3D reference map. The reference
        does not need to be perfect, but it should already represent the
        correct overall architecture of the particle. Severe structural
        inaccuracies or incorrect symmetry assumptions may lead the
        refinement toward incorrect solutions.

        Particle quality is equally important. Datasets containing large
        amounts of contamination, aggregation, damaged particles, or strong
        compositional heterogeneity may converge poorly or produce maps with
        limited interpretability. In many workflows, users perform extensive
        2D classification and particle cleaning before beginning 3D
        auto-refinement.

        Symmetry and Structural Assumptions

        The protocol supports symmetry-aware refinement, which can greatly
        improve the achievable resolution when the biological assembly truly
        obeys the specified symmetry. Applying the correct symmetry increases
        the effective number of observations and stabilizes alignment.

        However, incorrect symmetry assignment can introduce strong artifacts
        and misleading structural features. Biological users should therefore
        apply symmetry cautiously, particularly in systems with flexible
        domains, partial occupancy, or pseudo-symmetry. In uncertain cases,
        refinement without symmetry is often safer during exploratory stages.

        The protocol also supports symmetry relaxation strategies that allow
        deviations from perfect symmetry. These approaches can be useful for
        studying flexible or partially asymmetric assemblies while still
        benefiting from the overall symmetry constraints of the complex.

        Angular and Translational Refinement

        During refinement, the protocol searches for the most likely
        orientation and position of each particle projection relative to the
        evolving 3D map. The refinement strategy progressively transitions
        from broader searches toward increasingly fine angular and
        translational sampling as convergence improves.

        Early iterations typically focus on stabilizing the global alignment
        of the dataset, whereas later iterations concentrate on fine local
        adjustments that improve high-resolution detail. The balance between
        sampling precision and computational cost is important because overly
        coarse sampling may limit achievable resolution, while excessively
        fine sampling may increase runtime without meaningful improvement.

        For difficult datasets with significant conformational variability or
        low signal-to-noise ratio, refinement may converge slowly or become
        trapped in local optima. In such cases, improving particle cleaning,
        masking strategies, or reference quality often provides better
        results than simply increasing computational effort.

        Half-Map Strategy and FSC Validation

        A central feature of the protocol is the independent refinement of
        two particle half-sets. This gold-standard refinement strategy is
        designed to reduce overfitting and provide a reliable estimate of the
        true structural resolution.

        The protocol generates two independent half-maps together with a
        Fourier Shell Correlation curve used to estimate resolution. The FSC
        curve is one of the most important outputs because it reflects the
        consistency between the independently refined reconstructions.

        From a biological perspective, the reported resolution should always
        be interpreted together with visual inspection of the map quality.
        Nominal resolution alone does not guarantee that all regions of the
        structure are equally reliable. Flexible domains, membrane regions,
        or peripheral subunits often remain substantially worse resolved than
        the global FSC estimate suggests.

        Solvent Masks and Overfitting Control

        The protocol supports solvent masking strategies that focus the
        refinement on biologically relevant density while excluding solvent
        regions dominated by noise. Proper masking is essential for accurate
        FSC estimation and high-resolution refinement stability.

        A well-designed mask should closely follow the molecular envelope
        while remaining sufficiently soft at the edges to avoid artificial
        correlations. Excessively tight masks may inflate resolution
        estimates, whereas overly loose masks may reduce refinement quality.

        In practice, solvent masking becomes increasingly important at higher
        resolutions where small overfitting artifacts can significantly
        affect interpretation. Biological users should therefore verify that
        masks are physically meaningful and not artificially restrictive.

        Continuing and Extending Refinement

        The protocol supports continuation from previous refinement runs.
        This capability is especially useful when extending an existing
        refinement with improved parameters, additional particles, or more
        computational resources.

        Continuing refinement is common after intermediate particle cleaning,
        focused classification, or updated masking strategies. In many
        workflows, users alternate between refinement and classification
        steps in order to progressively improve homogeneity and map quality.

        Outputs and Biological Interpretation

        The protocol produces a refined 3D reconstruction, two independent
        half-maps, refined particle alignments, and Fourier Shell
        Correlation statistics. The refined particles preserve updated
        orientation and transformation information that can later be reused
        for classification, signal subtraction, focused refinement, or
        visualization.

        The final map should always be interpreted within the biological and
        experimental context of the sample. High nominal resolution does not
        necessarily imply structural correctness in poorly resolved regions,
        and apparent density features should be validated against known
        chemistry, symmetry, and biochemical evidence whenever possible.

        Practical Recommendations

        In routine cryo-EM processing, it is generally advisable to begin
        refinement with a conservative initial model and carefully cleaned
        particles. Excessive heterogeneity or inaccurate references are among
        the most common causes of unstable refinement behavior.

        Monitoring the FSC evolution and visual quality of intermediate maps
        is often more informative than relying solely on numerical metrics.
        If refinement stagnates or produces unrealistic density, users should
        reassess masking, symmetry assignment, particle quality, or angular
        sampling settings.

        For high-resolution studies, biological interpretation should focus
        not only on global resolution values but also on local map quality,
        conformational variability, and reproducibility across independent
        refinements.

        Final Perspective

        For most cryo-EM practitioners, 3D auto-refinement represents the
        critical transition from preliminary reconstruction to a biologically
        meaningful structural model. Careful dataset preparation, realistic
        symmetry assumptions, proper masking, and rigorous validation are the
        essential elements that determine whether the resulting map can
        support reliable structural and mechanistic conclusions.
    """
    _label = '3D auto-refine'
    _devStatus = PROD
    _possibleOutputs = outputs
    IS_CLASSIFY = False
    CHANGE_LABELS = ['rlnChangesOptimalOrientations',
                     'rlnChangesOptimalOffsets',
                     'rlnOverallAccuracyRotations',
                     'rlnOverallAccuracyTranslationsAngst']

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
        if self.relaxSymm.get():
            args['--relax_sym'] = self.relaxSymm.get()

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

            if self.useFinerSamplingFaster:
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

        self._defineOutputs(**{outputs.outputVolume.name: vol})
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(self.inputParticles, outImgSet)

        fsc = FSC(objLabel=self.getRunName())
        fn = self._getExtraPath("relion_model.star")
        table = Table(fileName=fn, tableName='model_class_1')
        resolution_inv = table.getColumnValues('rlnResolution')
        frc = table.getColumnValues('rlnGoldStandardFsc')
        fsc.setData(resolution_inv, frc)

        self._defineOutputs(**{outputs.outputFSC.name: fsc})
        self._defineSourceRelation(vol, fsc)

    # -------------------------- INFO functions -------------------------------
    def _validateNormal(self):
        errors = []

        if self.IS_3D and self.solventFscMask and not self.referenceMask.get():
            errors.append('When using solvent-corrected FSCs, '
                          'please provide a reference mask.')

        if self.numberOfMpi < 3:
            errors.append("3D auto-refine needs at least 3 MPI processes.")

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
            errors.append("You can continue only from the iteration %01d or less" % lastIter)

        if self.numberOfMpi < 3:
            errors.append("3D auto-refine needs at least 3 MPI processes.")
        
        return errors
    
    def _summaryNormal(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume not ready yet.")
            it = self._lastIter() or -1
            if it >= 1 and it > self._getContinueIter():
                table = Table(fileName=self._getFileName('half1_model', iter=it),
                              tableName='model_general')
                row = table[0]
                resol = float(row.rlnCurrentResolution)
                summary.append("Current resolution: *%0.2f A*" % resol)
        else:
            table = Table(fileName=self._getFileName('modelFinal'),
                          tableName='model_general')
            row = table[0]
            resol = float(row.rlnCurrentResolution)
            summary.append("Final resolution: *%0.2f A*" % resol)

        return summary
    
    def _summaryContinue(self):
        return ["Continue from iteration %01d" % self._getContinueIter()]

    # -------------------------- UTILS functions ------------------------------
    def _fillDataFromIter(self, imgSet, iteration):
        outImgsFn = self._getFileName('data', iter=iteration)
        imgSet.setAlignmentProj()
        self.reader = convert.createReader(alignType=ALIGN_PROJ,
                                           pixelSize=imgSet.getSamplingRate())

        mdIter = Table.iterRows('particles@' + outImgsFn, key='rlnImageId',
                                types=convert.LABELS_DICT)
        imgSet.copyItems(self._getInputParticles(), doClone=False,
                         updateItemCallback=self._updateParticle,
                         itemDataIterator=mdIter)

    def _updateParticle(self, particle, row):
        self.reader.setParticleTransform(particle, row)

        if getattr(self, '__updatingFirst', True):
            self.reader.createExtraLabels(particle, row, PARTICLE_EXTRA_LABELS)
            self.__updatingFirst = False
        else:
            self.reader.setExtraLabels(particle, row)
