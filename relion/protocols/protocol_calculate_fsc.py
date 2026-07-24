# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

import os
from emtable import Table

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtAnalysis3D
from pwem.objects import FSC, SetOfFSCs
from pwem.emlib.image import ImageHandler

from ..constants import (FSC_TYPE_OVERALL, FSC_TYPE_MODEL_MAP,
                         FSC_TYPE_WORK_FREE)
from .protocol_base import ProtRelionBase


class ProtRelionCalculateFSC(ProtAnalysis3D, ProtRelionBase):
    """
    Calculates Fourier Shell Correlation (FSC) curves between cryo-EM maps
    and atomic models using Relion image analysis tools. The protocol
    supports multiple validation strategies for estimating the consistency,
    reproducibility, and potential overfitting of reconstructed structures.

    AI Generated:

    Calculate FSC (ProtRelionCalculateFSC) — User Manual
        Overview

        The Calculate FSC protocol evaluates the agreement between cryo-EM
        reconstructions or between reconstructed maps and atomic models using
        Fourier Shell Correlation analysis. FSC is one of the most widely used
        validation metrics in structural biology because it provides a
        resolution-dependent estimate of similarity between two datasets.

        In practical cryo-EM workflows, FSC calculations are essential for
        validating reconstruction quality, estimating map resolution, and
        assessing whether an atomic model faithfully represents the experimental
        density. The protocol is designed to support both standard refinement
        validation and more advanced procedures aimed at detecting overfitting.

        Biological Purpose of FSC Analysis

        FSC analysis compares two independent representations of structural
        information across spatial frequencies. High FSC values indicate strong
        agreement, while decreasing correlation at higher frequencies reflects
        the loss of reliable signal at finer structural detail.

        For biological interpretation, FSC curves are commonly used to estimate
        the effective resolution of a reconstruction. They also help determine
        whether observed structural features are supported by experimental data
        or may instead arise from noise, masking artifacts, or refinement bias.

        The protocol supports several biologically meaningful FSC workflows,
        each adapted to different stages of cryo-EM analysis and validation.

        FSC Overall: Half-Map Consistency

        The FSC overall calculation measures the correlation between two
        independently refined half-maps. This is the standard procedure used
        in modern single-particle cryo-EM to estimate global map resolution.

        Biologically, this analysis evaluates the reproducibility of structural
        information recovered from independent subsets of particles. Strong
        agreement between half-maps indicates that the reconstruction contains
        reproducible signal rather than random noise.

        This mode is most appropriate after refinement or classification steps
        when independent half-maps are already available. The resulting FSC
        curve is typically used to determine the nominal resolution of the map
        using criteria such as the 0.143 threshold.

        FSC Model-Map: Agreement Between Structure and Density

        The FSC model-map mode evaluates the agreement between an atomic model
        and a reconstructed cryo-EM density map. This analysis is especially
        important after atomic model building and refinement.

        From a biological perspective, this calculation helps determine whether
        the structural model is genuinely supported by experimental density.
        Good agreement indicates that the atomic coordinates are consistent
        with observed structural features, while poor agreement may reveal
        local modeling errors, incorrect conformations, or overinterpretation
        of noisy regions.

        This mode is particularly valuable for validating flexible domains,
        ligand placement, membrane protein regions, or areas with heterogeneous
        local resolution. It is often used before deposition or publication as
        part of the final validation workflow.

        FSC Work and FSC Free: Detecting Overfitting

        The FSC work and FSC free calculations are designed to evaluate possible
        overfitting during model refinement. In this strategy, an atomic model
        refined against one half-map is independently compared against both the
        same half-map and the opposite half-map.

        FSC work reflects agreement between the model and the map used during
        refinement, while FSC free measures agreement against independent data.
        If FSC work remains substantially higher than FSC free at high spatial
        frequencies, this may indicate overfitting of noise rather than genuine
        structural signal.

        Biologically, this validation is especially important for high-resolution
        studies where refinement procedures can accidentally introduce features
        unsupported by experimental evidence. Careful interpretation of FSC work
        and FSC free curves helps ensure that reported atomic details are
        reliable and reproducible.

        Inputs and Data Preparation

        Depending on the selected validation strategy, the protocol accepts
        reconstructed maps, half-maps, or atomic coordinate models. The quality
        and consistency of these inputs strongly influence the interpretability
        of the resulting FSC curves.

        Half-maps should originate from truly independent refinements to avoid
        artificially inflated correlations. Similarly, atomic models should be
        refined carefully and should not contain unrealistic geometry or
        unsupported flexible conformations.

        When atomic models are used, they are converted into volumetric density
        representations compatible with FSC calculations. To obtain meaningful
        results, voxel size, box dimensions, and map sampling should remain
        consistent across all inputs.

        Interpretation of FSC Curves

        FSC curves must always be interpreted within their biological and
        experimental context. Smooth curves that gradually decay with frequency
        generally indicate stable and reliable reconstructions, whereas abrupt
        fluctuations or unusually elevated correlations may suggest masking
        artifacts, insufficient sampling, or refinement bias.

        Resolution thresholds should not be treated as absolute indicators of
        biological quality. Local flexibility, conformational heterogeneity,
        preferred particle orientation, and anisotropic resolution can all
        influence FSC behavior. Consequently, FSC analysis should be combined
        with visual inspection of maps and independent biological validation.

        Outputs and Their Use

        The protocol produces one or more FSC curves depending on the selected
        validation strategy. These outputs can be visualized, compared, or used
        as part of downstream validation and reporting workflows.

        In routine cryo-EM practice, FSC curves are commonly included in
        publications, structure depositions, and quality assessment reports.
        They provide a concise summary of reconstruction reproducibility and
        model consistency across spatial frequencies.

        Practical Recommendations

        For routine map resolution estimation, FSC overall calculations between
        independently refined half-maps are usually sufficient. When validating
        atomic models, model-map FSC analysis provides a direct assessment of
        structural consistency with the experimental density.

        In high-resolution studies or when extensive model refinement has been
        performed, FSC work and FSC free calculations are strongly recommended
        to detect potential overfitting. Large discrepancies between the two
        curves should prompt additional validation and careful visual inspection.

        Biological users should also remember that FSC values are sensitive to
        masking strategies, preprocessing operations, and map sharpening.
        Consistent preparation procedures and conservative interpretation are
        essential for obtaining biologically meaningful conclusions.

        Final Perspective

        FSC analysis is a central component of modern cryo-EM validation because
        it links structural interpretation to measurable experimental agreement.
        Whether estimating map resolution, validating atomic models, or checking
        for overfitting, careful FSC analysis provides critical confidence in
        the biological conclusions derived from cryo-EM reconstructions.
    """
    _label = 'calculate fsc'
    _devStatus = PROD
    _possibleOutputs = {
        'outputFSC': FSC,
        'outputSetOfFSCs': SetOfFSCs
    }

    def _createFilenameTemplates(self):
        """ Centralize how the files are called. """
        myDict = {
            'model': self._getTmpPath("input_model_final.mrc"),
            'model_half1': self._getTmpPath("input_model_half1.mrc"),
            'half1': self._getTmpPath("input_half1_unfil.mrc"),
            'half2': self._getTmpPath("input_half2_unfil.mrc"),
            'map': self._getTmpPath("input_map.mrc"),
            # outputs
            'fsc_overall': self._getExtraPath("fsc_overall.star"),
            'fsc_model-map': self._getExtraPath("fsc_model-map.star"),
            'fsc_work': self._getExtraPath("fsc_work.star"),
            'fsc_free': self._getExtraPath("fsc_free.star")
        }
        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('fscType', params.EnumParam, default=0,
                      choices=['FSC overall', 'FSC model-map',
                               'FSC work / FSC free'],
                      label="Select FSC type to compute",
                      help="1) FSC overall - between two half-maps\n"
                           "2) FSC model-map - between atomic model and refined map\n"
                           "3) FSC work - model refined against half-map 1, "
                           "compared to half-map 1\n"
                           "FSC free - model refined against half-map 1, "
                           "compared to half-map 2")

        form.addParam('half1', params.PointerParam, pointerClass='Volume',
                      condition="fscType!=%d" % FSC_TYPE_MODEL_MAP,
                      label="Input half map 1",
                      important=True, allowsNull=True)
        form.addParam('half2', params.PointerParam, pointerClass='Volume',
                      condition="fscType!=%d" % FSC_TYPE_MODEL_MAP,
                      label="Input half map 2",
                      important=True, allowsNull=True)

        form.addParam('map', params.PointerParam, pointerClass='Volume',
                      condition="fscType==%d" % FSC_TYPE_MODEL_MAP,
                      label="Final map",
                      important=True, allowsNull=True)

        form.addParam('model', params.PointerParam, pointerClass='AtomStruct',
                      condition="fscType==%d" % FSC_TYPE_MODEL_MAP,
                      label="Final atomic model",
                      important=True, allowsNull=True)

        form.addParam('model_half1', params.PointerParam, pointerClass='AtomStruct',
                      condition="fscType==%d" % FSC_TYPE_WORK_FREE,
                      label="Atomic model refined against half-map 1",
                      important=True, allowsNull=True)

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.calculateFSCStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, *args):
        """ Create links and convert pdb to mrc. """
        if self._getFSCType() == FSC_TYPE_OVERALL:
            self._createInputLink(self.half1.get().getFileName(), 'half1')
            self._createInputLink(self.half2.get().getFileName(), 'half2')

        elif self._getFSCType() == FSC_TYPE_MODEL_MAP:
            self._createInputLink(self.map.get().getFileName(), 'map')

            box = self.map.get().getDim()[0]
            apix = self.map.get().getSamplingRate()

            self._model2Map(self.model.get().getFileName(),
                            self._getFileName('model'), box, apix)

        else:
            self._createInputLink(self.half1.get().getFileName(), 'half1')
            self._createInputLink(self.half2.get().getFileName(), 'half2')

            box = self.half1.get().getDim()[0]
            apix = self.half1.get().getSamplingRate()

            self._model2Map(self.model_half1.get().getFileName(),
                            self._getFileName('model_half1'), box, apix)

    def calculateFSCStep(self):
        """ Run relion_image_handler. """
        if self._getFSCType() == FSC_TYPE_OVERALL:
            angpix = self.half1.get().getSamplingRate()
            params = self._getParams(self._getFileName('half1'),
                                     self._getFileName('half2'),
                                     self._getFileName('fsc_overall'), angpix)
        elif self._getFSCType() == FSC_TYPE_MODEL_MAP:
            angpix = self.map.get().getSamplingRate()
            params = self._getParams(self._getFileName('model'),
                                     self._getFileName('map'),
                                     self._getFileName('fsc_model-map'), angpix)
        else:
            # FSC work
            angpix = self.half1.get().getSamplingRate()
            params = self._getParams(self._getFileName('model_half1'),
                                     self._getFileName('half1'),
                                     self._getFileName('fsc_work'), angpix)
            self._runProgram('relion_image_handler', params)

            # FSC free
            params = self._getParams(self._getFileName('model_half1'),
                                     self._getFileName('half2'),
                                     self._getFileName('fsc_free'), angpix)

        self._runProgram('relion_image_handler', params)

    def createOutputStep(self):
        if self._getFSCType() == FSC_TYPE_OVERALL:
            fsc = self._getFSC(self._getFileName('fsc_overall'),
                               'FSC overall')
            self._defineOutputs(outputFSC=fsc)

        elif self._getFSCType() == FSC_TYPE_MODEL_MAP:
            fsc = self._getFSC(self._getFileName('fsc_model-map'),
                               'FSC model-map')
            self._defineOutputs(outputFSC=fsc)

        else:
            fscSet = self._createSetOfFSCs()
            fscw = self._getFSC(self._getFileName('fsc_work'),
                                'FSC work')
            fscf = self._getFSC(self._getFileName('fsc_free'),
                                'FSC free')
            fscSet.append(fscw)
            fscSet.append(fscf)
            self._defineOutputs(outputSetOfFSCs=fscSet)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []

        from pwem import Domain
        try:
            _ = Domain.importFromPlugin('eman2', doRaise=True)
        except:
            errors.append("EMAN2 is required to convert pdb to mrc map")
        return errors

    def _citations(self):
        return ['Amunts2014']

    def _summary(self):
        summary = []

        if hasattr(self, "outputFSC") or hasattr(self, "outputSetOfFSCs"):
            summary.append("FSC calculation ran with relion_image_handler")

        return summary

    # -------------------------- UTILS functions ------------------------------
    def _createInputLink(self, fn, link):
        if pwutils.getExt(fn) == ".mrc":
            pwutils.createAbsLink(os.path.abspath(fn), self._getFileName(link))
        else:
            ih = ImageHandler()
            ih.convert(fn, self._getFileName(link))

    def _getFSCType(self):
        return self.fscType.get()

    def _model2Map(self, model, map, box, apix):
        """ Convert model pdb to mrc map. """
        from pwem import Domain
        eman2 = Domain.importFromPlugin('eman2',
                                        errorMsg='EMAN2 is required to convert pdb to mrc map',
                                        doRaise=True)
        from pyworkflow.utils.process import runJob

        args = "%s %s --box %d --apix %0.3f" % (model, map, box, apix)
        runJob(self._log, eman2.Plugin.getProgram('e2pdb2mrc.py'), args,
               env=eman2.Plugin.getEnviron())

    def _getFSC(self, fn, label):
        fsc = FSC(objLabel=label)
        table = Table(fileName=fn, tableName='fsc')
        resolution_inv = table.getColumnValues('rlnResolution')
        frc = table.getColumnValues('rlnFourierShellCorrelation')
        fsc.setData(resolution_inv, frc)
        return fsc

    def _getParams(self, fn1, fn2, output, angpix):
        params = ' '.join([
            '--i %s ' % fn1,
            '--fsc %s' % fn2,
            '--angpix %0.3f' % angpix,
            '> %s' % output
        ])
        return params
