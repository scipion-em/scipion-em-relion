# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco     (josue.gomez-blanco@mcgill.ca) [1]
# *              J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Department of Anatomy and Cell Biology, McGill University
# * [2] SciLifeLab, Stockholm University
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
from enum import Enum
from emtable import Table

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume

import relion.convert as convert
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputVolume = Volume


class ProtRelionPostprocess(ProtAnalysis3D, ProtRelionBase):
    """
    Performs post-processing of cryo-EM reconstructions using Relion in order
    to improve map interpretability through masking, sharpening, filtering,
    and resolution estimation.

    AI Generated:

    Relion Postprocess (ProtRelionPostprocess) - User Manual
        Overview

        The Relion Postprocess protocol refines the final appearance and
        interpretability of reconstructed cryo-EM density maps after
        three-dimensional refinement. Its purpose is to produce a biologically
        meaningful final map by combining half-map information with solvent
        masking, Fourier Shell Correlation analysis, map sharpening, and
        detector correction procedures. This stage is often one of the final
        steps before visualization, atomic modeling, validation, or deposition.

        In practical cryo-EM workflows, reconstructed maps frequently contain
        residual noise, dampened high-resolution signal, or contrast imbalance.
        Post-processing improves the visibility of structural features such as
        secondary structure elements, side chains, ligand densities, or flexible
        regions. The protocol also estimates the final resolution using
        gold-standard FSC approaches, helping users assess the reliability of
        structural details.

        Inputs and Reconstruction Sources

        The protocol supports different starting points depending on the
        workflow used previously. Users may continue directly from Relion
        refinement protocols, including standard refinement and multi-body
        refinement approaches, or provide external half maps manually.

        In most biological workflows, the preferred input consists of two
        independent half maps generated through gold-standard refinement
        procedures. These half maps allow unbiased resolution estimation and
        reliable sharpening. Users may either provide a reconstructed volume
        already associated with half maps or directly supply the two half maps
        separately.

        When working with multi-body refinement results, individual flexible
        bodies may be processed independently. This is especially useful for
        macromolecular assemblies that contain large conformational variability,
        such as ribosomes, membrane complexes, or multi-domain proteins. In
        these situations, separate post-processing of each body often improves
        local interpretability and produces cleaner density maps.

        Solvent Masking and Biological Interpretation

        Solvent masking is one of the most influential aspects of post-processing
        because it determines which regions contribute to FSC estimation and
        sharpening. The protocol expects a soft solvent mask in which the
        molecular region is represented by high values while solvent regions are
        suppressed.

        Biologically meaningful masks should include the complete molecular
        envelope while avoiding excessive solvent or disconnected noise regions.
        Masks that are too tight may artificially inflate resolution estimates,
        whereas masks that are too loose may reduce sharpening effectiveness and
        obscure high-resolution detail.

        Soft mask edges are generally preferred because abrupt transitions may
        introduce Fourier artifacts. In practice, smooth masks are particularly
        important for flexible assemblies, membrane proteins, or elongated
        particles where density boundaries are not sharply defined.

        Resolution Estimation and FSC Weighting

        The protocol evaluates map quality using Fourier Shell Correlation
        between independent half maps. The resulting FSC curve provides an
        estimate of the spatial frequency range supported by reproducible signal.
        This information is then used to guide map filtering and sharpening.

        By default, FSC weighting is applied automatically. This strategy helps
        suppress noisy frequencies while preserving reliable structural detail.
        For most biological projects, FSC weighting provides a balanced final
        map suitable for interpretation and model building.

        In some cases, however, users may choose to bypass FSC-based weighting
        and instead apply an ad-hoc low-pass filter. This option can be useful
        when maps exhibit strong local-resolution variability, where certain
        regions are significantly better resolved than the global FSC estimate
        suggests. Biological caution is important in these situations because
        aggressive filtering or insufficient filtering may lead to either loss
        of signal or over-interpretation of noise.

        B-Factor Sharpening

        Cryo-EM maps often suffer from attenuation of high-frequency information.
        B-factor sharpening compensates for this effect and increases the
        visibility of fine structural features. The protocol supports both
        automatic B-factor estimation and user-defined sharpening values.

        Automatic estimation is generally recommended for routine workflows
        because it derives sharpening parameters directly from the Guinier
        behavior of the reconstruction. This approach often provides balanced
        enhancement without excessive amplification of noise.

        Manual B-factor specification may be useful for experienced users who
        wish to optimize map appearance for specific biological questions.
        Strong sharpening can improve visibility of side chains and secondary
        structure elements, but excessive sharpening may generate misleading
        high-frequency artifacts. Careful visual inspection and validation are
        therefore essential.

        Detector MTF Correction

        The protocol supports modulation transfer function correction using
        detector-specific MTF curves. MTF correction compensates for signal
        attenuation introduced by the detector and may improve recovery of
        high-resolution information.

        This option becomes particularly relevant in high-resolution cryo-EM
        projects where subtle structural features are important. Accurate MTF
        correction requires the appropriate detector curve and an accurate
        detector pixel size calibration. When properly applied, the resulting
        sharpened maps may display improved contrast and interpretability.

        Pixel Size Calibration

        The protocol allows the use of calibrated pixel sizes that differ from
        the original refinement values. This feature is biologically important
        when pixel size recalibration has been performed using atomic models,
        diffraction standards, or independent validation procedures.

        Accurate pixel calibration directly influences reported resolution,
        spatial measurements, and downstream structural interpretation.
        Therefore, users should ensure consistency between reconstruction,
        refinement, and post-processing parameters whenever recalibration is
        applied.

        Outputs and Interpretation

        The primary output is a sharpened and filtered density map optimized for
        biological interpretation. The protocol also provides quantitative
        information about final resolution and the sharpening parameters used
        during processing.

        The resulting map is typically suitable for visualization, segmentation,
        atomic model fitting, flexible fitting, or deposition into public
        databases. However, users should remember that map appearance depends
        strongly on masking and sharpening choices. Apparent high-resolution
        features should always be interpreted together with FSC validation and
        biological consistency.

        Practical Recommendations

        For most workflows, automatic B-factor estimation together with FSC
        weighting provides a reliable starting point. A carefully prepared soft
        solvent mask usually has the largest impact on obtaining stable and
        biologically meaningful results.

        When processing flexible assemblies or multi-body refinements, users
        should evaluate each region independently because local flexibility may
        strongly influence sharpening behavior and apparent resolution.

        Excessive sharpening should be avoided, especially when density maps are
        noisy or contain heterogeneous regions. Visual inspection alongside FSC
        validation remains one of the most important quality-control steps.

        Final Perspective

        Post-processing is not simply a cosmetic operation but a critical stage
        in cryo-EM structure determination that directly affects biological
        interpretation. Appropriate masking, careful sharpening, accurate
        resolution estimation, and conservative interpretation together produce
        density maps that more faithfully represent the underlying molecular
        structure.
    """
    _label = 'post-processing'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _getInputPath(self, *paths):
        return self._getPath('input', *paths)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'finalVolume': self._getInputPath("relion_class001.mrc"),
            'half1': self._getInputPath("relion_half1_class001_unfil.mrc"),
            'half2': self._getInputPath("relion_half2_class001_unfil.mrc"),
            'mask': self._getInputPath("input_mask.mrc"),
            'outputVolume': self._getExtraPath('postprocess.mrc')
        }
        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('relionInput', params.BooleanParam,
                      default=True,
                      label="Start from Relion refinement?",
                      help="Set to Yes if you wish to use as input "
                           "a Relion protocol. Otherwise set it to No")
        form.addParam('inputType', params.EnumParam, default=0,
                      display=params.EnumParam.DISPLAY_HLIST,
                      choices=['Volume with half-maps',
                               'Individual half-maps'],
                      label='Input type:',
                      condition='not relionInput')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Input volume with half-maps",
                      allowsNull=True,
                      condition="not relionInput and inputType==0")
        form.addParam('inputHalf1', params.PointerParam,
                      pointerClass='Volume',
                      label="Input half map 1",
                      important=True, allowsNull=True,
                      condition="not relionInput and inputType==1",
                      help='You might want to provide input half maps manually, '
                           'in case you did not use 3D auto-refine or multi-body protocol '
                           'that generates them automatically.')
        form.addParam('inputHalf2', params.PointerParam, pointerClass='Volume',
                      label="Input half map 2",
                      important=True, allowsNull=True,
                      condition="not relionInput and inputType==1",
                      help='You might want to provide input half maps manually, '
                           'in case you did not use 3D auto-refine or multi-body protocol '
                           'that generates them automatically.')

        form.addParam('protRefine', params.PointerParam,
                      pointerClass="ProtRefine3D, ProtRelionMultiBody",
                      condition="relionInput",
                      label='Select a previous refinement protocol',
                      help='Select any previous refinement protocol to get the '
                           '3D half maps. Note that it is recommended that the '
                           'refinement protocol uses a gold-standard method.')
        form.addParam('bodyNum', params.IntParam, default=1,
                      condition="relionInput",
                      label="Which body to process?",
                      help="Only relevant if input protocol is 3D multi-body.")
        form.addParam('solventMask', params.PointerParam,
                      pointerClass="VolumeMask",
                      label='Solvent mask',
                      help="Provide a soft mask where the protein is white "
                           "(1) and the solvent is black (0). Often, the "
                           "softer the mask the higher resolution estimates "
                           "you will get. A soft edge of 5-10 pixels is often "
                           "a good edge width.")
        form.addParam('calibratedPixelSize', params.FloatParam, default=0,
                      label='Calibrated pixel size (A)',
                      help="Provide the final, calibrated pixel size in "
                           "Angstroms. If 0, the input pixel size will be used. "
                           "This value may be different from the pixel-size "
                           "used thus far, e.g. when you have recalibrated "
                           "the pixel size using the fit to a PDB model. "
                           "The X-axis of the output FSC plot will use this "
                           "calibrated value.")

        form.addSection(label='Sharpening')
        group = form.addGroup('MTF')
        group.addParam('mtf', params.FileParam,
                       label='MTF of the detector',
                       help='User-provided STAR-file with the MTF-curve '
                            'of the detector. Use the wizard to load one '
                            'of the predefined ones provided at:\n'
                            '- [[https://www3.mrc-lmb.cam.ac.uk/relion/index.php/'
                            'FAQs#Where_can_I_find_MTF_curves_for_typical_detectors.3F]'
                            '[Relion\'s Wiki FAQs]]\n'
                            ' - [[https://www.gatan.com/techniques/cryo-em#MTF][Gatan\'s website]]\n\n'
                            'Relion param: *--mtf*')
        group.addParam('origPixelSize', params.FloatParam,
                       default=-1.0,
                       label='Original detector pixel size (A)',
                       help='This is the original pixel size (in Angstroms)'
                            ' in the raw (non-super-resolution!) micrographs')

        form.addParam('doAutoBfactor', params.BooleanParam, default=True,
                      label='Estimate B-factor automatically?',
                      help='If set to Yes, then the program will use the '
                           'automated procedure described by Rosenthal and '
                           'Henderson (2003, JMB) to estimate an overall '
                           'B-factor for your map, and sharpen it accordingly.')
        line = form.addLine('B-factor resolution (A): ',
                            condition='doAutoBfactor',
                            help='There are the frequency (in Angstroms), '
                                 'lowest and highest, that will be included in '
                                 'the linear fit of the Guinier plot as '
                                 'described in Rosenthal and Henderson '
                                 '(2003, JMB).')
        line.addParam('bfactorLowRes', params.FloatParam,
                      default=10.0, label='low')
        line.addParam('bfactorHighRes', params.FloatParam,
                      default=0.0, label='high')
        form.addParam('bfactor', params.FloatParam, default=-350,
                      condition='not doAutoBfactor',
                      label='Provide B-factor:',
                      help='User-provided B-factor (in A^2) for map '
                           'sharpening, e.g. -400. Use negative values for '
                           'sharpening. Be careful: if you over-sharpen\n'
                           'your map, you may end up interpreting noise for '
                           'signal!\n'
                           'Relion param: *--adhoc_bfac*')

        form.addSection(label='Filtering')
        form.addParam('skipFscWeighting', params.BooleanParam, default=False,
                      label='Skip FSC-weighting for sharpening?',
                      help='If set to No (the default), then the output map '
                           'will be low-pass filtered according to the '
                           'mask-corrected, gold-standard FSC-curve. '
                           'Sometimes, it is also useful to provide an ad-hoc '
                           'low-pass filter (option below), as due to local '
                           'resolution variations some parts of the map may '
                           'be better and other parts may be worse than the '
                           'overall resolution as measured by the FSC. In '
                           'such  cases, set this option to Yes and provide '
                           'an ad-hoc filter as described below.')
        form.addParam('lowRes', params.FloatParam, default=5,
                      condition='skipFscWeighting',
                      label='Ad-hoc low-pass filter (A):',
                      help='This option allows one to low-pass filter the map '
                           'at a user-provided frequency (in Angstroms). When '
                           'using a resolution that is higher than the '
                           'gold-standard FSC-reported resolution, take care '
                           'not to interpret noise in the map for signal.')
        form.addParam('filterEdgeWidth', params.IntParam, default=2,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Low-pass filter edge width:',
                      help='Width of the raised cosine on the low-pass filter '
                           'edge (in resolution shells)\n'
                           'Relion param: *--filter_edge_width*')
        form.addParam('randomizeAtFsc', params.FloatParam, default=0.8,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Randomize phases threshold',
                      help='Randomize phases from the resolution where FSC '
                           'drops below this value\n'
                           'Relion param: *--randomize_at_fsc*')
        form.addParam('forceMask', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Force mask?',
                      help='Use the mask even when the masked resolution '
                           'is worse than the unmasked resolution.')

        form.addParallelSection(threads=0, mpi=0)

    # -------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        if self.relionInput:
            objsId = [self.protRefine.get().getObjId()]
        else:
            if self.inputType.get() == 0:  # volume with half-maps
                objsId = [self.inputVolume.get().getObjId()]
            else:
                objsId = [self.inputHalf1.get().getObjId(),
                          self.inputHalf2.get().getObjId()]
        self._createFilenameTemplates()
        self._defineParamDict()
        self._insertFunctionStep(self.convertInputStep, objsId,
                                 needsGPU=False)
        self._insertFunctionStep(self.postProcessStep, self.paramDict,
                                 needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, objsId):
        pwutils.makePath(self._getInputPath())
        ih = ImageHandler()

        if self.relionInput:
            protRef = self.protRefine.get()

            if self._isInputMbody():
                for i, vol in enumerate(protRef.outputVolumes):
                    if i == self.bodyNum.get()-1:
                        outVol = vol
                        self.info("Using multi-body input: %s" % outVol.getHalfMaps())
                        break
            else:  # ProtRefine3D
                outVol = protRef.outputVolume

            newDim = outVol.getXDim()
            newPix = outVol.getSamplingRate()
            vols = outVol.getHalfMaps(asList=True)

            for vol, key in zip(vols, ['half1', 'half2']):
                ih.convert(vol, self._getFileName(key))
        else:
            if self.inputType.get() == 0:  # volume with half-maps
                [half1, half2] = self.inputVolume.get().getHalfMaps(asList=True)
                newDim = self.inputVolume.get().getXDim()
                newPix = self.inputVolume.get().getSamplingRate()
            else:
                half1 = self.inputHalf1.get()
                half2 = self.inputHalf2.get()
                newDim = half1.getXDim()
                newPix = half1.getSamplingRate()
                half1 = half1.getFileName()
                half2 = half2.getFileName()

            ih.convert(half1, self._getFileName('half1'))
            ih.convert(half2, self._getFileName('half2'))

        convert.convertMask(self.solventMask.get(),
                            self._getFileName('mask'), newPix, newDim)

    def postProcessStep(self, paramDict):
        params = ' '.join(['%s %s' % (k, str(v))
                           for k, v in paramDict.items()])
        self._runProgram('relion_postprocess', params)

    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getFileName('outputVolume'))
        if not self.relionInput:
            if self.inputType.get() == 0:
                vol = self.inputVolume
            else:
                vol = self.inputHalf1
        elif self._isInputMbody():
            vol = self.protRefine.get().outputVolumes
        else:
            vol = self.protRefine.get().outputVolume

        volume.setSamplingRate(self._getOutputPixelSize())
        self._defineOutputs(**{outputs.outputVolume.name: volume})
        self._defineSourceRelation(vol, volume)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        mtfFile = self.mtf.get()

        if mtfFile and not os.path.exists(mtfFile):
            errors.append("Missing MTF-file '%s'" % mtfFile)

        if (self.inputType.get() == 0 and
                self.inputVolume.hasValue() and
                not self.inputVolume.get().hasHalfMaps()):
            errors.append("Input volume is missing half-maps")

        return errors

    def _citations(self):
        return ['Chen2013']

    def _summary(self):
        summary = []
        postStarFn = self._getExtraPath("postprocess.star")
        if os.path.exists(postStarFn):
            table = Table(fileName=postStarFn, tableName='general')
            row = table[0]
            summary.append("Final resolution: *%0.2f A*" %
                           float(row.rlnFinalResolution))
            summary.append("B-factor: *%0.2f A\u00B2*" %
                           float(row.rlnBfactorUsedForSharpening))

        return summary

    # -------------------------- UTILS functions ------------------------------
    def _defineParamDict(self):
        """ Define all parameters to run relion_postprocess"""
        inputFn = self._getFileName('half1')

        self.paramDict = {'--i': inputFn,
                          '--o': self._getExtraPath('postprocess'),
                          '--angpix': self._getOutputPixelSize(),
                          '--filter_edge_width': self.filterEdgeWidth.get(),
                          '--randomize_at_fsc': self.randomizeAtFsc.get(),
                          '--mask': self._getFileName('mask')
                          }

        mtfFile = self.mtf.get()
        if mtfFile:
            self.paramDict['--mtf'] = mtfFile

        if self.doAutoBfactor:
            self.paramDict['--auto_bfac'] = ''
            self.paramDict['--autob_lowres'] = self.bfactorLowRes.get()
            self.paramDict['--autob_highres'] = self.bfactorHighRes.get()
        else:
            self.paramDict['--adhoc_bfac'] = self.bfactor.get()

        if self.skipFscWeighting:
            self.paramDict['--skip_fsc_weighting'] = ''
            self.paramDict['--low_pass'] = self.lowRes.get()

        if self.origPixelSize.get() != -1.0:
            self.paramDict['--mtf_angpix'] = self.origPixelSize.get()

        if self.forceMask:
            self.paramDict['--force_mask'] = ''

    def _isInputMbody(self):
        return self.protRefine.get().getClassName() == "ProtRelionMultiBody"

    def _getOutputPixelSize(self):
        """ Return the output pixel size, using the calibrated
        pixel size if non zero, or the input one. """
        if not self.relionInput:
            if self.inputType.get() == 0:
                volume = self.inputVolume.get()
            else:
                volume = self.inputHalf1.get()
        elif self._isInputMbody():
            volume = self.protRefine.get().outputVolumes
        else:
            volume = self.protRefine.get().outputVolume
        cps = self.calibratedPixelSize.get()
        return cps if cps > 0 else volume.getSamplingRate()
