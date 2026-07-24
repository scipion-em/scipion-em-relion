# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
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

from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD

import relion
import relion.convert
from .protocol_postprocess import ProtRelionPostprocess


class ProtRelionLocalRes(ProtRelionPostprocess):
    """
    Estimates local resolution variations within cryo-EM density maps using
    RELION-based postprocessing procedures. The protocol evaluates how map
    quality changes across different structural regions and generates locally
    filtered and sharpened reconstructions that help distinguish well-resolved
    domains from flexible or poorly ordered areas.

    AI Generated:

    Local Resolution Estimation (ProtRelionLocalRes) — User Manual
        Overview

        The Local Resolution protocol estimates spatial variations in map
        resolution across a cryo-EM reconstruction using RELION local
        resolution analysis tools. Instead of assigning a single global
        resolution value to the entire reconstruction, this protocol evaluates
        how structural quality changes from one region to another. This is
        particularly important for biological assemblies that contain flexible
        domains, mobile subunits, partially occupied regions, or compositional
        heterogeneity.

        In practical cryo-EM analysis, global Fourier shell correlation values
        often hide important local differences in map quality. Many complexes
        contain rigid cores with excellent structural detail together with
        peripheral regions that are substantially more flexible. Local
        resolution estimation helps users identify which regions are suitable
        for atomic interpretation and which regions should be interpreted with
        caution.

        Inputs and General Workflow

        The protocol requires two independently refined half maps originating
        from a previous gold-standard refinement procedure. These half maps are
        essential because the local resolution estimation relies on comparing
        independent signal contributions while minimizing overfitting effects.
        For best results, the input refinement should already be well converged
        and properly masked.

        During processing, the protocol evaluates local frequency content
        throughout the reconstruction by moving a soft spherical region across
        the map and estimating the local signal-to-noise characteristics. The
        final outputs include a local resolution map together with a locally
        filtered reconstruction in which each region is filtered according to
        its estimated local quality.

        Biological Interpretation of Local Resolution

        Local resolution maps provide direct biological insight into molecular
        flexibility and conformational variability. Regions with higher local
        resolution often correspond to rigid structural cores, stable secondary
        structure elements, tightly interacting interfaces, or highly occupied
        domains. Lower local resolution regions frequently indicate molecular
        motion, partial occupancy, conformational heterogeneity, or intrinsic
        disorder.

        For membrane proteins, local resolution differences commonly appear
        between the transmembrane core and flexible extracellular or cytosolic
        regions. In ribosomes or large macromolecular assemblies, peripheral
        domains and dynamic interaction partners often display lower local
        resolution than the central scaffold.

        From a biological perspective, local resolution estimation is therefore
        not only a technical quality-control procedure but also an indirect
        indicator of structural dynamics within the sample.

        Solvent Masking and Region Selection

        The protocol optionally accepts a solvent mask that defines the region
        of interest for analysis and visualization. Although the local
        resolution calculation itself is not strictly driven by the mask in the
        same way as refinement masking, the mask strongly influences histogram
        interpretation and helps isolate meaningful molecular regions from
        solvent background.

        In biological practice, masks should include all structurally relevant
        domains while avoiding excessive solvent regions. Poor masks may
        artificially distort local resolution statistics or complicate visual
        interpretation. For highly flexible complexes, broader masks are often
        preferable because overly restrictive masking may exclude biologically
        relevant motions.

        Sharpening and B-Factor Considerations

        The protocol supports application of a sharpening B-factor during local
        filtering. Negative B-factors enhance high-frequency information and
        can significantly improve map interpretability. However, excessive
        sharpening may amplify noise and create misleading structural features.

        In biological interpretation, sharpening should always be evaluated
        visually together with prior biochemical and structural knowledge.
        Strong sharpening may produce apparent side-chain features in regions
        that are not genuinely resolved. Conservative sharpening is generally
        safer for flexible regions or medium-resolution reconstructions.

        Pixel Size Calibration

        The protocol allows users to provide a calibrated pixel size that may
        differ from the original acquisition value. This becomes important when
        the microscope magnification has been refined using atomic models,
        diffraction standards, or post-acquisition calibration procedures.
        Accurate pixel size calibration improves the interpretation of
        resolution values and ensures consistency between reconstruction and
        atomic modeling workflows.

        Detector MTF Correction

        Advanced users may include a modulation transfer function correction
        for the detector. MTF correction compensates for detector-dependent
        attenuation of high-frequency information and may improve sharpening
        behavior and local resolution estimation accuracy.

        In most routine workflows, standard detector calibrations are
        sufficient. However, high-resolution projects aiming for atomic detail
        may benefit from carefully characterized detector MTF curves,
        particularly when comparing datasets acquired under different imaging
        conditions.

        Advanced Local Resolution Parameters

        The protocol exposes several advanced parameters controlling local
        sampling behavior, spherical mask dimensions, edge smoothing, phase
        randomization limits, and minimum allowed resolutions. These settings
        mainly influence the balance between spatial sensitivity and numerical
        stability.

        Smaller sampling regions may reveal fine local variations but can
        increase noise sensitivity. Larger regions provide smoother and more
        stable estimates but may obscure sharp transitions between rigid and
        flexible domains. For most biological applications, default parameters
        provide reliable results and should only be modified when specific map
        characteristics justify additional optimization.

        Outputs and Their Interpretation

        The protocol produces a local resolution map together with a locally
        filtered reconstruction. The local resolution map can be visualized as
        a colored overlay in molecular visualization software, allowing direct
        inspection of structural quality across the reconstruction.

        The locally filtered map is often more biologically interpretable than
        a globally filtered reconstruction because each region is filtered
        according to its own estimated quality. Well-resolved domains retain
        high-frequency detail, while poorly resolved regions remain smoother
        and less noisy.

        These outputs are particularly valuable during model building,
        validation, figure preparation, and interpretation of flexible
        assemblies.

        Practical Recommendations

        For most cryo-EM workflows, local resolution estimation should be
        performed after obtaining the final refined reconstruction and before
        extensive atomic interpretation. Users should visually compare the
        local resolution distribution with known flexible regions, biochemical
        expectations, and conformational variability observed during
        classification.

        Conservative sharpening is generally recommended during initial
        analysis. If certain regions appear artificially noisy or fragmented,
        reducing sharpening strength often improves interpretability.
        Similarly, local resolution values should always be interpreted
        together with map appearance rather than treated as absolute indicators
        of atomic accuracy.

        Final Perspective

        Local resolution estimation is one of the most informative quality
        assessment procedures in modern cryo-EM analysis. Beyond providing a
        technical characterization of reconstruction quality, it offers direct
        insight into structural flexibility, conformational variability, and
        molecular stability. Careful interpretation of local resolution maps
        helps users distinguish reliable structural features from uncertain
        regions and supports more accurate biological conclusions.
    """
    _label = 'local resolution'
    _devStatus = PROD
    relionInput = True

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
                 'half1': self._getInputPath("relion_half1_class001_unfil.mrc"),
                 'half2': self._getInputPath("relion_half2_class001_unfil.mrc"),
                 'outputVolume': self._getExtraPath('relion_locres_filtered.mrc'),
                 'resolMap': self._getExtraPath('relion_locres.mrc'),
                 'solventMask': self._getExtraPath('input_solvent_mask.mrc')
                 }

        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('protRefine', params.PointerParam,
                      pointerClass="ProtRefine3D",
                      label='Select a previous refinement protocol',
                      help='Select any previous refinement protocol to get the '
                           '3D half maps. Note that it is recommended that the '
                           'refinement protocol uses a gold-standard method.')
        form.addParam('solventMask', params.PointerParam,
                      pointerClass='VolumeMask', allowsNull=True,
                      label='User-provided solvent mask',
                      help='Provide a mask with values between 0 and 1 '
                           'around all domains of the complex. ResMap uses '
                           'this mask for local resolution calculation. '
                           'RELION does NOT use this mask for calculation, '
                           'but makes a histogram of local resolution '
                           'within this mask.')
        form.addParam('calibratedPixelSize', params.FloatParam, default=0.,
                      label='Calibrated pixel size (A)',
                      help="Provide the final, calibrated pixel size in "
                           "Angstroms. If 0, the input pixel size will be used. "
                           "This value may be different from the pixel-size "
                           "used thus far, e.g. when you have recalibrated "
                           "the pixel size using the fit to a PDB model. "
                           "The X-axis of the output FSC plot will use this "
                           "calibrated value.")
        form.addParam('bfactor', params.FloatParam, default=-100.,
                      label='Provide B-factor:',
                      help='Probably, the overall B-factor as was '
                           'estimated in the postprocess is a useful '
                           'value for here. Use negative values for '
                           'sharpening. Be careful: if you over-sharpen '
                           'your map, you may end up interpreting '
                           'noise for signal!')

        group = form.addGroup('MTF')
        group.addParam('mtf', params.FileParam,
                       label='MTF of the detector',
                       help='User-provided STAR-file with the MTF-curve '
                            'of the detector.'
                            'Relion param: <--mtf>')
        group.addParam('origPixelSize', params.FloatParam,
                       default=-1.0,
                       label='Original detector pixel size (A)',
                       help='This is the original pixel size (in Angstroms)'
                            ' in the raw (non-super-resolution!) micrographs')

        form.addSection(label='LocalRes')
        form.addParam('Msg', params.LabelParam,
                      label='Select Advanced level if you want to adjust the '
                            'parameters')
        form.addParam('locResSamp', params.IntParam, default=25,
                      label='Sampling rate (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Sampling rate (in Angstroms) with which to '
                           'sample the local-resolution map')
        form.addParam('locResMaskRad', params.IntParam, default=-1,
                      label='Mask radius (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Radius (in A) of spherical mask for '
                           'local-resolution map (default = 0.5*sampling)')
        form.addParam('locResEdgeWidth', params.IntParam, default=-1,
                      label='Edge width (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Width of soft edge (in A) on masks for '
                           'local-resolution map (default = sampling)')
        form.addParam('locResRand', params.FloatParam, default=25.0,
                      label='Randomize phases from (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Randomize phases from this resolution (in A)')
        form.addParam('locResMin', params.IntParam, default=50,
                      label='Lowest res limit (A)',
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Lowest local resolution allowed (in A)')

        form.addParallelSection(threads=0, mpi=1)
    
    # ------------------------- STEPS functions -------------------------------
    def convertInputStep(self, protId):
        pwutils.makePath(self._getInputPath())

        protRef = self.protRefine.get()
        vol = protRef.outputVolume
        newDim = vol.getXDim()
        newPix = vol.getSamplingRate()
        half1, half2 = vol.getHalfMaps().split(',')
        ih = ImageHandler()
        ih.convert(half1, self._getFileName("half1"))
        ih.convert(half2, self._getFileName("half2"))

        if self.solventMask.hasValue():
            relion.convert.convertMask(self.solventMask.get(),
                                       self._getFileName('solventMask'),
                                       newPix, newDim)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        mtfFile = self.mtf.get()

        if mtfFile and not os.path.exists(mtfFile):
            errors.append("Missing MTF-file '%s'" % mtfFile)

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volume not ready yet.")
        else:
            output = self.outputVolume
            summary.append("%s: Output volume was locally filtered "
                           "and sharpened" % self.getObjectTag(output))
        return summary
    
    # -------------------------- UTILS functions ------------------------------
    def _defineParamDict(self):
        """ Define all parameters to run relion_postprocess"""
        volume = self.protRefine.get().outputVolume
        # It seems that in Relion3 now the input should be the map
        # filename and not the prefix as before
        inputFn = self._getFileName('half1')
        cps = self.calibratedPixelSize.get()
        angpix = cps if cps > 0 else volume.getSamplingRate()

        self.paramDict = {'--i': inputFn,
                          '--o': self._getExtraPath('relion'),
                          '--angpix': angpix,
                          '--adhoc_bfac': self.bfactor.get(),
                          '--locres': '',
                          # Expert options
                          '--locres_sampling': self.locResSamp.get(),
                          '--locres_maskrad': self.locResMaskRad.get(),
                          '--locres_edgwidth': self.locResEdgeWidth.get(),
                          '--locres_randomize_at': self.locResRand.get(),
                          '--locres_minres': self.locResMin.get()
                          }

        mtfFile = self.mtf.get()
        if mtfFile:
            self.paramDict['--mtf'] = mtfFile
        if self.origPixelSize.get() != -1.0:
            self.paramDict['--mtf_angpix'] = self.origPixelSize.get()

        if self.solventMask.hasValue():
            self.paramDict['--mask'] = self._getFileName('solventMask')
