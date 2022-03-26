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

from pwem.emlib.image import ImageHandler
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD

import relion
import relion.convert
from .protocol_postprocess import ProtRelionPostprocess


class ProtRelionLocalRes(ProtRelionPostprocess):
    """ This protocol does local resolution estimation using Relion.

    This program basically performs a series of post-processing operations
    with a small soft, spherical mask that is moved over the entire map,
    while using phase-randomisation to estimate the convolution effects
    of that mask.
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
    
    # -------------------------- INSERT steps functions -----------------------

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
