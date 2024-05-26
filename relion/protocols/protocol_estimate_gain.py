# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
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
# ******************************************************************************

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.protocol.constants import STEPS_SERIAL
from pyworkflow.constants import PROD
from pwem.protocols import ProtProcessMovies
import pwem.emlib as emlib

import relion.convert as convert


class ProtRelionCompressEstimateGain(ProtProcessMovies):
    """
    Using *relion_estimate_gain* or *relion_compress_to_tiff* to estimate the gain reference from a set of movies.
    """
    _label = 'estimate gain reference'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtProcessMovies.__init__(self, **kwargs)
        self.isFloat32 = False
        self.isEER = False
        self.stepsExecutionMode = STEPS_SERIAL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff',
                               '.tif', '.eer', '.gain'] else 'mrc'

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {'input_star': self._getTmpPath('input_movies.star'),
                  'output_gain': self._getPath('gain_estimate.bin'),
                  'output_gain_extra': self._getPath('gain_estimate_reliablity.bin'),
                  'output_gain_mrc': self._getPath('gain_estimate.mrc'),
                  }
        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=pwutils.Message.LABEL_INPUT)
        form.addParam('inputMovies', params.PointerParam,
                      pointerClass='SetOfMovies',
                      important=True,
                      label=pwutils.Message.LABEL_INPUT_MOVS,
                      help='Select a set of movies to be used in the gain '
                           'estimation.')
        form.addParam('threshold', params.IntParam, default=50,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Pixel threshold',
                      help="*Only for 32-bit float MRC*\nNumber of success "
                           "needed to consider a pixel reliable. A pixel is "
                           "considered to be reliable when values which are "
                           "integer multiples of the current gain estimate "
                           "were observed at least --thresh times "
                           "(default 50) without being interrupted by mismatch.")
        form.addParam('maxFrames', params.IntParam, default=0,
                      label="Target number of frames to average",
                      help="Default 0 means use all.")
        form.addParam('random', params.BooleanParam, default=False,
                      label="Randomise input",
                      help="Randomise the order of input movies before "
                           "taking subset")
        form.addParam('eerSampling', params.EnumParam, default=1,
                      choices=['1x', '2x'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='EER upsampling',
                      help="EER upsampling (1x = physical 4K, "
                           "2x = super-resolution 8K)")

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- STEPS functions -------------------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep(self.convertInputStep,
                                 self.inputMovies.getObjId())
        self._insertFunctionStep(self.estimateGainStep)

    def convertInputStep(self, moviesId):
        self.info("Relion version:")
        self.runJob("relion_estimate_gain --version", "",
                    numberOfMpi=1)
        writer = convert.createWriter()
        writer.writeSetOfMovies(self.inputMovies.get(),
                                self._getFileName("input_star"))

    def estimateGainStep(self):
        args = [
            f"--i {self._getFileName('input_star')}",
            f"--j {self.numberOfThreads}",
        ]

        if self.isFloat32:
            program = 'relion_convert_to_tiff'
            args.extend([
                f"--thresh {self.threshold}",
                f"--o {self.getPath()}",
                "--estimate_gain"
            ])
        else:
            program = 'relion_estimate_gain'
            args.extend([
                f"--max_frames {self.maxFrames.get() if self.maxFrames.get() > 0 else -1}",
                f"--eer_upsampling {self.eerSampling.get()+1}",
                f"--o {self.getPath('gain_estimate.mrc')}"
            ])
            if self.random:
                args.append('--random')
            if self.isEER:
                args.append('--dont_invert')

        self.runJob(program, " ".join(args))

    def _stepsCheck(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('Gain estimation completed')
        return summary

    def _citations(self):
        return ['Zivanov2019']

    def _validate(self):
        errors = []
        inputMovies = self.inputMovies.get()
        firstMovie = inputMovies.getFirstItem()
        fn = firstMovie.getFileName()

        ih = emlib.image.ImageHandler()
        if "mrc" in pwutils.getExt(fn) and ih.getDataType(fn) == emlib.DT_FLOAT:
            self.isFloat32 = True
        if pwutils.getExt(fn) == ".eer":
            self.isEER = True

        return errors
