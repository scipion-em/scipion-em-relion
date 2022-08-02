# ******************************************************************************
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
# ******************************************************************************

import os

import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
from pwem.protocols import ProtAlignMovies
from pwem.objects import MovieAlignment
from pyworkflow.protocol import STEPS_PARALLEL


class ProtRelionCompressMovies(ProtAlignMovies):
    """
    Using *relion_convert_to_tiff* to compress a set of movies.
    """
    _label = 'compress movies'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        self.isEER = False
        self.stepsExecutionMode = STEPS_PARALLEL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif', '.eer'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        form.addParam('inputGainProt', params.PointerParam,
                      allowsNull=True,
                      pointerClass='ProtRelionCompressEstimateGain',
                      label='Input gain estimation protocol (optional)',
                      help='Provide an estimate gain reference protocol '
                           'from where the gain file will be taken.')

        group = form.addGroup("TIFF Options")
        group.addParam('compression', params.EnumParam, default=1,
                       choices=['none', 'auto', 'zip', 'lzw'],
                       label='Compression type')
        group.addParam('deflateLevel', params.IntParam, default=6,
                       label='Deflate level',
                       help="deflate level. 1 (fast) "
                            "to 9 (slowest but best compression)")

        form.addSection("EER")
        form.addParam('eerGroup', params.IntParam, default=32,
                      label='EER fractionation',
                      help="The number of hardware frames to group into one "
                           "fraction. This option is relevant only for Falcon4 "
                           "movies in the EER format. Falcon 4 operates at "
                           "248 frames/s.\nFractionate such that each fraction "
                           "has about 0.5 to 1.25 e/A2.")
        form.addParam('eerSampling', params.EnumParam, default=0,
                      choices=['1', '2'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='EER upsampling',
                      help="EER upsampling (1 = 4K or 2 = 8K). 8K rendering is not "
                           "recommended by Relion. See "
                           "https://relion.readthedocs.io/en/latest/Reference/MovieCompression.html")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions -------------------------------
    def _convertInputStep(self):
        self.info("Relion version:")
        self.runJob("relion_convert_to_tiff --version", "", numberOfMpi=1)
        # Create a local link to the input gain file if necessary
        if self.inputGainProt.get():
            inputGainProt = self.inputGainProt.get()
            inputGainProt._createFilenameTemplates()

            tmpGain = self._getTmpPath('gain_estimate.bin')
            pwutils.createAbsLink(os.path.abspath(inputGainProt._getFileName("output_gain")),
                                  tmpGain)
            pwutils.createAbsLink(os.path.abspath(inputGainProt._getFileName("output_gain_extra")),
                                  tmpGain.replace('.bin', '_reliablity.bin'))

        ProtAlignMovies._convertInputStep(self)

    def _processMovie(self, movie):
        fn = movie.getAttributeValue('_originalFileName', movie.getFileName())
        baseName = os.path.basename(fn)
        compression = self.getEnumText('compression')
        pwutils.createAbsLink(os.path.abspath(fn), self._getTmpPath(baseName))
        args = " --i %s --o ../extra/" % baseName
        args += " --compression %s" % compression
        # TODO: Check if deflateLevel is only valid for zip (deflate)
        if compression == 'zip':  # deflate
            args += " --deflate_level %d" % self.deflateLevel

        if self.inputGainProt.get():
            args += " --gain gain_estimate.bin"

        if self.isEER:
            args += " --eer_grouping %d" % self.eerGroup
            args += " --eer_upsampling %d" % (self.eerSampling.get() + 1)

        self.runJob('relion_convert_to_tiff', args, cwd=self._getTmpPath())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["Movies compressed by relion_convert_to_tiff, "
                   "compression type: %s" % self.getEnumText('compression')]

        return summary

    def _citations(self):
        return ['Zivanov2019']

    def _validate(self):
        errors = []
        firstMovie = self.inputMovies.get().getFirstItem()
        self.isEER = pwutils.getExt(firstMovie.getFileName()) == ".eer"

        errors.extend(ProtAlignMovies._validate(self))

        return errors

    def _warnings(self):
        warnings = []
        if self.isEER:
            warnings.append("Note that after compression into TIFF the "
                            "original EER gain reference should be inverted.")

        return warnings

    # ------------------------ Extra BASE functions ---------------------------
    def _getRelPath(self, baseName, refPath):
        return os.path.relpath(self._getExtraPath(baseName), refPath)

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _createOutputMovies(self):
        return True

    def _createOutputMicrographs(self):
        return False

    def _createOutputWeightedMicrographs(self):
        return False

    def _savePsSum(self):
        return self.getAttributeValue('savePSsum', False)

    def _createOutputMovie(self, movie):
        n = movie.getNumberOfFrames()
        newMovie = movie.clone()
        movieName = pwutils.replaceBaseExt(movie.getFileName(), 'tif')
        newMovie.setFileName(self._getExtraPath(movieName))
        newMovie.setAlignment(MovieAlignment(first=1, last=n,
                                             xshifts=[0]*n, yshifts=[0]*n))
        return newMovie

    def _updateOutputSet(self, outputName, outputSet,
                         state=pwobj.Set.STREAM_OPEN):
        """ Redefine this method to update gain file. """
        first = getattr(self, '_firstUpdate', True)
        if first and self.inputGainProt.get():
            outputSet.setGain(os.path.abspath(self._getExtraPath("gain-reference.mrc")))

        ProtAlignMovies._updateOutputSet(self, outputName, outputSet, state=state)
