# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# ******************************************************************************

import os
from itertools import izip
from math import ceil
import json

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol import STEPS_SERIAL

import relion
import relion.convert.metadata as md


class ProtRelionMotioncor(ProtAlignMovies):
    """
    Wrapper for the Relion's implementation of motioncor algorithm.
    """

    _label = 'motioncor'

    @classmethod
    def isDisabled(cls):
        return not relion.Plugin.isVersion3Active()

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):

        line = form.addLine('Frames for corrected SUM',
                             help='First and last frames to use in corrected '
                                  'average (starts counting at 1 and 0 as last '
                                  'means util the last frame in the movie). ')
        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')

        form.addParam('doDW', params.BooleanParam, default=False,
                      label='Do dose-weighting?',
                      help='If set to Yes, the averaged micrographs will be '
                           'dose-weighted. \n\n'
                           'NOTE: In Scipion the Voltage and and Dose '
                           'information is provided during import, so you '
                           'do not need to provide them anymore. ')

        form.addParam('saveNonDW', params.BooleanParam, default=False,
                      condition='doDW',
                      label='Save non-dose weighted as well?',
                      help='Aligned but non-dose weighted images are '
                           'sometimes useful in CTF estimation, although '
                           'there is no difference in most cases. Whichever '
                           'the choice, CTF refinement job is always done on '
                           'dose-weighted particles.')

        group = form.addGroup("Motion")

        group.addParam('bfactor', params.IntParam, default=150,
                       label='Bfactor',
                       help="The B-factor that will be applied to the "
                            "micrographs.")

        line = group.addLine('Number of patches',
                             help='Number of patches (in X and Y direction) to '
                                  'apply motion correction. \n')
        line.addParam('patchX', params.IntParam, default=1, label='X')
        line.addParam('patchY', params.IntParam, default=1, label='Y')

        group.addParam('groupFrames', params.IntParam, default=1,
                       label='Group frames',
                       help="Average together this many frames before "
                            "calculating the beam-induced shifts.")

        group.addParam('binFactor', params.FloatParam, default=1.,
                       label='Binning factor',
                       help='Bin the micrographs this much by a windowing '
                            'operation in the Fourier Tranform. Binning at '
                            'this level is hard to un-do later on, but may be '
                            'useful to down-scale super-resolution images. '
                            'Float-values may be used. Do make sure though '
                            'that the resulting micrograph size is even.')

        group.addParam('gainRot', params.EnumParam, default=0,
                       choices=['No rotation (0)',
                                ' 90 degrees (1)',
                                '180 degrees (2)',
                                '270 degrees (3)'],
                       label='Gain rotation',
                       help="Rotate the gain reference by this number times 90 "
                            "degrees clockwise in relion_display. This is the "
                            "same as -RotGain in MotionCor2. \n"
                            "Note that MotionCor2 uses a different convention "
                            "for rotation so it says 'counter-clockwise'.")

        group.addParam('gainFlip', params.EnumParam, default=0,
                       choices=['No flipping        (0)',
                                'Flip upside down   (1)',
                                'Flip left to right (2)'],
                       label='Gain flip',
                       help="Flip the gain reference after rotation. "
                            "This is the same as -FlipGain in MotionCor2. "
                            "0 means do nothing, 1 means flip Y (upside down) "
                            "and 2 means flip X (left to right).")

        form.addParam('extraParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="Extra parameters for Relion motion correction. ")

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Compute PSD (before/after)?",
                      help="If Yes, the protocol will compute for each movie "
                           "the average PSD before and after alignment, "
                           "for comparison")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail and keep it with the '
                           'micrograph object for visualization purposes. ')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions -------------------------------
    def _convertInputStep(self):
        self.info("Relion version:")
        self.runJob("relion_run_motioncorr --version", "", numberOfMpi=1)

        ProtAlignMovies._convertInputStep(self)

    def _processMovie(self, movie):
        movieFolder = self._getOutputMovieFolder(movie)
        inputStar = os.path.join(movieFolder,
                                 '%s_input.star' % self._getMovieRoot(movie))
        self.writeInputStar(inputStar, movie)

        pwutils.makePath(movieFolder, 'output')
        # The program will run in the movie folder, so let's put
        # the input files relative to that
        args = "--i %s --o output/ " % os.path.basename(inputStar)
        args += "--use_own "
        f0, fN = self._getRange(movie)
        args += "--first_frame_sum %d --last_frame_sum %d " % (f0, fN)
        args += "--bin_factor %f --bfactor %d " % (self.binFactor, self.bfactor)
        args += "--angpix %f " % (movie.getSamplingRate())
        args += "--patch_x %d --patch_y %d " % (self.patchX, self.patchY)
        args += "--j %d " % self.numberOfThreads

        inputMovies = self.inputMovies.get()
        if inputMovies.getGain():
            args += ' --gainref "%s"' % inputMovies.getGain()
            args += ' --gain_rot %d ' % self.gainRot
            args += ' --gain_flip %d ' % self.gainFlip

        if self.doDW:
            args += "--dose_weighting "
            if self.saveNonDW:
                args += " --save_noDW "

        preExp, dose = self._getCorrectedDose(self.inputMovies.get())
        voltage = movie.getAcquisition().getVoltage()
        args += "--voltage %d " % voltage
        args += "--dose_per_frame %f " % dose
        args += "--preexposure %f " % preExp

        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        self.runJob(self._getProgram(), args, cwd=movieFolder)

        self._computeExtra(movie)
        self._moveFiles(movie)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Zivanov2019']

    def _validate(self):
        # Check base validation before the specific ones for Motioncor
        errors = ProtAlignMovies._validate(self)
        return errors

    # ------------------------ Extra BASE functions ---------------------------
    def _getRelPath(self, baseName, refPath):
        return os.path.relpath(self._getExtraPath(baseName), refPath)

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _createOutputMovies(self):
        return True

    def _createOutputMicrographs(self):
        return not bool(self.doDW) or bool(self.saveNonDW)

    def _createOutputWeightedMicrographs(self):
        return bool(self.doDW)

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)
        self._setMotionValues(movie, mic)

    def _setMotionValues(self, movie, mic):
        """ Parse motion values from the 'corrected_micrographs.star' file
        generated for each movie. """
        fn = self._getMovieExtraFn(movie, 'corrected_micrographs.star')
        table = md.Table(fileName=fn)
        row = table[0]
        mic._rlnAccumMotionTotal = pwobj.Float(row.rlnAccumMotionTotal)
        mic._rlnAccumMotionEarly = pwobj.Float(row.rlnAccumMotionEarly)
        mic._rlnAccumMotionLate = pwobj.Float(row.rlnAccumMotionLate)

    def _getMovieShifts(self, movie, outStarFn=None):
        outStar = outStarFn or self._getMovieExtraFn(movie, '.star')
        first, last = self._getRange(movie)
        n = last - first + 1
        table = md.Table(fileName=outStar, tableName='global_shift')
        xShifts, yShifts = [], []

        for row in table:
            xShifts.append(float(row.rlnMicrographShiftX))
            yShifts.append(float(row.rlnMicrographShiftY))
            if len(xShifts) == n:
                break

        return xShifts, yShifts

    # --------------------------- UTILS functions -----------------------------
    def _getProgram(self, program='relion_run_motioncorr'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def writeInputStar(self, starFn, *images):
        """ Easy way to write a simple star file with a single micrographs.
        Used by the relion implementation of motioncor.
        """
        with open(starFn, 'w') as f:
            table = md.Table(columns=['rlnMicrographMovieName'])
            for img in images:
                table.addRow(os.path.basename(img.getFileName()))
            table.writeStar(f)

    def _getMovieOutFn(self, movie, suffix):
        movieBase = pwutils.removeBaseExt(movie.getFileName()).replace('.', '_')
        return os.path.join(self._getOutputMovieFolder(movie), 'output',
                            '%s%s' % (movieBase, suffix))

    def _getMovieExtraFn(self, movie, suffix):
        """ Return filenames in the extra directory with the prefix of this movie.
        Used to keep files associated with each micrograph.
        """
        movieBase = pwutils.removeBaseExt(movie.getFileName())
        return self._getExtraPath('%s%s' % (movieBase, suffix))

    def _getAbsPath(self, baseName):
        return os.path.abspath(self._getExtraPath(baseName))

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd_comparison', 'psd', extra=True)

    def _getPsdJpeg(self, movie):
        return self._getNameExt(movie, '_psd', 'jpeg', extra=True)

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = em.Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = em.Image(location=self._getPsdCorr(movie))
            mic.psdJpeg = em.Image(location=self._getPsdJpeg(movie))
        if self.doComputeMicThumbnail:
            mic.thumbnail = em.Image(
                location=self._getOutputMicThumbnail(movie))

    def _computeExtra(self, movie):
        """ Compute thumbnail, PSD and plots. """
        inputMovies = self.inputMovies.get()
        movieFolder = self._getOutputMovieFolder(movie)
        outMicFn = self._getMovieOutFn(movie, '.mrc')

        if self.doComputeMicThumbnail:
            self.computeThumbnail(outMicFn,
                                  outputFn=self._getOutputMicThumbnail(movie))

        if self.doComputePSD:
            #fakeShiftsFn = self.writeZeroShifts(movie)
            movieFn = movie.getFileName()
            aveMicFn = os.path.join(movieFolder,
                                    pwutils.removeBaseExt(movieFn) + "_tmp.mrc")
            self.averageMovie(movie, movieFn, aveMicFn,
                              binFactor=self.binFactor.get(),
                              dark=inputMovies.getDark(),
                              gain=inputMovies.getGain())

            self.computePSDs(movie, aveMicFn, outMicFn,
                             outputFnCorrected=self._getPsdJpeg(movie))

        self._saveAlignmentPlots(movie)

    def _moveFiles(self, movie):
        # It really annoying that Relion default names changes if you use DW or not
        # if use DW, the default name are DW and the others noDW
        if self.doDW:
            pwutils.moveFile(self._getMovieOutFn(movie, '.mrc'),
                             self._getExtraPath(self._getOutputMicWtName(movie)))
            if self.saveNonDW:
                pwutils.moveFile(self._getMovieOutFn(movie, '_noDW.mrc'),
                                 self._getExtraPath(self._getOutputMicName(movie)))
        else:
            pwutils.moveFile(self._getMovieOutFn(movie, '.mrc'),
                             self._getExtraPath(self._getOutputMicName(movie)))

        # Keep some local files of this movie in the extra folder
        for suffix in ['.star', '.log']:
            pwutils.moveFile(self._getMovieOutFn(movie, suffix),
                             self._getMovieExtraFn(movie, suffix))

        suffix = 'corrected_micrographs.star'
        fn = os.path.join(self._getOutputMovieFolder(movie), 'output', suffix)
        pwutils.moveFile(fn, self._getMovieExtraFn(movie, suffix))

    def _getRange(self, movie):
        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()
        first, last = self._getFrameRange(n, 'sum')

        if iniFrame != indxFrame:
            first -= iniFrame
            last -= iniFrame

        return first, last

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
        return movie.getNumberOfFrames()

    def _saveAlignmentPlots(self, movie):
        # Create plots and save as an image
        shiftsX, shiftsY = self._getMovieShifts(movie, self._getMovieOutFn(movie, '.star'))
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'align')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _createOutputMovie(self, movie):
        """ Overwrite this function to store the Relion's specific
        Motion model coefficients.
        """
        m = ProtAlignMovies._createOutputMovie(self, movie)
        # Load local motion values only if the patches are more than one
        if self.patchX.get() * self.patchY.get() > 1:
            table = md.Table(fileName=self._getMovieExtraFn(movie, '.star'),
                             tableName='local_motion_model')
            coeffs = [row.rlnMotionModelCoeff for row in table]
            m._rlnMotionModelCoeff = pwobj.String(json.dumps(coeffs))
        return m


def createGlobalAlignmentPlot(meanX, meanY, first):
    """ Create a plotter with the shift per frame. """
    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Global shift')
    ax.set_xlabel('Shift x (pixels)')
    ax.set_ylabel('Shift y (pixels)')

    i = first
    skipLabels = ceil(len(meanX)/10.0)
    labelTick = 1

    for x, y in izip(meanX, meanY):
        if labelTick == 1:
            ax.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    ax.plot(meanX, meanY, color='b')
    ax.plot(meanX, meanY, 'yo')

    plotter.tightLayout()

    return plotter
