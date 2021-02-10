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
from math import ceil
import json
import emtable as md

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
from pwem.protocols import ProtAlignMovies
from pwem.objects import Image
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol import STEPS_SERIAL

import relion
import relion.convert as convert
from relion.convert.convert31 import OpticsGroups


class ProtRelionMotioncor(ProtAlignMovies):
    """
    Wrapper for the Relion's implementation of motioncor algorithm.
    """

    _label = 'motion correction'

    def __init__(self, **kwargs):

        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        self.updatedSets = []

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

        form.addParam('doDW', params.BooleanParam, default=True,
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

        if self.IS_GT30():
            form.addParam('savePSsum', params.BooleanParam, default=False,
                          label='Save sum of power spectra?',
                          help='Sum of non-dose weighted power spectra '
                               'provides better signal for CTF estimation. '
                               'The power spectra can be used by CTFFIND4 '
                               'but not by GCTF.')
            form.addParam('dosePSsum', params.FloatParam, default=4.0,
                          condition='savePSsum',
                          label='Sum power spectra every e/A2',
                          help='McMullan et al. (Ultramicroscopy, 2015) '
                               'suggests summing power spectra every '
                               '4.0 e/A2 gives optimal Thon rings.')

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
                           'micrograph object for visualization purposes.\n\n'
                           '*IMPORTANT: this requires EMAN2 plugin and binaries.*')

        form.addSection("Motion")
        form.addParam('bfactor', params.IntParam, default=150,
                      label='Bfactor',
                      help="The B-factor that will be applied to the "
                           "micrographs.")

        line = form.addLine('Number of patches',
                            help='Number of patches (in X and Y direction) to '
                                 'apply motion correction. If <= 2 then '
                                 'only global correction will be done.')
        line.addParam('patchX', params.IntParam, default=1, label='X')
        line.addParam('patchY', params.IntParam, default=1, label='Y')

        form.addParam('groupFrames', params.IntParam, default=1,
                      label='Group frames',
                      help="Average together this many frames before "
                           "calculating the beam-induced shifts.")

        form.addParam('binFactor', params.FloatParam, default=1.,
                      label='Binning factor',
                      help='Bin the micrographs this much by a windowing '
                           'operation in the Fourier Tranform. Binning at '
                           'this level is hard to un-do later on, but may be '
                           'useful to down-scale super-resolution images. '
                           'Float-values may be used. Do make sure though '
                           'that the resulting micrograph size is even.')

        form.addParam('gainRot', params.EnumParam, default=0,
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

        form.addParam('gainFlip', params.EnumParam, default=0,
                      choices=['No flipping        (0)',
                               'Flip upside down   (1)',
                               'Flip left to right (2)'],
                      label='Gain flip',
                      help="Flip the gain reference after rotation. "
                           "This is the same as -FlipGain in MotionCor2. "
                           "0 means do nothing, 1 means flip Y (upside down) "
                           "and 2 means flip X (left to right).")

        form.addParam('defectFile', params.FileParam, allowsNull=True,
                      label='Defects file',
                      help='Location of a UCSF MotionCor2-style '
                           'defect text file or a defect map that '
                           'describe the defect pixels on the detector. '
                           'Each line of a defect text file should contain '
                           'four numbers specifying x, y, width and height '
                           'of a defect region. A defect map is an image '
                           '(MRC or TIFF), where 0 means good and 1 means '
                           'bad pixels. The coordinate system is the same '
                           'as the input movie before application of '
                           'binning, rotation and/or flipping.\n\n'
                           '_Note that the format of the defect text is '
                           'DIFFERENT from the defect text produced '
                           'by SerialEM!_\n One can convert a SerialEM-style '
                           'defect file into a defect map using IMOD '
                           'utilities e.g.:\n'
                           '*clip defect -D defect.txt -f tif movie.tif defect_map.tif*\n'
                           'See explanations in the SerialEM manual.\n'
                           'Leave empty if you do not have any defects, '
                           'or do not want to correct for defects on your detector.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions -------------------------------
    def _convertInputStep(self):
        self.info("Relion version:")
        self.runJob("relion_run_motioncorr --version", "", numberOfMpi=1)
        self.info("Detected version from config: %s" % relion.Plugin.getActiveVersion())

        ProtAlignMovies._convertInputStep(self)

    def _processMovie(self, movie):
        movieFolder = self._getOutputMovieFolder(movie)
        inputStar = os.path.join(movieFolder,
                                 '%s_input.star' % self._getMovieRoot(movie))
        pwutils.makePath(os.path.join(movieFolder, 'output'))

        og = OpticsGroups.fromImages(self.inputMovies.get())
        writer = convert.createWriter(optics=og)
        # Let's use only the basename, since we will launch the command
        # from the movieFolder
        movie.setFileName(os.path.basename(movie.getFileName()))
        writer.writeSetOfMovies([movie], inputStar)

        # The program will run in the movie folder, so let's put
        # the input files relative to that
        args = "--i %s --o output/ " % os.path.basename(inputStar)
        args += "--use_own "
        f0, fN = self._getRange(movie)
        args += "--first_frame_sum %d --last_frame_sum %d " % (f0, fN)
        args += "--bin_factor %f --bfactor %d " % (self.binFactor, self.bfactor)
        args += "--angpix %0.5f " % (movie.getSamplingRate())
        args += "--patch_x %d --patch_y %d " % (self.patchX, self.patchY)
        args += "--group_frames %d " % self.groupFrames
        args += "--j %d " % self.numberOfThreads

        inputMovies = self.inputMovies.get()
        if inputMovies.getGain():
            args += ' --gainref "%s" ' % inputMovies.getGain()
            args += ' --gain_rot %d ' % self.gainRot
            args += ' --gain_flip %d ' % self.gainFlip

        if self.IS_GT30():
            if self.defectFile.get():
                args += ' --defect_file "%s" ' % self.defectFile.get()

            if self._savePsSum():
                args += ' --grouping_for_ps %d ' % self._calcPsDose()

        if self.doDW:
            args += "--dose_weighting "
            preExp, dose = self._getCorrectedDose(self.inputMovies.get())
            args += "--dose_per_frame %f " % dose
            args += "--preexposure %f " % preExp

            if self.saveNonDW:
                args += " --save_noDW "

        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        try:
            self.runJob(self._getProgram(), args, cwd=movieFolder)

            self._computeExtra(movie)
            self._moveFiles(movie)
        except:
            print("ERROR: processing movie: ", movie.getFileName())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Zivanov2019']

    def _validate(self):
        # Check base validation before the specific ones for Motioncor
        errors = ProtAlignMovies._validate(self)

        if not relion.Plugin.getActiveVersion():
            errors.append("Could not detect the current Relion version. \n"
                          "RELION_HOME='%s'" % relion.Plugin.getHome())

        acq = self.inputMovies.get().getAcquisition()
        if self.doDW:
            dose = acq.getDosePerFrame()
            if dose is None or dose < 0.001:
                errors.append("Input movies do not contain the dose per frame, "
                              "dose-weighting can not be performed. ")

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

    def _savePsSum(self):
        return self.getAttributeValue('savePSsum', False)

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)
        self._setMotionValues(movie, mic)
        if self._savePsSum():
            outPs = self._getExtraPath(self._getOutputMicPsName(movie))
            mic._powerSpectra = Image(location=outPs)
            mic._powerSpectra.setSamplingRate(self._calcPSSampling())

    def _setMotionValues(self, movie, mic):
        """ Parse motion values from the 'corrected_micrographs.star' file
        generated for each movie. """
        fn = self._getMovieExtraFn(movie, 'corrected_micrographs.star')
        micsTableName = 'micrographs' if self.IS_GT30() else ''
        table = md.Table(fileName=fn, tableName=micsTableName)
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
            # Shifts are in pixels of the original (unbinned) movies
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
        mic.plotGlobal = Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = Image(location=self._getPsdCorr(movie))
            mic.psdJpeg = Image(location=self._getPsdJpeg(movie))
        if self.doComputeMicThumbnail:
            mic.thumbnail = Image(
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
            movieFn = os.path.join(movieFolder, movie.getFileName())
            aveMicFn = os.path.join(movieFolder,
                                    pwutils.removeBaseExt(movieFn) + "_tmp.mrc")
            self.averageMovie(movie, movieFn, aveMicFn,
                              binFactor=self.binFactor.get(),
                              dark=inputMovies.getDark(),
                              gain=inputMovies.getGain())

            self.computePSDImages(movie, aveMicFn, outMicFn,
                                  outputFnCorrected=self._getPsdJpeg(movie))

        self._saveAlignmentPlots(movie, inputMovies.getSamplingRate())

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

        if self._savePsSum():
            pwutils.moveFile(self._getMovieOutFn(movie, '_PS.mrc'),
                             self._getExtraPath(self._getOutputMicPsName(movie)))

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

    def _saveAlignmentPlots(self, movie, pixSize):
        # Create plots and save as an image
        shiftsX, shiftsY = self._getMovieShifts(movie, self._getMovieOutFn(movie, '.star'))
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'sum')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first, pixSize)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _createOutputMovie(self, movie):
        """ Overwrite this function to store the Relion's specific
        Motion model coefficients.
        """
        m = ProtAlignMovies._createOutputMovie(self, movie)
        # Load local motion values only if the patches are more than one
        if self.patchX.get() > 2 and self.patchY.get() > 2:
            try:
                table = md.Table(fileName=self._getMovieExtraFn(movie, '.star'),
                                 tableName='local_motion_model')
                coeffs = [row.rlnMotionModelCoeff for row in table]
            except:
                print("Failed to parse local motion from: %s" %
                      os.path.abspath(self._getMovieExtraFn(movie, '.star')))
                coeffs = []  # Failed to parse the local motion
            m._rlnMotionModelCoeff = pwobj.String(json.dumps(coeffs))
        return m

    def createOutputStep(self):
        # This method is re-implemented here because a bug in the base protocol
        # where the outputMicrographs is used without check if it is produced.
        # validate that we have some output movies
        if self._createOutputMovies():
            output = self.outputMovies
        elif self._createOutputMicrographs():
            output = self.outputMicrographs
        elif self._createOutputWeightedMicrographs():
            output = self.outputMicrographsDoseWeighted
        else:
            raise Exception("It does not seem like any output is produced!")

        inputSize = len(self.listOfMovies)
        outputSize = output.getSize()

        if outputSize == 0 and inputSize != 0:
            raise Exception("All movies failed, didn't create outputMicrographs."
                            "Please review movie processing steps above.")

        if outputSize < inputSize:
            self.warning(pwutils.yellowStr("WARNING - Failed to align %d movies."
                                           % (inputSize - outputSize)))

    def _updateOutputSet(self, outputName, outputSet,
                         state=pwobj.Set.STREAM_OPEN):
        """ Redefine this method to update optics info. """

        if outputName not in self.updatedSets:
            og = OpticsGroups.fromImages(outputSet)
            og.updateAll(rlnMicrographOriginalPixelSize=self.inputMovies.get().getSamplingRate())
            og.toImages(outputSet)
            self.updatedSets.append(outputName)

        ProtAlignMovies._updateOutputSet(self, outputName, outputSet,
                                         state=state)

    def _calcPsDose(self):
        _, dose = self._getCorrectedDose(self.inputMovies.get())
        dose_for_ps = round(self.dosePSsum.get() / dose)

        return 1 if dose_for_ps == 0 else dose_for_ps

    def _calcPSSampling(self):
        """ Adapted from relion 3.1 code. """
        movieSet = self.inputMovies.get()
        target_pixel_size = 1.4  # from CTFFIND 4.1
        ps = movieSet.getSamplingRate()
        ps_angpix = ps * self.binFactor.get()
        x, y, _ = movieSet.getDimensions()
        ps_size_square = min(x, y)
        if ps_angpix < target_pixel_size:
            nx_needed = ceil(ps_size_square * ps_angpix / target_pixel_size)
            nx_needed += nx_needed % 2
            ps_angpix = ps_size_square * ps_angpix / nx_needed
        return ps_angpix

    def _getOutputMicPsName(self, movie):
        """ Returns the name of the output PS
        (relative to micFolder)
        """
        return self._getMovieRoot(movie) + '_aligned_mic_PS.mrc'

    def _getFrameRange(self, n, prefix):
        # Reimplement this method to ignore prefix (called from base class)
        # and always use 'sum' as prefix
        return ProtAlignMovies._getFrameRange(self, n, 'sum')

    def IS_GT30(self):
        return relion.Plugin.IS_GT30()


def createGlobalAlignmentPlot(meanX, meanY, first, pixSize):
    """ Create a plotter with the shift per frame. """
    sumMeanX = []
    sumMeanY = []

    def px_to_ang(ax_px):
        y1, y2 = ax_px.get_ylim()
        x1, x2 = ax_px.get_xlim()
        ax_ang2.set_ylim(y1 * pixSize, y2 * pixSize)
        ax_ang.set_xlim(x1 * pixSize, x2 * pixSize)
        ax_ang.figure.canvas.draw()
        ax_ang2.figure.canvas.draw()

    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax_px = figure.add_subplot(111)
    ax_px.grid()
    ax_px.set_xlabel('Shift x (px)')
    ax_px.set_ylabel('Shift y (px)')

    ax_ang = ax_px.twiny()
    ax_ang.set_xlabel('Shift x (A)')
    ax_ang2 = ax_px.twinx()
    ax_ang2.set_ylabel('Shift y (A)')

    i = first
    # The output _rlnMicrographShiftX/Y shifts relative to the first frame.
    # Unit is pixels of the original (unbinned) movies (Takanori, 2018)
    skipLabels = ceil(len(meanX) / 10.0)
    labelTick = 1

    for x, y in zip(meanX, meanY):
        sumMeanX.append(x)
        sumMeanY.append(y)
        if labelTick == 1:
            ax_px.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    # automatically update lim of ax_ang when lim of ax_px changes.
    ax_px.callbacks.connect("ylim_changed", px_to_ang)
    ax_px.callbacks.connect("xlim_changed", px_to_ang)

    ax_px.plot(sumMeanX, sumMeanY, color='b')
    ax_px.plot(sumMeanX, sumMeanY, 'yo')
    ax_px.plot(sumMeanX[0], sumMeanY[0], 'ro', markersize=10, linewidth=0.5)
    ax_px.set_title('Global frame alignment')

    plotter.tightLayout()

    return plotter
