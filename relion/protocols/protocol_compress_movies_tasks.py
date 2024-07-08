# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@gmail.com) [1]
# *
# * [1] St.Jude Children's Research Hospital, Memphis, TN
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

from emtools.utils import Timer, Pretty
from emtools.jobs import Pipeline
from emtools.pwx import SetMonitor, BatchManager
from emtools.metadata import StarFile, Table

from pyworkflow import SCIPION_DEBUG_NOCLEAN
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pyworkflow.constants import BETA
from pwem.protocols import ProtProcessMovies
from pwem.objects import MovieAlignment, SetOfMovies, ImageDim, FramesRange
from pyworkflow.protocol import STEPS_SERIAL


class ProtRelionCompressMoviesTasks(ProtProcessMovies):
    """
    Using *relion_convert_to_tiff* to compress a set of movies.
    """
    _label = 'compress movies (tasks)'
    _devStatus = BETA

    def __init__(self, **kwargs):
        ProtProcessMovies.__init__(self, **kwargs)
        self.isEER = False
        # We don't need parallelization at the steps level
        # we will use Pipeline/Tasks
        self.stepsExecutionMode = STEPS_SERIAL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff',
                               '.tif', '.eer', '.gain'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        ProtProcessMovies._defineParams(self, form)

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
                           "fraction. This option is relevant only for Falcon "
                           "movies in the EER format. Fractionate such "
                           "that each fraction has about 0.5 to 1.25 e/A2.")
        form.addParam('eerSampling', params.EnumParam, default=0,
                      choices=['1x', '2x'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='EER upsampling',
                      help="EER upsampling (1 = 4K or 2 = 8K). See "
                           "https://relion.readthedocs.io/en/latest/Reference/MovieCompression.html")

        form.addParallelSection(threads=4, mpi=0)

        self._defineStreamingParams(form)
        # Make default 1 minute for sleeping when no new input movies
        form.getParam('streamingSleepOnWait').setDefault(30)

    # We are not using the steps mechanism for parallelism from Scipion
    def _stepsCheck(self):
        pass

    @classmethod
    def worksInStreaming(cls):
        return True
    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self.samplingRate = self.inputMovies.get().getSamplingRate()
        self._insertFunctionStep(self._processAllMoviesStep)

    def _linkGain(self):
        gainFile = None
        inputGain = None

        if self.inputGainProt.get():
            inputGainProt = self.inputGainProt.get()
            inputGainProt._createFilenameTemplates()
            inputGainBin = inputGainProt._getFileName("output_gain")
            inputGainMrc = inputGainProt._getFileName("output_gain_mrc")
            inputDefects = inputGainProt._getFileName("output_gain_extra")

            if os.path.exists(inputGainMrc):
                inputGain = inputGainMrc
                gainBase = os.path.basename(inputGain)
                gainFile = self._getPath(gainBase)
            elif os.path.exists(inputGainBin):  # .bin gain files
                inputGain = inputGainBin
                gainBase = os.path.basename(inputGain)
                gainFile = self._getPath(gainBase)
                pwutils.createLink(inputDefects,
                                   gainFile.replace('.bin', '_reliablity.bin'))

        elif gainFn := self.inputMovies.get().getGain():
            inputGain = gainFn
            gainBase = os.path.basename(gainFn)
            gainFile = self._getPath(gainBase)

        if inputGain and os.path.exists(inputGain) and not os.path.exists(gainFile):
            pwutils.createLink(inputGain, gainFile)

        return gainFile

    def _processAllMoviesStep(self):
        self.info("Relion version:")
        self._runProgram('--version')

        moviesMtr = SetMonitor(SetOfMovies,
                               self.inputMovies.get().getFileName(),
                               blacklist=getattr(self, 'outputMovies', None))
        moviesIter = moviesMtr.iterProtocolInput(self, 'movies',
                                                 waitSecs=self.streamingSleepOnWait.get())
        batchMgr = BatchManager(self.streamingBatchSize.get(), moviesIter,
                                self._getTmpPath())

        self._outputMovies = None
        self._gainFile = self._linkGain()
        self.cmd = self._getCmd()

        pipe = Pipeline()
        g = pipe.addGenerator(batchMgr.generate)
        outputQueue = None
        for i in range(self.numberOfThreads.get()):
            proc = pipe.addProcessor(g.outputQueue, self._processBatch,
                                     outputQueue=outputQueue)
            outputQueue = proc.outputQueue

        pipe.addProcessor(outputQueue, self._outputFromBatch)
        pipe.run()

        for batch in batchMgr.generate():
            self._processBatch(batch)

        self._updateOutputSet('outputMovies', self._outputMovies,
                              pwobj.Set.STREAM_CLOSED)

    def _processBatch(self, batch):
        try:
            self.info(pwutils.cyanStr(f">>> Processing batch {batch['path']}"))
            batchPath = batch['path']
            starFn = os.path.join(batchPath, 'movies.star')
            with StarFile(starFn, 'w') as sf:
                t = Table(['rlnMicrographMovieName'])
                for movie in batch['items']:
                    fn = movie.getFileName()
                    bn = os.path.basename(fn)
                    pwutils.createLink(fn, os.path.join(batchPath, bn))
                    t.addRowValues(bn)
                sf.writeTable('movies', t)

            self._runProgram(self.cmd, cwd=batch['path'])

            # Check resulting files and update movies
            for movie in batch['items']:
                fn = movie.getFileName()
                tifBn = pwutils.replaceExt(os.path.basename(fn), 'tif')
                outputFn = os.path.join(batchPath, tifBn)
                dstFn = self._getExtraPath(tifBn)
                if os.path.exists(outputFn):
                    pwutils.moveFile(outputFn, dstFn)
                    movie.setFileName(dstFn)
                else:
                    movie.setFileName(None)

            gain = 'gain-reference.mrc'
            outputGain = os.path.join(batchPath, gain)
            newGain = self._getExtraPath(gain)
            if os.path.exists(outputGain) and not os.path.exists(newGain):
                pwutils.moveFile(outputGain, newGain)

            # Clean batch folder if not in debug mode
            if not pwutils.envVarOn(SCIPION_DEBUG_NOCLEAN):
                os.system('rm -rf %s' % batchPath)

        except Exception as e:
            eStr = str(e)
            self.error("ERROR: relion_convert_to_tiff has failed for batch %s. --> %s\n"
                       % (batch['id'], eStr))
            batch['error'] = eStr
            import traceback
            traceback.print_exc()

        return batch

    def _outputFromBatch(self, batch):
        # First time we are running this function for this execution
        firstOutput = False

        if self._outputMovies is None:
            outputMovies = getattr(self, 'outputMovies', None)
            if outputMovies is None:  # there is no previous output
                outputMovies = self._createSetOfMovies()
                outputMovies.setStreamState(pwobj.Set.STREAM_OPEN)
                outputMovies.copyInfo(self.inputMovies.get())
                m = batch['items'][0]
                dim = m.getDim()
                outputMovies.setDim(dim)  # Clear image dim
                framesRange = [1, dim[2], 1]
                acq = outputMovies.getAcquisition()
                newDose = acq.getDosePerFrame() * self.eerGroup.get()
                acq.setDosePerFrame(newDose)
                outputMovies.setFramesRange(framesRange)
                outputGain = self._getExtraPath('gain-reference.mrc')
                if os.path.exists(outputGain):
                    outputMovies.setGain(outputGain)
                firstOutput = True
            else:
                outputMovies.loadAllProperties()

            self._outputMovies = outputMovies
        else:
            outputMovies = self._outputMovies

        outputMovies.enableAppend()
        acq = outputMovies.getAcquisition()
        framesRange = outputMovies.getFramesRange()

        for movie in batch['items']:
            if movie.getFileName():
                # Fix acq and frames range based on grouping and new dose
                movie.setAcquisition(acq)
                movie.setFramesRange(framesRange)
                self._outputMovies.append(movie)

        self._updateOutputSet('outputMovies', outputMovies,
                              pwobj.Set.STREAM_OPEN)
        if firstOutput:
            self._defineSourceRelation(self.inputMovies, outputMovies)

    def _processMovie(self, movie):
        raise Exception("Not processing individual movies.")

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

        errors.extend(ProtProcessMovies._validate(self))

        return errors

    def _warnings(self):
        warnings = []
        if self.isEER:
            warnings.append("Note that after compression into TIFF the "
                            "original EER gain reference should be inverted.")

        return warnings

    # --------------------------- UTILS functions -----------------------------
    def _getCmd(self):
        """ Set return a command string that will be used for each batch. """
        compression = self.getEnumText('compression')
        cmd = " --i movies.star --o ./ "
        cmd += " --compression %s" % compression

        # TODO: Check if deflateLevel is only valid for zip (deflate)
        if compression == 'zip':  # deflate
            cmd += " --deflate_level %d" % self.deflateLevel

        # Gain file is expected at the run working folder
        # so, two levels up from tmp batch folder
        if self._gainFile:
            cmd += " --gain ../../" + os.path.basename(self._gainFile)

        if self.isEER:
            cmd += " --eer_grouping %d" % self.eerGroup
            cmd += " --eer_upsampling %d" % (self.eerSampling.get() + 1)

        return cmd

    def _runProgram(self, cmd, **kwargs):
        # We are using Scipion parallelization in batches, so not using MPI here
        self.runJob('relion_convert_to_tiff', cmd, numberOfMpi=1, **kwargs)
