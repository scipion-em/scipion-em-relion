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
    Compresses cryo-EM movie datasets into TIFF format using RELION
    movie conversion utilities. The protocol is intended to reduce
    storage requirements while preserving the information needed for
    downstream motion correction, refinement, and high-resolution
    reconstruction workflows.

    AI Generated:

    Compress Movies Tasks (ProtRelionCompressMoviesTasks) — User Manual
        Overview

        The Compress Movies Tasks protocol converts cryo-EM movies into
        compressed TIFF representations optimized for efficient storage
        and data management. Its main objective is to decrease disk
        usage while maintaining compatibility with standard RELION-based
        image processing pipelines.

        In modern cryo-EM facilities, movie datasets often occupy many
        terabytes of storage due to the large number of recorded frames
        and the increasing detector sizes used during acquisition.
        Compression therefore becomes an important practical step for
        long-term storage, data transfer, and computational efficiency.
        This protocol is especially useful in large screening projects,
        facility-scale processing environments, and workflows involving
        extensive archival requirements.

        Inputs and General Workflow

        The protocol requires a set of input movies that will be
        converted into compressed TIFF files. The resulting movies
        preserve the temporal and structural information needed for
        downstream cryo-EM analysis while reducing the physical storage
        footprint.

        The protocol supports streaming execution, allowing movies to be
        processed continuously as they become available. This capability
        is particularly valuable in automated acquisition pipelines
        where compression can occur in parallel with ongoing microscope
        data collection.

        Biological users should understand that compression changes the
        storage representation of the data but is not intended to alter
        the biological interpretation of the recorded signal. The output
        movies remain suitable for standard downstream procedures such
        as motion correction, CTF estimation, particle extraction, and
        reconstruction.

        TIFF Compression Strategies

        Several compression strategies are available depending on the
        balance desired between storage efficiency and processing speed.
        Lossless approaches are generally preferred in cryo-EM because
        they preserve the original detector information without
        introducing irreversible distortions.

        Automatic compression modes are convenient for routine
        workflows, while explicit ZIP or LZW compression can be useful
        when users wish to optimize storage efficiency or compatibility
        with specific infrastructures. Higher compression levels may
        provide smaller file sizes but can increase computational time
        during conversion.

        In most biological workflows, moderate compression settings
        offer the best compromise between storage reduction and runtime
        efficiency.

        Gain Reference Handling

        The protocol can incorporate detector gain references generated
        either externally or from dedicated gain estimation workflows.
        Proper gain handling is biologically important because detector
        response variations influence image uniformity and quantitative
        signal interpretation.

        When gain references are available, they are associated with the
        compressed dataset to preserve compatibility with downstream
        processing steps. This is particularly important in high-
        resolution cryo-EM workflows where subtle detector corrections
        can influence the final reconstruction quality.

        EER Movie Support

        The protocol provides dedicated support for Falcon EER movie
        formats. EER data contain highly fractionated temporal
        information that enables flexible dose grouping during
        processing. The protocol allows users to define grouping and
        sampling strategies appropriate for their imaging conditions.

        From a biological perspective, frame grouping affects the
        balance between temporal resolution and signal-to-noise ratio.
        Excessively fine grouping may preserve motion information but
        increase noise, whereas overly coarse grouping may obscure rapid
        beam-induced motion. Practical grouping choices usually depend
        on total dose, detector performance, and specimen sensitivity.

        The protocol also supports different EER upsampling modes,
        allowing users to adapt the effective detector sampling to their
        experimental requirements and computational constraints.

        Streaming and Parallel Processing

        One of the major strengths of this protocol is its ability to
        process movies in parallel batches while supporting continuous
        streaming operation. This design is particularly useful in
        modern cryo-EM infrastructures where datasets are generated
        continuously during automated acquisition sessions.

        By organizing movies into processing batches, the protocol
        improves throughput and allows efficient utilization of
        computational resources. This is especially beneficial when
        handling very large datasets acquired over extended microscope
        sessions.

        Outputs and Their Interpretation

        The protocol produces a new set of compressed TIFF movies ready
        for downstream cryo-EM analysis. The resulting dataset preserves
        acquisition metadata, frame organization, and experimental
        information required for later processing stages.

        When gain references are included, they are propagated together
        with the compressed movie set to maintain consistency throughout
        the workflow. For EER datasets, the output reflects the selected
        grouping and sampling strategy.

        Biological interpretation of the movies remains unchanged, but
        the storage representation becomes more efficient and practical
        for long-term computational workflows.

        Practical Recommendations

        In routine cryo-EM practice, movie compression is most valuable
        immediately after acquisition or before transferring datasets to
        long-term storage systems. Early compression can substantially
        reduce storage demands while maintaining workflow compatibility.

        For high-resolution projects, users should generally prefer
        conservative, lossless compression strategies to avoid any risk
        of compromising downstream refinement quality. It is also
        advisable to validate compressed outputs on a subset of movies
        before converting very large datasets.

        Users working with EER data should carefully select frame
        grouping values that match specimen motion behavior and exposure
        conditions. Appropriate grouping often improves downstream
        motion correction stability and computational efficiency.

        Final Perspective

        Efficient movie storage management has become increasingly
        important as cryo-EM datasets continue to grow in size and
        complexity. The Compress Movies Tasks protocol provides a
        practical solution for reducing storage requirements while
        preserving compatibility with modern RELION-based workflows.
        Careful selection of compression settings, gain handling, and
        EER grouping strategies helps ensure reliable and efficient
        downstream cryo-EM analysis.
    """
    _label = 'compress movies (tasks)'
    _devStatus = BETA
    # We don't need parallelization at the steps level
    # we will use Pipeline/Tasks
    stepsExecutionMode = STEPS_SERIAL

    def __init__(self, **kwargs):
        ProtProcessMovies.__init__(self, **kwargs)
        self.isEER = False

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
        self._insertFunctionStep(self._processAllMoviesStep, needsGPU=False)

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
