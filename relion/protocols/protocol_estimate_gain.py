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
    Estimates a detector gain reference from a collection of cryo-EM movies
    using Relion gain estimation utilities. The protocol generates a gain map
    that can later be applied during motion correction and preprocessing to
    compensate for detector-specific intensity variations and improve the
    consistency of downstream image analysis.

    AI Generated:

    Estimate Gain Reference (ProtRelionCompressEstimateGain) — User Manual
        Overview

        The Estimate Gain Reference protocol computes a detector gain
        correction map directly from a set of experimental movies. In cryo-EM
        workflows, gain correction is an essential preprocessing step because
        direct electron detectors often exhibit pixel-to-pixel sensitivity
        differences that can introduce systematic artifacts into recorded data.
        By estimating a reliable gain reference from the acquired movies, the
        protocol helps normalize detector response and improve the quality of
        subsequent image processing tasks such as motion correction, CTF
        estimation, particle extraction, and reconstruction.

        This protocol is particularly useful when a calibrated gain reference
        is unavailable, incomplete, or suspected to be inaccurate. Biological
        users commonly apply it during the early stages of data processing,
        especially for newly acquired datasets or when working with detector
        formats that require specialized correction procedures.

        Input Movies and Data Requirements

        The protocol requires a set of cryo-EM movies acquired under relatively
        stable detector conditions. Ideally, all movies should originate from
        the same imaging session and detector configuration so that the
        estimated gain reflects a consistent detector response. Mixing movies
        collected under substantially different acquisition settings may reduce
        the reliability of the resulting gain map.

        The method supports several common movie formats used in cryo-EM,
        including MRC, TIFF, and EER data. Different detector formats may
        require different estimation strategies, particularly when dealing with
        floating-point images or electron event representation data. The
        protocol automatically adapts to the characteristics of the input data
        in order to generate an appropriate gain estimation workflow.

        Gain Estimation Strategy

        The protocol derives the detector gain by analyzing statistical
        properties across many movie frames. The underlying assumption is that
        detector-specific intensity variations remain relatively stable across
        acquisitions, allowing the protocol to estimate a correction pattern
        from accumulated observations.

        In practical biological workflows, using a sufficiently large and
        representative dataset improves the robustness of the gain estimate.
        Small datasets may produce noisier gain references, while datasets with
        strong contamination, severe ice gradients, or excessive artifacts may
        bias the estimation process.

        For floating-point movie formats, the protocol can apply additional
        reliability criteria when evaluating detector pixels. These controls
        help identify stable detector behavior and reduce the influence of
        inconsistent measurements. Advanced users may adjust reliability
        thresholds in challenging datasets, although the default settings are
        generally appropriate for routine processing.

        Working With EER Data

        The protocol supports EER movie formats commonly generated by modern
        direct electron detectors. In these cases, users may select different
        upsampling modes depending on whether they wish to estimate the gain at
        physical or super-resolution sampling levels.

        From a biological perspective, the selected sampling mode should remain
        consistent with the intended downstream processing workflow. If movies
        will later be processed in super-resolution mode, estimating the gain
        under matching conditions helps maintain consistency throughout the
        preprocessing pipeline.

        Frame Selection and Dataset Sampling

        Users may control how many frames contribute to the estimation process.
        Using all available frames generally produces the most stable gain map,
        especially for large datasets. However, limiting the number of frames
        can reduce runtime and computational cost during exploratory analyses
        or rapid facility screening workflows.

        The protocol also allows randomized selection of input movies before
        processing. Randomization can help reduce potential acquisition-order
        bias in very large datasets, particularly when imaging conditions drift
        gradually during data collection sessions.

        Computational Considerations

        Gain estimation is primarily an input-output and CPU-oriented task
        rather than a GPU-intensive operation. Nevertheless, the runtime may
        increase substantially for very large movie collections or high-
        resolution detector formats. Efficient storage access and sufficient
        threading resources can significantly improve performance.

        For large-scale facility processing, it is often practical to begin
        with a moderate subset of movies to validate the estimation procedure
        before computing a final gain reference from the complete dataset.

        Outputs and Their Interpretation

        The main output of the protocol is an estimated detector gain reference
        that can be used in subsequent preprocessing steps. Depending on the
        acquisition format and estimation strategy, additional reliability maps
        or intermediate correction files may also be generated.

        Biologically, the gain reference itself does not contain structural
        information about the specimen. Instead, it represents a calibration
        component that improves the accuracy and consistency of all downstream
        image analyses. Poor gain correction can introduce subtle image
        artifacts that propagate through the entire reconstruction workflow,
        potentially affecting map quality and interpretability.

        Practical Recommendations

        In routine cryo-EM practice, gain estimation works best when performed
        on datasets with stable detector behavior and sufficient sampling of
        experimental conditions. Users should visually inspect corrected
        micrographs after applying the estimated gain to confirm that detector
        artifacts have been adequately removed.

        When processing EER datasets, maintaining consistency between gain
        estimation settings and downstream reconstruction settings is highly
        recommended. For exploratory analyses or troubleshooting sessions,
        smaller subsets of movies may be sufficient, but final production
        workflows generally benefit from larger and more representative input
        collections.

        Final Perspective

        Accurate gain correction is a foundational step in cryo-EM image
        preprocessing because it directly influences the quality of all later
        analyses. Reliable detector normalization improves motion correction,
        enhances CTF estimation stability, and contributes to cleaner particle
        images and higher-quality reconstructions. Although often considered a
        technical preprocessing task, gain estimation plays an important role
        in ensuring the biological reliability of the final structural results.
    """
    _label = 'estimate gain reference'
    _devStatus = PROD
    stepsExecutionMode = STEPS_SERIAL


def __init__(self, **kwargs):
    ProtProcessMovies.__init__(self, **kwargs)
    self.isFloat32 = False
    self.isEER = False


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
                             self.inputMovies.getObjId(),
                             needsGPU=False)
    self._insertFunctionStep(self.estimateGainStep, needsGPU=False)


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
            f"--eer_upsampling {self.eerSampling.get() + 1}",
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
