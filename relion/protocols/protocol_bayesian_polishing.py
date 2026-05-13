# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# * Authors:     Grigory Sharov     (gsharov@mrc-lmb.cam.ac.uk) [2]
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

import os
import json
from enum import Enum
from emtable import Table

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtParticles
import pwem.emlib.metadata as md
from pwem.constants import ALIGN_PROJ
from pwem.objects import SetOfParticles

import relion.convert as convert
from .protocol_base import ProtRelionBase


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtRelionBayesianPolishing(ProtParticles, ProtRelionBase):
    """
    Performs Bayesian polishing of cryo-EM particle images using Relion in
    order to correct beam-induced particle motion at the individual particle
    level. The protocol improves particle alignment quality and enhances the
    high-resolution signal of reconstructed maps by estimating smooth and
    spatially coherent motion trajectories across movie frames.

    AI Generated:

    Bayesian Polishing (ProtRelionBayesianPolishing) — User Manual
        Overview

        The Bayesian Polishing protocol refines particle motion correction
        using the Bayesian framework implemented in Relion. Its primary goal
        is to compensate for beam-induced movement that occurs during cryo-EM
        image acquisition, thereby improving the quality of particle averages
        and the final reconstructed density maps. Unlike global motion
        correction approaches that operate at the micrograph level, this
        protocol estimates motion trajectories for each individual particle,
        allowing more accurate correction of local movements.

        In practical cryo-EM workflows, Bayesian polishing is commonly applied
        after an initial high-resolution refinement and often after CTF
        refinement. The protocol is especially valuable for datasets targeting
        near-atomic resolution, where subtle motion-related blurring can limit
        map interpretability. By improving frame weighting and trajectory
        estimation, polishing frequently leads to measurable gains in final FSC
        resolution and map sharpness.

        Inputs and Experimental Context

        The protocol requires aligned movies, refined particles, and a
        postprocessing result containing the relevant FSC information and mask.
        The aligned movies provide the temporal information needed to estimate
        particle trajectories, while the refined particles define the particle
        orientations and positions already obtained in previous refinement
        steps. The postprocessed reconstruction contributes the FSC-based
        weighting scheme used during frame combination.

        From a biological perspective, the quality of the input refinement is
        critical. Poor particle alignments or inaccurate CTF estimation will
        propagate into polishing and reduce the benefit of the procedure.
        Therefore, Bayesian polishing is typically considered an advanced
        refinement step rather than an initial preprocessing operation.

        Training and Motion Statistics Estimation

        The protocol supports two operational modes. The first mode estimates
        optimal motion regularization parameters from a representative subset
        of particles. This training stage evaluates the statistical properties
        of particle motion in the dataset and determines appropriate levels of
        smoothness and coherence for the trajectories.

        Biologically, this step adapts the polishing behavior to the specific
        characteristics of the experiment. Datasets collected with unstable ice,
        significant beam-induced motion, or large particles may require
        different regularization strengths than highly stable datasets.
        Estimating these parameters automatically often provides more reliable
        results than manually selected values.

        The second mode performs the actual polishing using either the trained
        parameters or user-provided values. During this process, the protocol
        estimates per-particle motion tracks and generates improved particle
        images suitable for downstream refinement.

        Motion Regularization Parameters

        Bayesian polishing relies on three biologically meaningful motion
        regularization parameters that control the expected behavior of particle
        trajectories over time and across the micrograph.

        The velocity parameter controls how rapidly particles are allowed to
        move between frames. Smaller values favor shorter and smoother motion
        trajectories, while larger values permit greater flexibility. In highly
        stable datasets, stronger regularization often improves consistency,
        whereas datasets with substantial beam-induced drift may benefit from
        more relaxed motion estimates.

        The divergence parameter determines how similar neighboring particle
        trajectories are expected to be within the same micrograph. Lower values
        enforce spatial coherence and are useful when particle movement is
        dominated by collective ice motion. Larger values allow particles to
        move more independently, which may better describe heterogeneous ice
        behavior or local deformations.

        The acceleration parameter controls how strongly abrupt changes in
        particle direction are penalized. Smaller values favor smoother and
        straighter trajectories over time, while larger values allow more
        irregular motion patterns. Excessively weak regularization may lead to
        overfitting noise rather than genuine particle movement.

        Frame Selection and Dose Considerations

        The protocol allows selection of the movie frames included in motion
        estimation and frame combination. This is biologically important because
        early movie frames often contain the highest-resolution information but
        may also experience stronger beam-induced movement. Later frames are
        generally more stable but increasingly affected by radiation damage.

        Choosing an appropriate frame range therefore depends on the balance
        between motion correction and radiation preservation. In many practical
        datasets, excluding severely damaged late frames improves the quality of
        the polished particles and downstream reconstructions.

        Extraction and Rescaling Strategy

        The protocol optionally allows particles to be extracted using a
        different box size or rescaled pixel size during polishing. These
        options are important for balancing computational cost and achievable
        resolution.

        Larger extraction boxes preserve more peripheral signal and are useful
        for large macromolecular assemblies or flexible complexes. However,
        excessively large boxes increase computational demands and may introduce
        unnecessary background noise. Downsampling can accelerate processing,
        particularly during exploratory refinements, although aggressive scaling
        may limit high-resolution recovery.

        The selected extraction and scaling strategy should remain consistent
        with the expected biological resolution and particle size.

        EER Data and Motion Metadata

        The protocol supports advanced movie metadata, including motion models,
        hot-pixel information, gain references, and EER acquisition schemes.
        This ensures compatibility with modern direct electron detector
        workflows and allows polishing to incorporate detailed acquisition
        information when available.

        Accurate handling of detector corrections and optics information is
        particularly important in high-resolution cryo-EM projects, where even
        small calibration inconsistencies may reduce reconstruction quality.

        Outputs and Their Interpretation

        After completion, the protocol produces a polished set of particles with
        improved frame alignment and optimized weighting. These particles are
        intended for subsequent refinement, classification, or reconstruction
        steps.

        Biologically, successful polishing is often reflected by improved map
        sharpness, enhanced side-chain visibility, clearer secondary structure,
        and higher FSC resolution estimates. However, the magnitude of the
        improvement depends strongly on the original data quality, particle
        behavior, and accuracy of the preceding refinement steps.

        In training mode, the protocol also reports optimized motion
        regularization parameters. These values can guide future polishing runs
        on the same dataset or similar acquisitions.

        Practical Recommendations

        In routine cryo-EM workflows, Bayesian polishing is typically performed
        after obtaining a stable high-resolution refinement and after CTF
        refinement has converged. Running the training stage first is generally
        recommended because it adapts the motion model to the specific dataset.

        For datasets with strong beam-induced motion, careful optimization of
        frame ranges and motion regularization parameters can substantially
        improve reconstruction quality. Conversely, for highly stable datasets,
        conservative polishing settings may already provide near-optimal
        results.

        Users should visually inspect the resulting reconstructions after
        polishing rather than relying exclusively on nominal FSC improvements.
        Genuine biological improvement is best assessed through clearer density
        features and improved interpretability of flexible or previously blurred
        regions.

        Final Perspective

        Bayesian polishing represents one of the most important refinement
        stages in modern high-resolution cryo-EM workflows. By modeling
        particle-specific motion in a statistically robust manner, it improves
        the consistency and interpretability of reconstructed structures.
        Careful preparation of the input data, appropriate motion
        regularization, and biologically informed frame selection are key
        factors for obtaining reliable and meaningful improvements in final
        cryo-EM maps.
    """

    _label = 'bayesian polishing'
    _devStatus = PROD
    _possibleOutputs = outputs

    OP_TRAIN = 0
    OP_POLISH = 1

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_mics': self._getPath('input_corrected_micrographs.star'),
            'input_particles': self._getPath('input_particles.star'),
            'bfactors': self._getExtraPath('bfactors.star'),
            'shiny': self._getExtraPath('shiny.star'),
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMovies', params.PointerParam, pointerClass='SetOfMovies',
                      important=True,
                      label='Input ALIGNED movies',
                      help='Provide a set of movies that have at '
                           'least global alignment information.')
        form.addParam('inputParticles', params.PointerParam,
                      important=True, pointerCondition='hasAlignmentProj',
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles from 3D auto-refine '
                           'or CTF refinement.')
        form.addParam('inputPostprocess', params.PointerParam,
                      important=True,
                      label='Input Postprocess',
                      pointerClass='ProtRelionPostprocess',
                      help='Select a PostProcess job. The mask used for this '
                           'postprocessing will be applied to the unfiltered '
                           'half-maps and should encompass the entire complex. '
                           'The resulting FSC curve will be used for weighting '
                           'the different frequencies.')
        line = form.addLine('Movie frames',
                            help='First and last frames to take into account '
                                 'in motion fit and combination step '
                                 '(first frame is 1, last is 0).')
        line.addParam('frame0', params.IntParam, default=1,
                      label='first')
        line.addParam('frameN', params.IntParam, default=0,
                      label='last')

        form.addParam('extrSize', params.IntParam, default=-1,
                      label="Extraction size (px in unbinned movie)",
                      help="Size of the extracted particles in the "
                           "unbinned original movie (in pixels). "
                           "This should be an even number.")
        form.addParam('rescaledSize', params.IntParam, default=-1,
                      label="Re-scaled size (px)",
                      help="The re-scaled value needs to be an even number.")

        form.addParam('saveFloat16', params.BooleanParam, default=True,
                      label="Write output in float16?",
                      help="Relion can write output images in float16 "
                           "MRC (mode 12) format to save disk space. "
                           "By default, float32 format is used.")

        form.addSection(label='Train or Polish')
        form.addParam('operation', params.EnumParam, default=1,
                      choices=['Train optimal parameters',
                               'Perform particle polishing'],
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Operation',
                      help="If *train optimal parameters* , then "
                           "relion_motion_refine will estimate optimal "
                           "parameter values for the three sigma values above "
                           "on a subset of the data (determined by the minimum "
                           "number of particles to be used below).\n\n"
                           "If *perform particle polishing* then "
                           "relion_motion_refine will be run to estimate "
                           "per-particle motion-tracks using the parameters "
                           "below, and polished particles will be generated. ")

        condTrain = "operation==%s" % self.OP_TRAIN
        group = form.addGroup('Train', condition=condTrain)
        group.addParam('fractionFourierPx', params.FloatParam, default=0.5,
                       label='Fraction of Fourier pixels for testing',
                       help="This fraction of Fourier pixels (at higher "
                            "resolution) will be used for evaluation of the "
                            "parameters (test set), whereas the rest (at lower "
                            "resolution) will be used for parameter estimation "
                            "itself (work set).")
        group.addParam('numberOfParticles', params.IntParam, default=10000,
                       label='Use this many particles',
                       help='Use at least this many particles for the '
                            'meta-parameter optimisation. The more particles '
                            'the more expensive in time and computer memory '
                            'the calculation becomes, but the better the results '
                            'may get.')

        condPolish = "operation==%s" % self.OP_POLISH
        group = form.addGroup('Polish', condition=condPolish)
        group.addParam('sigmaVel', params.FloatParam, default=0.2,
                       label='Sigma for velocity (A/dose)',
                       help='Standard deviation for the velocity regularisation. '
                            'Smaller values requires the tracks to be shorter.')
        group.addParam('sigmaDiv', params.FloatParam, default=5000,
                       label='Sigma for divergence (A)',
                       help='Standard deviation for the divergence of tracks '
                            'across the micrograph. Smaller values requires '
                            'the tracks to be spatially more uniform in a '
                            'micrograph.')
        group.addParam('sigmaAcc', params.FloatParam, default=2,
                       label='Sigma for acceleration (A/dose)',
                       help='Standard deviation for the acceleration '
                            'regularisation. Smaller values requires the '
                            'tracks to be straighter.')
        line = group.addLine("Resolution for B-factor fit (A)",
                             help='The minimum and maximum spatial frequencies '
                                  '(in Angstrom) used in the B-factor fit.'
                                  'If a negative value is given as the maximum,'
                                  'it is determined from the input FSC curve.')
        line.addParam('minResBfactor', params.FloatParam, default=20, label='min')
        line.addParam('maxResBfactor', params.FloatParam, default=-1, label='max')

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- STEPS functions -------------------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep(self.convertInputStep,
                                 self.inputMovies.get().getObjId(),
                                 self.inputParticles.get().getObjId(),
                                 self.inputPostprocess.get().getObjId(),
                                 needsGPU=False)
        self._insertFunctionStep(self.trainOrPolishStep, self.operation.get(),
                                 needsGPU=False)
        if self.operation == self.OP_POLISH:
            self._insertFunctionStep(self.createOutputStep, 3, needsGPU=False)

    def convertInputStep(self, movId, partId, postId):
        inputMovies = self.inputMovies.get()
        inputParts = self.inputParticles.get()
        imgStar = self._getFileName('input_particles')
        inputPartsFolder = self._getInputPath('particles')
        pwutils.makePath(inputPartsFolder)

        self.info("Converting set from '%s' into '%s'" %
                  (inputParts.getFileName(), imgStar))

        # Create the first row, later only the movieName will be updated
        xdim, ydim, ndim = inputMovies.getDim()
        acq = inputMovies.getAcquisition()
        doseRate = acq.getDosePerFrame()
        firstMovie = inputMovies.getFirstItem()
        a0, aN = firstMovie.getAlignment().getRange()
        moviesPixelSize = inputMovies.getSamplingRate()
        binningFactor = inputParts.getSamplingRate() / moviesPixelSize

        og = convert.OpticsGroups.fromImages(inputMovies)
        writer = convert.createWriter(optics=og)
        writer.writeSetOfMicrographs(inputMovies,
                                     self._getFileName('input_mics'),
                                     postprocessImageRow=self._updateMic)

        # Handle EER case
        isEER = False
        if og.hasColumn('rlnEERGrouping'):
            isEER = True
            eerGrouping = og.first().rlnEERGrouping
            eerSampling = og.first().rlnEERUpsampling
            ndim //= eerGrouping
            doseRate *= eerGrouping

        generalCols = ['rlnImageSizeX',
                       'rlnImageSizeY',
                       'rlnImageSizeZ',
                       'rlnMicrographMovieName',
                       'rlnMicrographBinning',
                       'rlnMicrographOriginalPixelSize',
                       'rlnMicrographDoseRate',
                       'rlnMicrographPreExposure',
                       'rlnVoltage',
                       'rlnMicrographStartFrame',
                       'rlnMotionModelVersion',
                       'rlnMicrographGainName',
                       'rlnMicrographDefectFile']
        if isEER:
            generalCols.extend(['rlnEERGrouping', 'rlnEERUpsampling'])

        tableGeneral = Table(columns=generalCols)
        tableShifts = Table(columns=['rlnMicrographFrameNumber',
                                     'rlnMicrographShiftX',
                                     'rlnMicrographShiftY'])
        tableCoeffs = Table(columns=['rlnMotionModelCoeffsIdx',
                                     'rlnMotionModelCoeff'])
        tablePixels = Table(columns=['rlnCoordinateX',
                                     'rlnCoordinateY'])

        if not isEER:
            tableGeneral.addRow(xdim, ydim, ndim, 'movieName',
                                binningFactor, moviesPixelSize,
                                doseRate, acq.getDoseInitial(),
                                acq.getVoltage(), a0, 0, '""', '""')
        else:
            tableGeneral.addRow(xdim, ydim, ndim, 'movieName',
                                binningFactor, moviesPixelSize,
                                doseRate, acq.getDoseInitial(),
                                acq.getVoltage(), a0, 0, '""', '""',
                                eerGrouping, eerSampling)
        row = tableGeneral[0]

        for movie in inputMovies:
            movieStar = self._getMovieStar(movie)
            ogId = movie.getAttributeValue('_rlnOpticsGroup', 1)
            gainFn = og[ogId].get('rlnMicrographGainName', None)
            defectFn = og[ogId].get('rlnMicrographDefectFile', None)

            with open(movieStar, 'w') as f:
                coeffs = json.loads(movie.getAttributeValue('_rlnMotionModelCoeff', '[]'))
                motionMode = 1 if coeffs else 0
                hotpix = json.loads(movie.getAttributeValue('_rlnHotPixels', '[]'))

                # Update some params in the general table
                replaceDict = {'rlnMicrographMovieName': movie.getFileName(),
                               'rlnMotionModelVersion': motionMode}
                if gainFn:
                    replaceDict['rlnMicrographGainName'] = gainFn
                if defectFn:
                    replaceDict['rlnMicrographDefectFile'] = defectFn

                tableGeneral[0] = row._replace(**replaceDict)
                tableGeneral.writeStar(f, tableName='general', singleRow=True)
                # Write shifts
                tableShifts.clearRows()
                alignment = movie.getAlignment()
                shiftsX, shiftsY = alignment.getShifts()
                a0, aN = alignment.getRange()
                empty = -9999.000
                for i in range(1, a0):
                    tableShifts.addRow(i, empty, empty)
                # Adjust the shifts to be relative to the first frame
                # so let's add the opposite value
                xoff, yoff = -shiftsX[0], -shiftsY[0]
                for i in range(a0, aN + 1):
                    tableShifts.addRow(i, shiftsX[i-a0] + xoff,
                                       shiftsY[i-a0] + yoff)
                for i in range(aN + 1, ndim + 1):
                    tableShifts.addRow(i, empty, empty)
                tableShifts.writeStar(f, tableName='global_shift')

                # Write coefficients
                tableCoeffs.clearRows()
                if coeffs:
                    for i, c in enumerate(coeffs):
                        tableCoeffs.addRow(i, c)
                    tableCoeffs.writeStar(f, tableName='local_motion_model')

                # Write hot pixels
                tablePixels.clearRows()
                if hotpix:
                    for coord in hotpix:
                        tablePixels.addRow(coord[0], coord[1])
                    tablePixels.writeStar(f, tableName='hot_pixels')

        convert.writeSetOfParticles(inputParts, imgStar,
                                    outputDir=inputPartsFolder,
                                    alignType=ALIGN_PROJ,
                                    fillMagnification=True)

    def trainOrPolishStep(self, operation):
        postProt = self.inputPostprocess.get()
        args = "--i %s " % self._getFileName('input_particles')
        args += "--o %s " % self._getExtraPath()
        postStar = postProt._getExtraPath('postprocess.star')
        args += "--f %s " % postStar
        args += "--angpix_ref %0.5f " % postProt.outputVolume.getSamplingRate()
        args += "--corr_mic %s " % self._getFileName('input_mics')
        args += "--first_frame %d --last_frame %d " % (self.frame0, self.frameN)

        if self.extrSize.get() != -1:
            args += "--window %d " % self.extrSize.get()
        if self.rescaledSize.get() != -1:
            args += "--scale %d " % self.rescaledSize.get()

        if self.operation == self.OP_TRAIN:
            args += "--min_p %d " % self.numberOfParticles
            args += "--eval_frac %0.3f " % self.fractionFourierPx
            args += "--align_frac %0.3f " % self.fractionFourierPx
            args += "--params3 "
        else:  # OP_POLISH
            args += "--s_vel %0.3f " % self.sigmaVel
            args += "--s_div %0.3f " % self.sigmaDiv
            args += "--s_acc %0.3f " % self.sigmaAcc
            args += "--bfac_minfreq %0.3f " % self.minResBfactor
            args += "--bfac_maxfreq %0.3f " % self.maxResBfactor
            args += "--combine_frames "

        if self.saveFloat16:
            args += "--float16 "

        args += "--j %d " % self.numberOfThreads

        self.runJob(self._getProgram('relion_motion_refine'), args)

    def createOutputStep(self, id=1):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        pixSize = self._getOutputPixSize()
        outImgSet.setSamplingRate(pixSize)

        outImgsFn = md.MetaData('particles@' + self._getFileName('shiny'))
        rowIterator = md.SetMdIterator(outImgsFn, sortByLabel=md.RLN_IMAGE_ID,
                                       keyLabel=md.RLN_IMAGE_ID,
                                       updateItemCallback=self._updatePtcl)
        outImgSet.copyItems(imgSet,
                            updateItemCallback=rowIterator.updateItem)

        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(self.inputParticles, outImgSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        def _params(label, *params):
            summary.append('%s params:' % label)
            summary.append('    Sigma for velocity: *%0.3f*' % params[0])
            summary.append('    Sigma for divergence: *%0.1f*' % params[1])
            summary.append('    Sigma for acceleration: *%0.2f*' % params[2])

        if self.operation != self.OP_TRAIN:
            _params('Input', self.sigmaVel, self.sigmaDiv, self.sigmaAcc)
        else:
            outputFn = None
            for fn in ['opt_params.txt', 'opt_params_all_groups.txt']:
                if os.path.exists(self._getExtraPath(fn)):
                    outputFn = self._getExtraPath(fn)

            if outputFn is None:
                summary.append('Output is not ready yet.')
            else:
                with open(outputFn) as f:
                    line = [float(x) for x in f.readline().split()]
                    _params('Output', *line)

        return summary

    def _validate(self):
        errors = []

        win = self.extrSize.get()
        scale = self.rescaledSize.get()
        if win * scale <= 0:
            errors.append("Please specify both the extraction box size and "
                          "the downsampled size, or leave both the default (-1)")
        if win != -1 and scale != -1:
            if win % 2 != 0:
                errors.append("ERROR: The extraction box size must be an "
                              "even number")
            if scale % 2 != 0:
                errors.append("ERROR: The downsampled box size must be an "
                              "even number")
            if scale > win:
                errors.append("ERROR: The downsampled box size cannot be "
                              "larger than the extraction size")

        if self.operation == self.OP_TRAIN and self.numberOfMpi > 1:
            errors.append("MPI is not supported for parameters estimation.")

        return errors

    def _warnings(self):
        warnings = ['If you have provided a gain reference or defects file during '
                    'movie import or motion correction, please *make sure to '
                    'run first "assign optics groups" protocol for aligned '
                    'movies*, specifying the gain file etc. Currently, Scipion '
                    'has no other way of knowing if you have e.g. rotated the '
                    'gain during motion correction.\n\nOutput movies then can be '
                    'used in this polishing protocol.']

        return warnings

    # -------------------------- UTILS functions ------------------------------
    def _getInputPath(self, *paths):
        return self._getPath('input', *paths)

    def _getMovieStar(self, movie):
        return self._getInputPath(pwutils.replaceBaseExt(movie.getMicName(),
                                                         'star'))

    def _updatePtcl(self, particle, row):
        newLoc = convert.relionToLocation(row.getValue('rlnImageName'))
        particle.setLocation(newLoc)

    def _updateMic(self, mic, row):
        row['rlnMicrographName'] = os.path.basename(mic.getMicName())
        row['rlnMicrographMetadata'] = self._getMovieStar(mic)

    def _getOutputPixSize(self):
        parts = self.inputParticles.get()
        movies = self.inputMovies.get()

        if self.rescaledSize.get() == -1:
            # no scale or window, return particle pix size
            return parts.getSamplingRate()
        else:
            if self.rescaledSize.get() == self.extrSize.get():
                # window only, return movie pix size
                return movies.getSamplingRate()
            else:
                # rescale and window
                return movies.getSamplingRate() * self.extrSize.get() / self.rescaledSize.get()
