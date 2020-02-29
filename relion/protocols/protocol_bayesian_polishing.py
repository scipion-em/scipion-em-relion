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
import json
from io import open

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pwem.protocols import ProtParticles
import pwem.emlib.metadata as md
from pwem.constants import ALIGN_PROJ

import relion.convert as convert
from ..convert.metadata import Table


class ProtRelionBayesianPolishing(ProtParticles):
    """
    Wrapper protocol for the Relion's Bayesian Polishing.

    As of release 3.0, Relion also implements a new Bayesian approach to beam
    induced motion correction. This approach aims to optimise a regularised
    likelihood, which allows us to associate with each hypothetical set of
    particle trajectories a prior likelihood that favors spatially coherent
    and temporally smooth motion without imposing any hard constraints.
    The smoothness prior term requires three parameters that describe the
    statistics of the observed motion. To estimate the prior that yields the
    best motion tracks for this particular dataset, we can first run the
    program in 'training mode'. Once the estimates have been obtained, one
    can then run the program again to fit tracks for the motion of all
    particles in the data set and to produce adequately weighted averages of
    the aligned movie frames.

    """

    _label = 'bayesian polishing'

    OP_TRAIN = 0
    OP_POLISH = 1

    def _defineParams(self, form):
        form.addSection(label='Input')

        # TODO: conditions on movies?
        form.addParam('inputMovies', params.PointerParam, pointerClass='SetOfMovies',
                      important=True,
                      label='Input movies',
                      help='Provide a set of movies that have at '
                           'least global alignment information.')
        # TODO: conditions on particles?
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
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
                           'the different frequencies. ')
        line = form.addLine('Movie frames',
                            help='First and last frames to take into account '
                                 'in motion fit and combination step. '
                                 '(starts counting at 1 and 0 as last '
                                 'means util the last frame in the movie).')
        line.addParam('frame0', params.IntParam, default=1,
                      label='first')
        line.addParam('frameN', params.IntParam, default=0,
                      label='last')

        form.addSection(label='Train or Polish')
        form.addParam('operation', params.EnumParam, default=0,
                      choices=[' Train optimal parameters ',
                               ' Perform particle polishing '],
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
        self._insertFunctionStep('convertInputStep',
                                 self.inputMovies.get().getObjId(),
                                 self.inputParticles.get().getObjId(),
                                 self.inputPostprocess.get().getObjId())
        self._insertFunctionStep('trainOrPolishStep', self.operation.get())
        if self.operation == self.OP_POLISH:
            self._insertFunctionStep('createOutputStep', 3)

    def _writeMovieStar(self, movie, starFn):
        """ Write the required """

    def convertInputStep(self, movId, partId, postId):
        inputMovies = self.inputMovies.get()
        inputParts = self.inputParticles.get()
        imgStar = self._getPath('input_particles.star')
        inputPartsFolder = self._getInputPath('particles')
        pwutils.makePath(inputPartsFolder)

        self.info("Converting set from '%s' into '%s'" %
                  (inputParts.getFileName(), imgStar))

        tableMovies = Table(columns=['rlnMicrographName',
                                     'rlnMicrographMetadata'])
        tableGeneral = Table(columns=['rlnImageSizeX',
                                      'rlnImageSizeY',
                                      'rlnImageSizeZ',
                                      'rlnMicrographMovieName',
                                      'rlnMicrographBinning',
                                      'rlnMicrographOriginalPixelSize',
                                      'rlnMicrographDoseRate',
                                      'rlnMicrographPreExposure',
                                      'rlnVoltage',
                                      'rlnMicrographStartFrame',
                                      'rlnMotionModelVersion'])
        tableShifts = Table(columns=['rlnMicrographFrameNumber',
                                     'rlnMicrographShiftX',
                                     'rlnMicrographShiftY'])
        tableCoeffs = Table(columns=['rlnMotionModelCoeffsIdx',
                                     'rlnMotionModelCoeff'])

        # Create the first row, later only the movieName will be updated
        xdim, ydim, ndim = inputMovies.getDim()
        acq = inputMovies.getAcquisition()
        firstMovie = inputMovies.getFirstItem()
        a0, aN = firstMovie.getAlignment().getRange()
        moviesPixelSize = inputMovies.getSamplingRate()
        binningFactor = inputParts.getSamplingRate() / moviesPixelSize
        hasLocal = firstMovie.hasAttribute('_rlnMotionModelCoeff')
        motionMode = 1 if hasLocal else 0

        tableGeneral.addRow(xdim, ydim, ndim, 'movieName',
                            binningFactor, moviesPixelSize,
                            acq.getDosePerFrame(), acq.getDoseInitial(),
                            acq.getVoltage(), a0, motionMode)
        row = tableGeneral[0]

        for movie in inputMovies:
            movieFn = movie.getFileName()
            movieBase = os.path.basename(movieFn)
            movieStar = self._getInputPath(pwutils.replaceBaseExt(movieFn,
                                                                  'star'))
            tableMovies.addRow(movieBase, movieStar)
            with open(movieStar, 'w') as f:
                # Update Movie name
                tableGeneral[0] = row._replace(rlnMicrographMovieName=movieFn)
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
                if hasLocal:
                    coeffs = movie.getAttributeValue('_rlnMotionModelCoeff', '')
                    tableCoeffs.clearRows()
                    for i, c in enumerate(json.loads(coeffs)):
                        tableCoeffs.addRow(i, c)
                    tableCoeffs.writeStar(f, tableName='local_motion_model')

        with open(self._getPath('input_corrected_micrographs.star'), 'w') as f:
            tableMovies.writeStar(f)

        convert.writeSetOfParticles(inputParts, imgStar,
                                    inputPartsFolder,
                                    alignType=ALIGN_PROJ,
                                    fillMagnification=True,
                                    fillRandomSubset=True)

    def trainOrPolishStep(self, operation):
        args = "--i %s " % self._getPath('input_particles.star')
        args += "--o %s " % self._getExtraPath()
        postStar = self.inputPostprocess.get()._getExtraPath('postprocess.star')
        args += "--f %s " % postStar
        postprocessTuple = convert.getVolumesFromPostprocess(postStar)
        args += "--m1 %s --m2 %s --mask %s " % postprocessTuple
        args += "--corr_mic %s " % self._getPath('input_corrected_micrographs.star')
        args += "--first_frame %d --last_frame %d " % (self.frame0, self.frameN)

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

        args += "--j %d " % self.numberOfThreads

        prog = "relion_motion_refine" + ("_mpi" if self.numberOfMpi > 1 else "")
        self.runJob(prog, args)

    def createOutputStep(self, id=1):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)

        outImgsFn = md.MetaData(self._getExtraPath('shiny.star'))
        rowIterator = md.SetMdIterator(outImgsFn, sortByLabel=md.RLN_IMAGE_ID,
                                       keyLabel=md.RLN_IMAGE_ID,
                                       updateItemCallback=self._updateItem)
        outImgSet.copyItems(imgSet,
                            updateItemCallback=rowIterator.updateItem)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

    def _updateItem(self, particle, row):
        newLoc = convert.relionToLocation(row.getValue('rlnImageName'))
        particle.setLocation(newLoc)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if self.operation != self.OP_TRAIN:
            summary.append('Sigma for velocity: %0.3f' % self.sigmaVel)
            summary.append('Sigma for divergence: %0.1f' % self.sigmaDiv)
            summary.append('Sigma for acceleration: %0.2f' % self.sigmaAcc)
        else:
            if pwutils.exists(self._getExtraPath('opt_params.txt')):
                f = open(self._getExtraPath('opt_params.txt'))
                line = list(map(float, f.readline().split()))
                summary.append('Sigma for velocity: %0.3f' % line[0])
                summary.append('Sigma for divergence: %0.1f' % line[1])
                summary.append('Sigma for acceleration: %0.2f' % line[2])
                f.close()
            else:
                summary.append('Output is not ready yet.')

        return summary

    def _validate(self):
        errors = []

        if self.operation == self.OP_TRAIN and self.numberOfMpi > 1:
                errors.append("Parameter estimation is not supported in MPI mode.")
        return errors

    def _getInputPath(self, *paths):
        return self._getPath('input', *paths)
