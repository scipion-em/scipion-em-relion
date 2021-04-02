# **************************************************************************
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
# **************************************************************************

from os.path import relpath

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pyworkflow.protocol import STEPS_SERIAL
from pwem.protocols import ProtParticlePickingAuto

from .protocol_autopick import ProtRelionAutopickBase


class ProtRelionAutopickLoG(ProtRelionAutopickBase):
    """
    This Relion protocol uses 'relion_autopick' program for the
    Laplacian of Gaussian (LoG) option.
    """
    _label = 'auto-picking LoG'
    _devStatus = PROD

    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    # -------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      label='Input micrographs', important=True,
                      help='Select the input micrographs. '
                           'If using the *Optimize* mode, just a subset of '
                           'micrographs are used to compute the FOM maps. '
                           'If in *Compute* mode, all micrographs will be '
                           'auto-picked.')

        form.addParam('boxSize', params.IntParam,
                      label='Box size (px)',
                      help="Box size in pixels.")

        group = form.addGroup('Laplacian of Gaussian')

        line = group.addLine('Diameter for LoG filter (A)',
                             help="Min and Max of LoG filter")

        line.addParam('minDiameter', params.IntParam, default=200,
                      label='Min')

        line.addParam('maxDiameter', params.IntParam, default=250,
                      label='Max')

        group.addParam('areParticlesWhite', params.BooleanParam,
                       default=False,
                       label='Are the particles white?',
                       help='Set this option to No if the particles are black, '
                            'and to Yes if the particles are white.')

        group.addParam('maxResolution', params.FloatParam, default=20,
                       label='Maximum resolution to consider (A)',
                       help='The Laplacian-of-Gaussian filter will be applied '
                            'to downscaled micrographs with the corresponding '
                            'size. Give a negative value to skip downscaling.')

        group.addParam('threshold', params.FloatParam, default=0,
                       label='Adjust default threshold (stddev)',
                       help='Use this to pick more (negative number -> lower '
                            'threshold) or less (positive number -> higher '
                            'threshold) particles compared to the default '
                            'setting.')

        group.addParam('threshold2', params.FloatParam, default=999,
                       label='Upper threshold (stddev)',
                       help='Use this to discard picks with LoG thresholds '
                            'that are this many standard deviations above '
                            'the average, e.g. to avoid high contrast '
                            'contamination like ice and ethane droplets. '
                            'Good values depend on the contrast of '
                            'micrographs and need to be interactively '
                            'explored; for low contrast micrographs, '
                            'values of ~ 1.5 may be reasonable, but the '
                            'same value will be too low for high-contrast '
                            'micrographs.')

        form.addParam('extraParams', params.StringParam, default='',
                      label='Additional arguments:',
                      help='In this box command-line arguments may be provided '
                           'that are not generated by the GUI. This may be '
                           'useful for testing developmental options and/or '
                           'expert use of the program. \n'
                           'The command "relion_autopick" will print a list '
                           'of possible options.')

        self._defineStreamingParams(form)

        form.addParallelSection(threads=0, mpi=4)

    # -------------------------- STEPS functions ------------------------------

    def getAutopickParams(self):
        """ Return the autopicking parameters except for the interactive ones. """
        params = ' --pickname autopick --odir "./" --LoG --shrink 0'
        params += ' --lowpass %0.3f' % self.maxResolution

        # Add extra params is any
        params += ' %s' % self.extraParams

        return params

    def _getPickArgs(self):
        """ This method is used by base class to provide the arguments to
        new inserted steps.
        """
        args = [
            self.getAutopickParams(),
            self.minDiameter.get(),
            self.maxDiameter.get(),
            self.threshold.get(),
            self.threshold2.get()
        ]

        return args

    def _pickMicrographsFromStar(self, micStarFile, cwd, params,
                                 minDiameter, maxDiameter, threshold, threshold2):
        """ Launch the 'relion_autopick' for micrographs in the inputStarFile.
         If the input set of complete, the star file will contain all the
         micrographs. If working in streaming, it will be only one micrograph.
        """
        params += ' --i %s' % relpath(micStarFile, cwd)
        params += ' --LoG_diam_min %0.3f' % minDiameter
        params += ' --LoG_diam_max %0.3f' % maxDiameter
        params += ' --LoG_adjust_threshold %0.3f' % threshold
        params += ' --LoG_upper_threshold %0.3f' % threshold2

        program = self._getProgram('relion_autopick')

        self.runJob(program, params, cwd=cwd)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        return summary

    # -------------------------- UTILS functions -------------------------------
    def getBoxSize(self):
        """ Return a reasonable box-size in pixels. """
        return self.boxSize.get()
