# **************************************************************************
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
# **************************************************************************

import pyworkflow.em as em
import pyworkflow.protocol.params as params

from .protocol_base import ProtRelionBase


class ProtRelionAssignOpticsGroup(ProtRelionBase):
    """ Assign Optics Group name and related parameters to an input set.
     Input set can be: movies, micrographs or particles.
    """
    _label = 'assign optics group'
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam,
                      pointerClass='SetOfMovies,SetOfMicrographs,SetOfParticles',
                      label="Input set", important=True,
                      help='Select the input set (Movies, Micrographs or '
                           'Particles) to assign Optics Group parameters.')

        form.addParam('opticsGroupName', params.StringParam,
                      default='opticsGroup1',
                      label='Optics group name',
                      help='Relion-specific option. Name of this optics group. '
                           'Each group of movies with different '
                           'optics characteristics for CTF refinement '
                           'should have a unique name.')

        form.addParam('mtfFile', params.FileParam,
                      label='MTF-curve file',
                      help='User-provided STAR-file with the MTF-curve '
                           'of the detector. Use the wizard to load one '
                           'of the predefined ones provided at:\n'
                           '- [[https://www3.mrc-lmb.cam.ac.uk/relion/index.php/'
                           'FAQs#Where_can_I_find_MTF_curves_for_typical_detectors.3F]'
                           '[Relion\'s Wiki FAQs]]\n'
                           ' - [[http://www.gatan.com/K3][Gatan\'s website]]\n\n'
                           'Relion param: *--mtf*')

        line = form.addLine('Beam tilt (mrad)',
                            help='Known beam tilt in the X/Y-direction (in mrad). '
                                 'Set to zero if unknown.')
        line.addParam('beamTiltX', params.FloatParam, default=0.,
                      label='X')
        line.addParam('beamTiltY', params.FloatParam, default=0.,
                      label='Y')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',
                                 self.inputSet.get().getObjId())

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self, inputId):
        inputSet = self.inputSet.get()
        
        if isinstance(inputSet, em.SetOfMovies):
            outputSet = self._createSetOfMovies()
            outputName = 'outputMovies'
        elif isinstance(inputSet, em.SetOfMicrographs):
            outputSet = self._createSetOfMicrographs()
            outputName = 'outputMicrographs'
        elif isinstance(inputSet, em.SetOfParticles):
            outputSet = self._createSetOfParticles()
            outputName = 'outputParticles'
        else:
            raise Exception("Invalid input of type %s, expecting:\n"
                            "SetOfMovies, SetOfMicrographs or SetOfParticles"
                            % inputSet.getClassName())

        # Copy general info from input set
        outputSet.copyInfo(inputSet)
        # Update the acquisition object with new parameters from input
        acq = outputSet.getAcquisition()
        for attr in ['opticsGroupName', 'mtfFile', 'beamTiltX', 'beamTiltY']:
            setattr(acq, attr, getattr(self, attr).clone())

        # Copy items from input and set the new Optics Group name
        self._outputAcquisition = acq
        outputSet.copyItems(inputSet, updateItemCallback=self._setOpticsGroupName)
        self._defineOutputs(**{outputName: outputSet})
        self._defineTransformRelation(inputSet, outputSet)
    
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        validateMsgs = []
        return validateMsgs
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return []

    # --------------------------- UTILS functions -----------------------------
    def _setOpticsGroupName(self, item, row):
        # Take advantage of the fact that the 'item' object is the
        # used during the iteration in the copyItems, this can be a wrong
        # assumtion of the implementation changes, but I will take the risk
        # now for the sake of performance and only set the value once
        if getattr(self, '__firstTime', True):
            item.setAcquisition(self._outputAcquisition)
            self.__firstTime = False

