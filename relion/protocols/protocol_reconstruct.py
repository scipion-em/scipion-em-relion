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

from pyworkflow.protocol.params import (PointerParam, FloatParam,  
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pwem.objects import Volume
from pwem.protocols import ProtReconstruct3D
from pwem.constants import ALIGN_PROJ

import relion.convert as convert

# TODO: Check if we can centralize this, and how it combines with related functions


class ProtRelionReconstruct(ProtReconstruct3D):
    """ This protocol reconstructs a volume using Relion.

    Reconstruct a volume from a given set of particles.
    The alignment parameters will be converted to a Relion star file
    and used as direction projections to reconstruct.
    """
    _label = 'reconstruct'

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",
                      help='Select the input images from the project.')
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group",
                      help='See [[Relion Symmetry][http://www2.mrc-lmb.cam.ac.uk/'
                           'relion/index.php/Conventions_%26_File_formats#Symmetry]] '
                           'page for a description of the symmetry format '
                           'accepted by Relion')
        form.addParam('maxRes', FloatParam, default=-1,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).')
        form.addParam('pad', FloatParam, default=2,
                      label="Padding factor")
        
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Extra parameters: ', 
                      help='Extra parameters to *relion_reconstruct* program:\n'
                      """
                        --subtract ():\t\tSubtract projections of this map from the images used for reconstruction
                        --NN (false):\t\tUse nearest-neighbour instead of linear interpolation before gridding correction
                        --blob_r (1.9):\t\tRadius of blob for gridding interpolation
                        --blob_m (0):\t\tOrder of blob for gridding interpolation
                        --blob_a (15):\t\tAlpha-value of blob for gridding interpolation
                        --iter (10):\t\tNumber of gridding-correction iterations
                        --refdim (3):\t\tDimension of the reconstruction (2D or 3D)
                        --angular_error (0.):\t\tApply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles
                        --shift_error (0.):\t\tApply random deviations with this standard deviation (in pixels) to each of the 2 translations
                        --fom_weighting (false):\t\tWeight particles according to their figure-of-merit (_rlnParticleFigureOfMerit)
                        --fsc ():\t\tFSC-curve for regularized reconstruction
                        ...
                      """)
        form.addSection('CTF')
        form.addParam('doCTF', BooleanParam, default=False,
                      label='Apply CTF correction?')
        form.addParam('ctfIntactFirstPeak', BooleanParam, default=False,
                      condition='doCTF',
                      label='Leave CTFs intact until first peak?')
        
        form.addParallelSection(threads=1, mpi=1)
        # TODO: Add an option to allow the user to
        # decide if copy binary files or not
            
    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        # self._initialize()
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')
        self._insertReconstructStep()
        self._insertFunctionStep('createOutputStep')

    def _getProgram(self, program='relion_reconstruct'):
        """ Get the program name depending on the MPI use or not. """
        if self.numberOfMpi > 1:
            program += '_mpi'
        return program

    def _insertReconstructStep(self):
        imgSet = self.inputParticles.get()
        imgStar = self._getFileName('input_particles.star')
        
        params = ' --i %s' % imgStar
        params += ' --o %s' % self._getPath('output_volume.vol')
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --angpix %0.3f' % imgSet.getSamplingRate()
        params += ' --maxres %0.3f' % self.maxRes.get()
        params += ' --pad %0.3f' % self.pad.get()

        # TODO: Test that the CTF part is working
        if self.doCTF:
            params += ' --ctf'
            if self.ctfIntactFirstPeak:
                params += ' --ctf_intact_first_peak'

        if self.extraParams.hasValue():
            params += " " + self.extraParams.get()

        self._insertFunctionStep('reconstructStep', params)

    # -------------------------- STEPS functions ------------------------------
    def reconstructStep(self, params):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        self.runJob(self._getProgram(), params)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_particles.star': self._getTmpPath('input_particles.star'),
            'output_volume': self._getPath('output_volume.vol')
            }
        self._updateFilenamesDict(myDict)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self.inputParticles.get()
        imgStar = self._getFileName('input_particles.star')

        # Pass stack file as None to avoid write the images files
        convert.writeSetOfParticles(
            imgSet, imgStar, self._getTmpPath(), alignType=ALIGN_PROJ)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        volume = Volume()
        volume.setFileName(self._getFileName('output_volume'))
        volume.setSamplingRate(imgSet.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputParticles, volume)
    
    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        errors = []

        if self.numberOfMpi > 1 and self.numberOfThreads > 1:
            errors.append('Relion reconstruct can run either with mpi or '
                          'threads, not both!')
        return errors
    
    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        return []
