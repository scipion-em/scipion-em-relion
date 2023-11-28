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

import os

import pyworkflow.utils as pwutils
from pwem.protocols import ProtParticlePickingAuto

import relion.convert as convert
from .protocol_base import ProtRelionBase


class ProtRelionAutopickBase(ProtParticlePickingAuto, ProtRelionBase):
    """ Base class for auto-picking protocols in Relion.
    """
    _label = None

    def _pickMicrograph(self, mic, *args):
        """ This method should be invoked only when working in streaming mode.
        """
        self._pickMicrographList([mic], *args)

    def _pickMicrographList(self, micList, *args):
        if not micList:
            return

        micsDir = self._createTmpMicsDir(micList)
        micStar = os.path.join(micsDir, 'input_micrographs.star')
        writer = convert.createWriter(rootDir=micsDir, outputDir=micsDir)
        writer.writeSetOfMicrographs(micList, micStar)
        self._pickMicrographsFromStar(micStar, micsDir, *args)
        # Move coordinates files to tmp
        os.system('mv %s/*autopick.star %s/' % (micsDir, self._getTmpPath()))

    def _createSetOfCoordinates(self, micSet, suffix=''):
        """ Override this method to set the box size. """
        coordSet = ProtParticlePickingAuto._createSetOfCoordinates(
            self, micSet, suffix=suffix)
        coordSet.setBoxSize(self.getBoxSize())
        return coordSet

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        """ Parse back the output star files and populate the SetOfCoordinates.
        """
        template = self._getTmpPath("mic_%06d_autopick.star")
        starFiles = [template % mic.getObjId() for mic in micList]
        convert.readSetOfCoordinates(coordSet, starFiles, micList)

    def _pickMicrographsFromStar(self, micStar, micsDir, *args):
        """ Should be defined in subclasses. """
        pass

    def getBoxSize(self):
        """ Return a reasonable box-size in pixels. """
        return None

    def getInputMicrographsPointer(self):
        return self.inputMicrographs

    def getInputMicrographs(self):
        return self.getInputMicrographsPointer().get()

    def getMicrographList(self):
        """ Return the list of micrographs
        that will be used for optimizing the parameters or the picking.
        """
        return self.getInputMicrographs()

    def getCoordsDir(self):
        return self._getTmpPath('xmipp_coordinates')

    def _writeXmippCoords(self, coordSet):
        micSet = self.getInputMicrographs()
        coordPath = self._getTmpPath('xmipp_coordinates')
        pwutils.cleanPath(coordPath)
        pwutils.makePath(coordPath)
        micPath = micSet.getFileName()
        convert.writeSetOfCoordinatesXmipp(coordPath, coordSet, ismanual=False)
        return micPath, coordPath

    def writeXmippOutputCoords(self):
        return self._writeXmippCoords(self.outputCoordinates)

    def writeXmippCoords(self):
        """ Write the SetOfCoordinates as expected by Xmipp
        to be displayed with its GUI.
        """
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(self.getBoxSize())
        starFiles = [self._getExtraPath(pwutils.removeBaseExt(mic.getFileName())
                                        + '_autopick.star') for mic in micSet]
        convert.readSetOfCoordinates(coordSet, starFiles)
        return self._writeXmippCoords(coordSet)

    def __getMicListPrefix(self, micList):
        n = len(micList)
        if n == 0:
            raise ValueError("Empty micrographs list!")
        micsPrefix = 'mic_%06d' % micList[0].getObjId()
        if n > 1:
            micsPrefix += "-%06d" % micList[-1].getObjId()
        return micsPrefix

    def _getMicStarFile(self, micList):
        return self._getTmpPath(self.__getMicListPrefix(micList) + '.star')

    def _createTmpMicsDir(self, micList):
        """ Create a temporary path to work with a list of micrographs. """
        micsDir = self._getTmpPath(self.__getMicListPrefix(micList))
        pwutils.makePath(micsDir)
        return micsDir
