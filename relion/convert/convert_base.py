# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
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
"""
New conversion functions dealing with Relion3.1 new star files format.
"""

import os

import pyworkflow.em as pwem
import pyworkflow.utils as pwutils

from ..constants import *
from .metadata import Table, Column


class WriterBase:
    """ Helper class to convert from Scipion SetOfImages subclasses
    into Relion>3.1 star files (and binaries if conversion needed).
    """
    def __init__(self, **kwargs):
        """
        Create a new instance with some configuration parameters.

        Keyword Args:
            rootDir: Specify a directory that will be used as "root"
                for setting the path values in the star file pointing
                to the binaries
            convertPolicy: By default, conversion of binary files will
                create symbolic links if no conversion is required.
                By passing convertPolicy=CONVERT_ALWAYS, it will force
                the conversion.
            useBaseName: (bool) By default the writer will use the id to
                generate shorter names. If this option is True, then
                the images base name will be used instead. This option
                might be useful in export protocols.

        """
        self._optics = None
        # Not used now
        # self.convertPolicy = kwargs.get('convertPolicy', CONVERT_IF_NEEDED)
        self.rootDir = kwargs.get('rootDir', None)
        self.outputDir = kwargs.get('outputDir', None)
        self.useBaseName = kwargs.get('useBaseName', False)
        self.extensions = kwargs.get('extensions', ['mrc'])
        self._ih = pwem.convert.ImageHandler()  # used to convert images

    def _convert(self, image):
        imageFn = image.getFileName()

        if self.outputDir is None:
            return imageFn

        ext = pwutils.getExt(imageFn)

        if ext in self.extensions:
            finalExt = ext
            convertFunc = pwutils.createLink
        else:
            finalExt = self.extensions[0]
            convertFunc = self._ih.convert

        if self.useBaseName:
            newName = pwutils.replaceBaseExt(image.getFileName(), finalExt)
        else:
            newName = "%s_%06d.%s" % (self._prefix, image.getObjId(), finalExt)

        newPath = os.path.join(self.outputDir, newName)
        convertFunc(imageFn, newPath)
        # If there is a rootDir defined, we should return the path relative
        # to that location, probably to be used from the star file
        if self.rootDir is not None:
            newPath = os.path.relpath(newPath, self.rootDir)

        return newPath

    def writeSetOfMovies(self, moviesIterable, starFile):
        pass

    def writeSetOfMicrographs(self, micsIterable, starFile):
        pass

    def _writeSetOfMoviesOrMics(self, imgIterable,
                                starFile, tableName, imgLabelName):
        pass

    def _createTableFromDict(self, rowDict):
        """ Helper function to create a Table instance from
        an input dict with keys as columns names and type
        the type of the values in the dict.
        """
        return Table(columns=[
            Column(k, type=type(v)) for k, v in rowDict.iteritems()])

    def _micToRow(self, mic, row):
        row['rlnImageId'] = mic.getObjId()

        if mic.hasCTF():
            self._ctfToRow(mic.getCTF(), row)

    def _ctfToRow(self, ctf, row):
        row['rlnCtfImage'] = ctf.getPsdFile()
        dU, dV, dAngle = ctf.getDefocus()
        row['rlnDefocusU'] = dU
        row['rlnDefocusV'] = dV
        # FIXME Check how astigmatism is defined in Relion
        row['rlnCtfAstigmatism'] = dU / dV
        row['rlnDefocusAngle'] = dAngle
        row['rlnCtfFigureOfMerit'] = ctf.getFitQuality()
        row['rlnCtfMaxResolution'] = ctf.getResolution()

        phaseShift = ctf.getPhaseShift()

        if phaseShift is not None:
            row['rlnCtfPhaseShift'] = phaseShift

