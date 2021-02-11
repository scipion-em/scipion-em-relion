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
"""
New conversion functions dealing with Relion3.1 new star files format.
"""

import os
from emtable import Table

import pyworkflow.utils as pwutils
from pyworkflow.object import ObjectWrap

from pwem.constants import ALIGN_NONE
from pwem.emlib.image import ImageHandler


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
        self._optics = kwargs.get('optics', None)
        # Not used now
        # self.convertPolicy = kwargs.get('convertPolicy', self.CONVERT_IF_NEEDED)
        self.rootDir = None
        self.outputDir = None
        self.outputStack = None
        self.useBaseName = False
        self.extensions = ['mrc']
        self.update(['rootDir', 'outputDir', 'userBaseName', 'extensions'],
                    **kwargs)
        self._ih = ImageHandler()  # used to convert images
        self._filesDict = {}  # used to map file names (converted or linked)
        self._dimensionality = 2
        self._imageSize = None

    def update(self, attrsList, **kwargs):
        """ Update the some attributes with values from kwargs. """
        for attr in attrsList:
            if attr in kwargs:
                setattr(self, attr, kwargs[attr])

    def writeSetOfMovies(self, moviesIterable, starFile, **kwargs):
        pass

    def writeSetOfMicrographs(self, micsIterable, starFile, **kwargs):
        pass

    def _writeSetOfMoviesOrMics(self, imgIterable,
                                starFile, tableName, imgLabelName, **kwargs):
        pass

    def writeSetOfParticles(self, partsSet, starFile, **kwargs):
        """ Convert a set of particles into a star file and maybe binary files.

        Params:
            partsSet: input particles set
            starFile: the filename of the star file that will be written.

        Keyword Arguments:
            blockName: The name of the data block (default particles)
            fillMagnification: If True set magnification values (default False)

            outputStack: A file name to write all particles. If this option is
                passed, then the outputDir will be ignored.

        """
        pass

    def _convert(self, image):
        imageFn = image.getFileName()

        if self.outputDir is None:
            return imageFn

        ext = pwutils.getExt(imageFn)[1:]
        if ext in self.extensions:
            finalExt = ext
            convertFunc = pwutils.createAbsLink
        else:
            finalExt = self.extensions[0]
            convertFunc = self._ih.convert

        if self.useBaseName:
            newName = pwutils.replaceBaseExt(image.getFileName(), finalExt)
        else:
            newName = "%s_%06d.%s" % (self._prefix, image.getObjId(), finalExt)

        newPath = os.path.join(self.outputDir, newName)
        convertFunc(os.path.abspath(imageFn), newPath)
        # If there is a rootDir defined, we should return the path relative
        # to that location, probably to be used from the star file
        if self.rootDir is not None:
            newPath = os.path.relpath(newPath, self.rootDir)

        return newPath

    def _createTableFromDict(self, rowDict):
        """ Helper function to create a Table instance from
        an input dict with keys as columns names and type
        the type of the values in the dict.
        """
        return Table(columns=[
            Table.Column(k, type=type(v)) for k, v in rowDict.items()])

    def _micToRow(self, mic, row):
        row['rlnImageId'] = mic.getObjId()

    def _ctfToRow(self, ctf, row):
        psd = ctf.getPsdFile()
        if psd:
            row['rlnCtfImage'] = psd
        dU, dV, dAngle = ctf.getDefocus()
        row['rlnDefocusU'] = dU
        row['rlnDefocusV'] = dV
        row['rlnCtfAstigmatism'] = abs(dU-dV)
        row['rlnDefocusAngle'] = dAngle
        row['rlnCtfFigureOfMerit'] = ctf.getFitQuality() or 0
        row['rlnCtfMaxResolution'] = ctf.getResolution() or 0

        phaseShift = ctf.getPhaseShift()

        if phaseShift is not None:
            row['rlnPhaseShift'] = phaseShift


class ReaderBase:
    """ Helper class to grab information from star file rows
     and fill the required values in Scipion objects
     (e.g particles, micrographs, etc)
    """
    def __init__(self, **kwargs):
        """
        """
        self._alignType = kwargs.get('alignType', ALIGN_NONE)
        self._pixelSize = kwargs.get('pixelSize', 1.0)
        self._invPixelSize = 1. / self._pixelSize
        self._extraLabels = []

    def readSetOfParticles(self, starFile, partsSet, **kwargs):
        """ Convert a star file into a set of particles.

        Params:
            starFile: the filename of the star file
            partsSet: output particles set

        Keyword Arguments:
            blockName: The name of the data block (default particles)
            alignType:
            removeDisabled:

        """
        pass

    def setParticleTransform(self, particle, row):
        """ Set the transform values from the row. """
        pass

    def createExtraLabels(self, item, row, extraLabels):
        """ Create new Objects for each extra label if contained in
        the columnObj. It will set the self._extraLabels property.
        Args:
            item: Object item that will have new extra labels objects
            row: column object that should have a method hasColumn
            extraLabels: list of label names that will be set if present
                in columnObj
        """
        self._extraLabels = [l for l in extraLabels if row.hasColumn(l)]
        for label in self._extraLabels:
            setattr(item, '_' + label,
                    ObjectWrap(getattr(row, label)))

    def setExtraLabels(self, item, row):
        """ Set values for already computed extraLabels with
        self.createExtraLabels. """
        for label in self._extraLabels:
            getattr(item, '_' + label).set(getattr(row, label))
