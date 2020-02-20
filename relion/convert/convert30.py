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
from io import open

from ..constants import *
from .convert_base import WriterBase, ReaderBase
from .convert_deprecated import rowToAlignment


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfImages subclasses
    with star file format previous to Relion>3.1, but providing the same
     interface as the new Writer class.
    """

    def writeSetOfMovies(self, moviesIterable, starFile):
        self._writeSetOfMoviesOrMics(moviesIterable, starFile,
                                     'movies', 'rlnMicrographMovieName')

    def writeSetOfMicrographs(self, micsIterable, starFile):
        self._writeSetOfMoviesOrMics(micsIterable, starFile,
                                     'micrographs', 'rlnMicrographName')

    def writeSetOfParticles(self, partsSet, starFile, **kwargs):
        # FIXME: Remove deprecated import
        from .convert_deprecated import _writeSetOfParticles
        _writeSetOfParticles(partsSet, starFile, self.outputDir, **kwargs)

    def _writeSetOfMoviesOrMics(self, imgIterable,
                                starFile, tableName, imgLabelName):
        """ Internal function to write either movies or micrographs
        star files. Input can be any iterable of these type of images (e.g
        set, list, etc).
        """
        # Process the first item and create the table based
        # on the generated columns
        self._imgLabelName = imgLabelName
        self._prefix = tableName[:3]

        micRow = OrderedDict()
        micRow[imgLabelName] = ''  # Just to add label, proper value later
        iterMics = iter(imgIterable)
        mic = next(iterMics)
        self._micToRow(mic, micRow)

        micsTable = self._createTableFromDict(micRow)

        while mic is not None:
            micRow[imgLabelName] = self._convert(mic)
            self._micToRow(mic, micRow)
            micsTable.addRow(**micRow)
            mic = next(iterMics, None)

        with open(starFile, 'w') as f:
            f.write("# Star file generated with Scipion\n")
            micsTable.writeStar(f, tableName=tableName)

    def _micToRow(self, mic, row):
        WriterBase._micToRow(self, mic, row)

        # Add now the Acquisition labels
        acq = mic.getAcquisition()
        row.update({
            'rlnVoltage': acq.getVoltage(),
            'rlnSphericalAberration': acq.getSphericalAberration(),
            'rlnAmplitudeContrast': acq.getAmplitudeContrast()
        })


class Reader(ReaderBase):
    def setParticleTransform(self, particle, row):
        """ Set the transform values from the row. """
        particle.setTransform(rowToAlignment(row, self._alignType))
