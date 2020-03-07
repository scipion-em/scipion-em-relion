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
from io import open
import numpy as np
from collections import OrderedDict

import pwem
import pwem.convert.transformations as tfs

from .convert_base import WriterBase, ReaderBase
from .convert_utils import convertBinaryFiles, locationToRelion


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfImages subclasses
    into Relion>3.1 star files (and binaries if conversion needed).
    """
    def writeSetOfMovies(self, moviesIterable, starFile):
        self._writeSetOfMoviesOrMics(moviesIterable, starFile,
                                     'movies', 'rlnMicrographMovieName')

    def writeSetOfMicrographs(self, micsIterable, starFile):
        self._writeSetOfMoviesOrMics(micsIterable, starFile,
                                     'micrographs', 'rlnMicrographName')

    def _writeSetOfMoviesOrMics(self, imgIterable,
                                starFile, tableName, imgLabelName):
        """ This function can be used to write either movies or micrographs
        star files. Input can be any iterable of these type of images (e.g
        set, list, etc).
        """
        # Process the first item and create the table based
        # on the generated columns
        self._imgLabelName = imgLabelName
        self._imgLabelPixelSize = 'rlnMicrographPixelSize'

        self._prefix = tableName[:3]
        self._optics = OrderedDict()
        micRow = OrderedDict()
        micRow[imgLabelName] = ''  # Just to add label, proper value later
        iterMics = iter(imgIterable)
        mic = next(iterMics)
        self._imageSize = mic.getXDim()
        self._micToRow(mic, micRow)

        opticsTable = self._createTableFromDict(list(self._optics.values())[0])
        micsTable = self._createTableFromDict(micRow)

        while mic is not None:
            micRow[imgLabelName] = self._convert(mic)
            self._micToRow(mic, micRow)
            micsTable.addRow(**micRow)
            mic = next(iterMics, None)

        for opticsDict in self._optics.values():
            opticsTable.addRow(**opticsDict)

        with open(starFile, 'w') as f:
            f.write("# Star file generated with Scipion\n")
            f.write("# version 30001\n")
            opticsTable.writeStar(f, tableName='optics')
            f.write("# version 30001\n")
            micsTable.writeStar(f, tableName=tableName)

    def _getOpticsGroupNumber(self, img):
        """ Get the optics group number based on acquisition.
        Params:
            img: input image, movie, micrograph or particle
        """
        # Add now the new Optics Group stuff
        acq = img.getAcquisition()
        ogName = acq.opticsGroupName.get() or 'DefaultOpticsGroup'
        ps = img.getSamplingRate()

        if ogName not in self._optics:
            ogNumber = len(self._optics) + 1
            self._optics[ogName] = {
                'rlnOpticsGroupName': ogName,
                'rlnOpticsGroup': ogNumber,
                #'rlnMtfFileName': acq.mtfFile.get() or 'No-MTF',
                # FIXME: Check when we need to update the following
                'rlnMicrographOriginalPixelSize': ps,
                self._imgLabelPixelSize: ps,
                'rlnVoltage': acq.getVoltage(),
                'rlnSphericalAberration': acq.getSphericalAberration(),
                'rlnAmplitudeContrast': acq.getAmplitudeContrast(),
                'rlnBeamTiltX': acq.beamTiltX.get() or 0.,
                'rlnBeamTiltY': acq.beamTiltY.get() or 0.,
                'rlnImageDimensionality': self._dimensionality,
                'rlnImageSize': self._imageSize,
            }
            mtfFile = acq.mtfFile.get()
            if mtfFile is not None:
                self._optics[ogName]['rlnMtfFileName'] = mtfFile
        else:
            ogNumber = self._optics[ogName]['rlnOpticsGroup']

        return ogNumber

    def _setAttributes(self, obj, row, attributes):
        for attr in attributes:
            attrLabel = '_%s' % attributes
            if hasattr(obj, attrLabel):
                row[attr] = obj.getAttributeValue(attrLabel)

    def _micToRow(self, mic, row):
        WriterBase._micToRow(self, mic, row)
        row['rlnOpticsGroup'] = self._getOpticsGroupNumber(mic)

    def _align2DToRow(self, alignment, row):
        matrix = alignment.getMatrix()
        shifts = tfs.translation_from_matrix(matrix)
        shifts *= self._pixelSize
        angles = -np.rad2deg(tfs.euler_from_matrix(matrix, axes='szyz'))
        row['rlnOriginXAngst'], row['rlnOriginYAngst'] = shifts[:2]
        row['rlnAnglePsi'] = -(angles[0] + angles[2])

    def _alignProjToRow(self, alignment, row):
        matrix = np.linalg.inv(alignment.getMatrix())
        shifts = -tfs.translation_from_matrix(matrix)
        shifts *= self._pixelSize
        angles = -np.rad2deg(tfs.euler_from_matrix(matrix, axes='szyz'))
        row['rlnOriginXAngst'], row['rlnOriginYAngst'], row['rlnOriginZAngst'] = shifts
        row['rlnAngleRot'], row['rlnAngleTilt'], row['rlnAnglePsi'] = angles

    def _partToRow(self, part, row):
        row['rlnImageId'] = part.getObjId()

        # Add coordinate information
        coord = part.getCoordinate()
        if coord is not None:
            x, y = coord.getPosition()
            row['rlnCoordinateX'] = x
            row['rlnCoordinateY'] = y
            # Add some specify coordinate attributes
            self._setAttributes(coord, row, ['rlnClassNumber',
                                             'rlnAutopickFigureOfMerit',
                                             'rlnAnglePsi'])
            micName = coord.getMicName()
            if micName:
                row['rlnMicrographName'] = str(micName.replace(" ", ""))
            else:
                if coord.getMicId():
                    row['rlnMicrographName'] = str(coord.getMicId())

        index, fn = part.getLocation()
        if self.outputStack:
            row['rlnOriginalParticleName'] = locationToRelion(index, fn)
            index, fn = self._counter, self._relOutputStack
            if self._counter > 0:
                self._ih.convert(part, (index, self.outputStack))
        else:
            if self.outputDir is not None:
                fn = self._filesDict.get(fn, fn)

        row['rlnImageName'] = locationToRelion(index, fn)

        if self._setRandomSubset:
            row['rlnRandomSubset'] = part._rlnRandomSubset.get()

        # Set CTF values
        if self._setCtf:
            self._ctfToRow(part.getCTF(), row)

        # Set alignment if necessary
        if self._setAlign:
            self._setAlign(part.getTransform(), row)

        # Set additional labels if present
        self._setAttributes(part, row, self._extraLabels)

        # Add now the new Optics Group stuff
        row['rlnOpticsGroup'] = self._getOpticsGroupNumber(part)

        self._counter += 1

    def writeSetOfParticles(self, partsSet, starFile, **kwargs):
        # Process the first item and create the table based
        # on the generated columns
        self._imgLabelPixelSize = 'rlnImagePixelSize'

        self._optics = OrderedDict()
        partRow = OrderedDict()
        firstPart = partsSet.getFirstItem()

        # Convert binaries if required
        self.outputStack = kwargs.get('outputStack', None)
        if self.outputStack:
            self._relOutputStack = os.path.relpath(self.outputStack,
                                                   os.path.dirname(starFile))
        if self.outputDir is not None:
            forceConvert = kwargs.get('forceConvert', False)
            self._filesDict = convertBinaryFiles(partsSet, self.outputDir,
                                                 forceConvert=forceConvert)

        # Compute some flags from the first particle...
        # when flags are True, some operations will be applied to all particles
        self._preprocessImageRow = kwargs.get('preprocessImageRow', None)
        self._setRandomSubset = (kwargs.get('fillRandomSubset') and
                                 firstPart.hasAttribute('_rlnRandomSubset'))

        self._setCtf = kwargs.get('writeCtf', True) and firstPart.hasCTF()

        alignType = kwargs.get('alignType', partsSet.getAlignment())

        if alignType == pwem.ALIGN_2D:
            self._setAlign = self._align2DToRow
        elif alignType == pwem.ALIGN_PROJ:
            self._setAlign = self._alignProjToRow
        elif alignType == pwem.ALIGN_3D:
            raise Exception(
                "3D alignment conversion for Relion not implemented. "
                "It seems the particles were generated with an incorrect "
                "alignment type. You may either re-launch the protocol that "
                "generates the particles with angles or set 'Consider previous"
                " alignment?' to No")
        elif alignType == pwem.ALIGN_NONE:
            self._setAlign = None
        else:
            raise Exception("Invalid value for alignType: %s" % alignType)

        self._extraLabels = kwargs.get('extraLabels', [])
        self._extraLabels.extend(['rlnParticleSelectZScore',
                                  'rlnMovieFrameNumber'])
        self._postprocessImageRow = kwargs.get('postprocessImageRow', None)

        self._imageSize = firstPart.getXDim()
        self._pixelSize = firstPart.getSamplingRate() or 1.0

        self._counter = 0  # Mark first conversion as special one
        self._partToRow(firstPart, partRow)
        opticsTable = self._createTableFromDict(list(self._optics.values())[0])
        partsTable = self._createTableFromDict(partRow)
        partsTable.addRow(**partRow)

        with open(starFile, 'w') as f:
            # Write particles table
            f.write("# Star file generated with Scipion\n")
            f.write("# version 30001\n")
            # Write header first
            partsTable.writeStar(f, tableName='particles', writeRows=False)
            # Write all rows
            for part in partsSet:
                self._partToRow(part, partRow)
                partsTable.writeStarLine(f, partRow.values())

            # Write Optics at the end
            for opticsDict in self._optics.values():
                opticsTable.addRow(**opticsDict)
            f.write("\n# version 30001\n")
            opticsTable.writeStar(f, tableName='optics')


class Reader(ReaderBase):

    ALIGNMENT_LABELS = [
        "rlnOriginXAngst",
        "rlnOriginYAngst",
        "rlnOriginZAngst",
        "rlnAngleRot",
        "rlnAngleTilt",
        "rlnAnglePsi",
    ]

    def __init__(self, **kwargs):
        """
        """
        ReaderBase.__init__(self, **kwargs)
        self._first = False

    def setParticleTransform(self, particle, row):
        """ Set the transform values from the row. """
        self._pixelSize = particle.getSamplingRate()
        self._invPixelSize = 1. / self._pixelSize

        if (self._alignType == pwem.ALIGN_NONE or
            not row.containsAny(self.ALIGNMENT_LABELS)):
            self.setParticleTransform = self.__setParticleTransformNone
        else:
            # Ensure the Transform object exists
            self._angles = np.zeros(3)
            self._shifts = np.zeros(3)

            particle.setTransform(pwem.objects.Transform())

            if self._alignType == pwem.ALIGN_2D:
                self.setParticleTransform = self.__setParticleTransform2D
            elif self._alignType == pwem.ALIGN_PROJ:
                self.setParticleTransform = self.__setParticleTransformProj
            else:
                raise Exception("Unexpected alignment type: %s"
                                % self._alignType)

        # Call again the modified function
        self.setParticleTransform(particle, row)

    def __setParticleTransformNone(self, particle, row):
        particle.setTransform(None)

    def __setParticleTransform2D(self, particle, row):
        angles = self._angles
        shifts = self._shifts
        ips = self._invPixelSize

        def _get(label):
            return float(row.getValue(label, 0.))

        shifts[0] = _get('rlnOriginXAngst') * ips
        shifts[1] = _get('rlnOriginYAngst') * ips
        angles[2] = _get('rlnAnglePsi')
        radAngles = -np.rad2deg(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = shifts[:3]
        particle.getTransform().setMatrix(M)

    def __setParticleTransformProj(self, particle, row):
        angles = self._angles
        shifts = self._shifts
        ips = self._invPixelSize

        def _get(label):
            return float(row.getValue(label, 0.))

        shifts[0] = _get('rlnOriginXAngst') * ips
        shifts[1] = _get('rlnOriginYAngst') * ips
        shifts[2] = _get('rlnOriginZAngst') * ips

        angles[0] = _get('rlnAngleRot')
        angles[1] = _get('rlnAngleTilt')
        angles[2] = _get('rlnAnglePsi')

        radAngles = -np.deg2rad(angles)

        # TODO: jmrt: Maybe we should test performance and consider if keeping
        # TODO: the matrix and not creating one everytime will make things faster
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = -shifts[:3]
        M = np.linalg.inv(M)
        particle.getTransform().setMatrix(M)
