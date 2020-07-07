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
import relion
from os.path import abspath
from pwem.convert.transformations import translation_from_matrix
from relion.convert import Table
from .convert_base import WriterBase, ReaderBase


class Writer(WriterBase):
    """ Helper class to convert from Scipion SetOfImages subclasses
    with star file format previous to Relion>3.1, but providing the same
     interface as the new Writer class.
    """

    # def writeSetOfMovies(self, moviesIterable, starFile):
    #     self._writeSetOfMoviesOrMics(moviesIterable, starFile,
    #                                  'movies', 'rlnMicrographMovieName')
    #
    # def writeSetOfMicrographs(self, micsIterable, starFile):
    #     self._writeSetOfMoviesOrMics(micsIterable, starFile,
    #                                  'micrographs', 'rlnMicrographName')

    def writeSetOfSubtomograms(self, subtomoSet, subtomosStar, **kwargs):
        # # FIXME: Remove deprecated import
        # from .convert_deprecated import _writeSetOfParticles
        # _writeSetOfParticles(partsSet, starFile, **kwargs)
        tomoTable = self._createStarTomoTable()
        for subtomo in subtomoSet:
            angles, shifts = self._getTransformInfoFromSubtomo(subtomo)
            magn = subtomo.getAcquisition().getMagnification()
            rlnMicrographName = subtomo.getVolName()
            rlnCoordinateX = subtomo.getCoordinate3D().getX()
            rlnCoordinateY = subtomo.getCoordinate3D().getY()
            rlnCoordinateZ = subtomo.getCoordinate3D().getZ()
            rlnImageName = subtomo.getFileName()
            rlnCtfImage = abspath(self._getCTFFileFromSubtomo(subtomo))
            rlnMagnification = magn if magn else 10000  # 64000
            rlnDetectorPixelSize = subtomo.getSamplingRate()
            rlnAngleRot = angles[0]
            rlnAngleTilt = angles[1]
            rlnAnglePsi = angles[2]
            rlnOriginX = shifts[0]
            rlnOriginY = shifts[1]
            rlnOriginZ = shifts[2]
            # Add row to the table which will be used to generate the STAR file
            tomoTable.addRow(rlnMicrographName,
                             rlnCoordinateX,
                             rlnCoordinateY,
                             rlnCoordinateZ,
                             rlnImageName,
                             rlnCtfImage,
                             rlnMagnification,
                             rlnDetectorPixelSize,
                             rlnAngleRot,
                             rlnAngleTilt,
                             rlnAnglePsi,
                             rlnOriginX,
                             rlnOriginY,
                             rlnOriginZ
                             )
        # Write the STAR file
        if relion.Plugin.IS_30():
            tomoTable.write(subtomosStar)

    @ staticmethod
    def _createStarTomoTable():
        return Table(columns=['rlnMicrographName',
                              'rlnCoordinateX',
                              'rlnCoordinateY',
                              'rlnCoordinateZ',
                              'rlnImageName',
                              'rlnCtfImage',
                              'rlnMagnification',
                              'rlnDetectorPixelSize',
                              'rlnAngleRot',
                              'rlnAngleTilt',
                              'rlnAnglePsi',
                              'rlnOriginX',
                              'rlnOriginY',
                              'rlnOriginZ',
                              ])

    @ staticmethod
    def _getCTFFileFromSubtomo(subtomo):
        return subtomo.getCoordinate3D()._3dcftMrcFile.get()


    @staticmethod
    def _getTransformInfoFromSubtomo(subtomo):
        angles = [0, 0, 0]
        shifts = [0, 0, 0]
        T = subtomo.getTransform()
        if T:  # Alignment performed before

            # TODO: check if matrix must be inverted to get the correct angles
            M = subtomo.getTransform().getMatrix()

            from relion.convert import geometryFromMatrix
            calcInv = False
            _, angles = geometryFromMatrix(M, calcInv)
            shifts = translation_from_matrix(M)
            if calcInv:
                shifts = -shifts

            # # Direct
            # angles = -rad2deg(euler_from_matrix(M, axes='szyz'))
            # shifts = translation_from_matrix(M)

            # # Inverse
            # angularPart = array(M, dtype=float64, copy=False)[:3, :3]
            # shifts = -translation_from_matrix(M)
            # angularPart = linalg.inv(angularPart)
            # angles = -rad2deg(euler_from_matrix(angularPart, axes='szyz'))
        return angles, shifts



#
# class Reader(ReaderBase):
#
#     def readSetOfParticles(self, starFile, partsSet, **kwargs):
#         """ Convert a star file into a set of particles.
#
#         Params:
#             starFile: the filename of the star file
#             partsSet: output particles set
#
#         Keyword Arguments:
#             blockName: The name of the data block (default particles)
#             alignType:
#             removeDisabled:
#
#         """
#         readSetOfParticles(starFile, partsSet, **kwargs)
#
#     def setParticleTransform(self, particle, row):
#         """ Set the transform values from the row. """
#         particle.setTransform(rowToAlignment(row, self._alignType))
