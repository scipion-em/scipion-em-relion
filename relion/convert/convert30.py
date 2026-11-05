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

import os
import numpy as np
from emtable import Table

import pwem
from pwem.constants import ALIGN_NONE
from pwem.objects import Particle, CTFModel, Coordinate, Acquisition
import pwem.convert.transformations as tfs

from relion.constants import PARTICLE_EXTRA_LABELS, LABELS_DICT
from .convert_base import ReaderBase
from .convert_utils import relionToLocation


class Reader(ReaderBase):

    ALIGNMENT_LABELS = [
        "rlnOriginX",
        "rlnOriginY",
        "rlnOriginZ",
        "rlnAngleRot",
        "rlnAngleTilt",
        "rlnAnglePsi",
    ]

    CTF_LABELS = [
        "rlnDefocusU",
        "rlnDefocusV",
        "rlnDefocusAngle",
    ]

    COORD_LABELS = [
        "rlnCoordinateX",
        "rlnCoordinateY",
        # extra labels below
        "rlnAutopickFigureOfMerit",
        "rlnClassNumber",
        "rlnAnglePsi"
    ]

    # read by the old xmipp-based reader through IMAGE_EXTRA_LABELS
    EXTRA_LABELS = [
        "rlnParticleSelectZScore",
        "rlnMovieFrameNumber",
        "rlnParticleId",
    ]

    def readSetOfParticles(self, starFile, partsSet, **kwargs):
        """ Convert a star file into a set of particles.

        Params:
            starFile: the filename of the star file
            partsSet: output particles set

        Keyword Arguments:
            alignType: alignment type
            removeDisabled: drop the rows marked as disabled (default True)

        """
        self._preprocessImageRow = kwargs.get('preprocessImageRow', None)
        self._postprocessImageRow = kwargs.get('postprocessImageRow', None)
        self._alignType = kwargs.get('alignType', ALIGN_NONE)
        self._pixelSize = kwargs.get('pixelSize', 1.0)
        self._invPixelSize = 1. / self._pixelSize

        # 3.0 star files have a single table, and its name is anything:
        # relion wrote data_images, Scipion data_Particles and the
        # tutorial data_noname, so take the first one whatever it is,
        # which is what md.MetaData(fn) used to do here
        partsReader = Table.Reader(starFile, types=LABELS_DICT)
        rows = iter(partsReader)

        if kwargs.get('removeDisabled', True):
            rows = (row for row in rows
                    if getattr(row, 'rlnEnabled', 1) != 0)

        firstRow = next(rows, None)
        if firstRow is None:
            raise ValueError("No particles found in: %s" % starFile)

        self._setClassId = firstRow.hasColumn('rlnClassNumber')
        self._setCtf = all(firstRow.hasColumn(l) for l in self.CTF_LABELS)
        self._setCoord = all(firstRow.hasColumn(l)
                             for l in self.COORD_LABELS[:2])
        self._setMicId = firstRow.hasColumn('rlnMicrographId')

        particle = Particle()
        if self._setCtf:
            particle.setCTF(CTFModel())

        if kwargs.get('readAcquisition', True):
            acq = Acquisition()
            self.rowToAcquisition(firstRow, acq)
            acq.setMagnification(kwargs.get('magnification', 10000))
            partsSet.setAcquisition(acq)

        extraLabels = (kwargs.get('extraLabels', []) +
                       PARTICLE_EXTRA_LABELS + self.EXTRA_LABELS)
        self.createExtraLabels(particle, firstRow, extraLabels)
        self._partExtraLabels = list(self._extraLabels)

        coordExtraLabels = (kwargs.get('coordExtraLabels', []) +
                            self.COORD_LABELS[2:])
        self._coordExtraLabels = [l for l in dict.fromkeys(coordExtraLabels)
                                  if firstRow.hasColumn(l)]

        self._rowToPart(firstRow, particle)
        partsSet.append(particle)

        for row in rows:
            self._rowToPart(row, particle)
            partsSet.append(particle)

        partsSet.setHasCTF(self._setCtf)
        partsSet.setAlignment(self._alignType)

    def _rowToPart(self, row, particle):
        particle.setObjId(getattr(row, 'rlnImageId', None))

        if self._preprocessImageRow:
            self._preprocessImageRow(particle, row)

        index, filename = relionToLocation(row.rlnImageName)
        particle.setLocation(index, filename)

        if self._setClassId:
            particle.setClassId(row.rlnClassNumber)

        if self._setCtf:
            self.rowToCtf(row, particle.getCTF())

        if self._setMicId:
            particle.setMicId(row.rlnMicrographId)

        self.setParticleTransform(particle, row)
        self._extraLabels = self._partExtraLabels
        self.setExtraLabels(particle, row)

        if self._setCoord:
            coord = Coordinate()
            self.rowToCoord(row, coord)
            if self._coordExtraLabels:
                partExtraLabels = self._extraLabels
                self.createExtraLabels(coord, row, self._coordExtraLabels)
                self._extraLabels = partExtraLabels
            particle.setCoordinate(coord)

        if self._postprocessImageRow:
            self._postprocessImageRow(particle, row)

    @staticmethod
    def rowToCoord(row, coord):
        """ Create a Coordinate from the row. """
        coord.setPosition(row.rlnCoordinateX, row.rlnCoordinateY)
        micName = getattr(row, 'rlnMicrographName', None)
        if micName is not None:
            coord.setMicName(os.path.basename(micName))

    @staticmethod
    def rowToCtf(row, ctf):
        """ Fill a CTFModel from the row. """
        ctf.setDefocusU(row.rlnDefocusU)
        ctf.setDefocusV(row.rlnDefocusV)
        ctf.setDefocusAngle(row.rlnDefocusAngle)
        ctf.setResolution(row.get('rlnCtfMaxResolution', 0))
        ctf.setFitQuality(row.get('rlnCtfFigureOfMerit', 0))

        if row.hasColumn('rlnPhaseShift'):
            ctf.setPhaseShift(row.rlnPhaseShift)
        ctf.standardize()

        if row.hasColumn('rlnCtfImage'):
            ctf.setPsdFile(row.rlnCtfImage)

    @staticmethod
    def rowToAcquisition(row, acq):
        """ Before 3.1 the acquisition was repeated on every particle
        row, there was no optics table to take it from. """
        for label, setter in [('rlnAmplitudeContrast',
                               acq.setAmplitudeContrast),
                              ('rlnSphericalAberration',
                               acq.setSphericalAberration),
                              ('rlnVoltage', acq.setVoltage)]:
            if row.hasColumn(label):
                setter(getattr(row, label))

    def setParticleTransform(self, particle, row):
        """ Set the transform values from the row. """

        if ((self._alignType == pwem.ALIGN_NONE) or
                not row.hasAnyColumn(self.ALIGNMENT_LABELS)):
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
                raise TypeError("Unexpected alignment type: %s"
                                % self._alignType)

        # Call again the modified function
        self.setParticleTransform(particle, row)

    def __setParticleTransformNone(self, particle, row):
        particle.setTransform(None)

    def __setParticleTransform2D(self, particle, row):
        angles = self._angles
        shifts = self._shifts

        def _get(label):
            return float(getattr(row, label, 0.))

        shifts[0] = _get('rlnOriginX')
        shifts[1] = _get('rlnOriginY')
        angles[2] = -_get('rlnAnglePsi')
        radAngles = -np.deg2rad(angles)
        M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
        M[:3, 3] = shifts[:3]
        particle.getTransform().setMatrix(M)

    def __setParticleTransformProj(self, particle, row):
        angles = self._angles
        shifts = self._shifts

        def _get(label):
            return float(getattr(row, label, 0.))

        shifts[0] = _get('rlnOriginX')
        shifts[1] = _get('rlnOriginY')
        shifts[2] = _get('rlnOriginZ')

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
