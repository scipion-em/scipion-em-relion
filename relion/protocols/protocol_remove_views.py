# ******************************************************************************
# *
# * Authors:     Grigory Sharov     (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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
# ******************************************************************************

import numpy as np
from random import SystemRandom
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.protocols import ProtParticles
from pwem.objects import SetOfParticles
import pwem.convert.transformations as tfs


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtRelionRemovePrefViews(ProtParticles):
    """ Protocol to remove preferential views from a particle set.

    Inspired by https://github.com/leschzinerlab/Relion

    """
    _label = 'remove preferential views'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles.')
        form.addParam('numToRemove', params.IntParam,
                      default=0,
                      label="Number of particles to remove",
                      help="Number of particles to remove "
                           "WITHIN the limits below.")

        group = form.addGroup("Remove views WITHIN the limits below")
        line = group.addLine('Rot')
        line.addParam('rotMin', params.IntParam, default=-180,
                      label='min')
        line.addParam('rotMax', params.IntParam, default=180,
                      label='max')

        line = group.addLine('Tilt')
        line.addParam('tiltMin', params.IntParam, default=0,
                      label='min')
        line.addParam('tiltMax', params.IntParam, default=180,
                      label='max')

        group.addParam('removePsi', params.BooleanParam, default=False,
                       label="Remove views with specific in-plane rotation?",
                       help="Particle orientation on Euler sphere is "
                            "defined by rot and tilt angles. Psi is for "
                            "in-plane rotation only. Select *Yes* "
                            "if you want to provide psi limits.")
        line = group.addLine('Psi', condition='removePsi')
        line.addParam('psiMin', params.IntParam, default=-180,
                      condition='removePsi', label='min')
        line.addParam('psiMax', params.IntParam, default=180,
                      condition='removePsi', label='max')

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.processAnglesStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    def processAnglesStep(self):
        """ Write down the ids of the particles to leave out of the output. """
        ptclToRemove = self.numToRemove.get()
        rmin, rmax = self.rotMin.get(), self.rotMax.get()
        tmin, tmax = self.tiltMin.get(), self.tiltMax.get()
        pmin, pmax = self.psiMin.get(), self.psiMax.get()

        # The angles are computed right here, and the ids are written to a
        # file, because nothing kept in self survives a continue: with the
        # angles in a step of their own, continuing a run skipped it and
        # left this one matching nothing, quietly copying the whole input
        # to the output.
        matched = []
        for part in self.inputParticles.get().iterItems():
            matrix = np.linalg.inv(part.getTransform().getMatrix())
            rot, tilt, psi = -np.rad2deg(tfs.euler_from_matrix(matrix,
                                                               axes='szyz'))
            if self.removePsi:
                within = (self.withinLimits(rot, rmin, rmax) and
                          self.withinLimits(tilt, tmin, tmax) and
                          self.withinLimits(psi, pmin, pmax))
            else:
                # check only rot & tilt
                within = (self.withinLimits(rot, rmin, rmax) and
                          self.withinLimits(tilt, tmin, tmax))

            if within:
                matched.append(part.getObjId())

        self.info("Found %d particles within the specified limits."
                  % len(matched))

        if ptclToRemove > len(matched):
            self.info("Number to remove (%d) is more than maximum available (%d). "
                      "Removing %d particles..." % (ptclToRemove,
                                                    len(matched),
                                                    len(matched)))
        else:
            self.info("Randomly removing %d particles within "
                      "specified limits..." % ptclToRemove)
            matched = SystemRandom().sample(matched, ptclToRemove)

        self._writeRemovedIds(matched)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        self.removedIds = self._readRemovedIds()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        imgSet.setAlignmentProj()
        outImgSet.copyItems(imgSet, updateItemCallback=self._removeViews)
        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(self.inputParticles, outImgSet)

    def _removeViews(self, item, row):
        if item.getObjId() in self.removedIds:
            setattr(item, "_appendItem", False)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, "outputParticles"):
            summary.append("Input particles: %d" % self.inputParticles.get().getSize())
            # inputParticles is a Pointer and needs the get(), outputParticles
            # is the output Set itself and chokes on it
            summary.append("Output particles: %d" % self.outputParticles.getSize())

        return summary

    def _validate(self):
        errors = []

        if self.numToRemove > self.inputParticles.get().getSize():
            errors.append("You cannot remove more particles "
                          "than provided in the input.")

        if (self.rotMin < -180 or self.psiMin < -180 or
                self.rotMax > 180 or self.tiltMax > 180 or
                self.psiMax > 180 or self.tiltMin < 0):
            errors.append("Angles must be within the following limits:\n\n"
                          "-180 < rot < 180\n0 < tilt < 180\n-180 < psi < 180\n")

        return errors

    # -------------------------- UTILS functions ------------------------------
    def withinLimits(self, value, min, max):
        return min < float(value) < max

    def _getRemovedIdsFile(self):
        return self._getExtraPath("removed_ids.txt")

    def _writeRemovedIds(self, ids):
        with open(self._getRemovedIdsFile(), "w") as f:
            f.write("".join("%d\n" % objId for objId in ids))

    def _readRemovedIds(self):
        """ Ids of the particles to leave out, as a set, since the output
        step asks for every single input particle. """
        with open(self._getRemovedIdsFile()) as f:
            return {int(line) for line in f if line.strip()}
