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

import os
import numpy as np
from random import sample
from enum import Enum

import pyworkflow.protocol.params as params
from pyworkflow.constants import BETA
from pwem.protocols import EMProtocol
from pwem.objects import SetOfParticles
import pwem.convert.transformations as tfs

from emtools.utils import Process


class outputs(Enum):
    outputParticles = SetOfParticles


class ProtParticlesToScratch(EMProtocol):
    """ Protocol to copy particles to Scratch
    """
    _label = 'particles to scratch'
    _devStatus = BETA
    _possibleOutputs = outputs

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles.')

    # -------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    def createOutputStep(self):
        projName = self.getProject().shortName
        logger = Process.Logger()
        imgSet = self.inputParticles.get()
        fnMap = {}
        dirMap = {}

        def _updatePath(item, row):
            fn = item.getFileName()
            dn, base = os.path.split(fn)
            dstDir = os.path.join(os.environ['SCIPION_SCRATCH'], projName, dn)

            if fn not in fnMap:
                fnMap[fn] = True
                if dn not in dirMap:
                    logger.mkdir(dstDir)
                    dirMap[dn] = dstDir
                logger.cp(fn, dstDir)

            item.setFileName(os.path.join(dstDir, base))

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        imgSet.setAlignment(imgSet.getAlignment())
        outImgSet.copyItems(imgSet, updateItemCallback=_updatePath)
        self._defineOutputs(**{outputs.outputParticles.name: outImgSet})
        self._defineTransformRelation(self.inputParticles, outImgSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if hasattr(self, "outputParticles"):
            summary.append("Input particles: %d" % self.inputParticles.get().getSize())
            summary.append("Output particles: %d" % self.outputParticles.get().getSize())

        return summary

    def _validate(self):
        errors = []
        return errors
