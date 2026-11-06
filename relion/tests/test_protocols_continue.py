# **************************************************************************
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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""Continue processes again what failed, on real runs.

An item that cannot be read fails alone and the run still ends well.
Once it is readable again, Continue processes it and leaves the others
alone.
"""

import os
import shutil

from pyworkflow.protocol.constants import MODE_RESUME
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr

from ..protocols import (ProtRelionAutopickLoG, ProtRelionExtractParticles,
                         ProtRelionMotioncor)
from .test_protocols_base import TestRelionBase, CPUS

MIC_NAMES = ('006.mrc', '008.mrc', '016.mrc')
MOVIE_NAMES = ('20170629_00021_frameImage.tiff',
               '20170629_00022_frameImage.tiff',
               '20170629_00023_frameImage.tiff')


def _break(fileName):
    os.remove(fileName)

    with open(fileName, 'wb') as handle:
        handle.write(b'not an image')


def _countsBy(outputSet, attribute):
    return {int(agg[attribute]): agg['count']
            for agg in outputSet.aggregate(["count"], attribute, [attribute])}


class TestRelionContinueAfterFailedItems(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.micsSource = DataSet.getDataSet('relion_tutorial').getFile('micrographs')
        cls.micsDir = cls._copyInto('mics', cls.micsSource, MIC_NAMES)
        cls.protImportMics = cls.runImportMics(
            os.path.join(cls.micsDir, '*.mrc'), samplingRate=7.08)

    @classmethod
    def _copyInto(cls, folder, sourceDir, names):
        targetDir = os.path.abspath(cls.proj.getTmpPath(folder))
        os.makedirs(targetDir)

        for name in names:
            shutil.copy(os.path.join(sourceDir, name), targetDir)

        return targetDir

    def _restore(self, folder, sourceDir, name):
        target = os.path.join(folder, name)
        os.remove(target)
        shutil.copy(os.path.join(sourceDir, name), target)

    def _continue(self, prot):
        prot.runMode.set(MODE_RESUME)
        self.launchProtocol(prot)

    def _pick(self, label):
        prot = self.newProtocol(ProtRelionAutopickLoG, objLabel=label,
                                inputMicrographs=self.protImportMics.outputMicrographs,
                                boxSize=64, minDiameter=260, maxDiameter=360,
                                streamingBatchSize=1)
        self.launchProtocol(prot)
        return prot

    def testAFailedMicrographIsPickedOnContinue(self):
        broken = MIC_NAMES[1]
        _break(os.path.join(self.micsDir, broken))

        print(magentaStr("\n==> Testing relion - LoG with a broken micrograph:"))
        prot = self._pick('LoG broken mic')
        picked = _countsBy(prot.outputCoordinates, '_micId')
        self.assertEqual(len(MIC_NAMES) - 1, len(picked),
                         "Only the broken micrograph is left out.")

        self._restore(self.micsDir, self.micsSource, broken)
        print(magentaStr("\n==> Continue once the micrograph is readable:"))
        self._continue(prot)

        pickedAgain = _countsBy(prot.outputCoordinates, '_micId')
        self.assertEqual(len(MIC_NAMES), len(pickedAgain),
                         "Continue picks the micrograph that failed.")

        for micId, count in picked.items():
            self.assertEqual(count, pickedAgain[micId],
                             "A micrograph picked before is not picked twice.")

    def testAFailedExtractionIsDoneOnContinue(self):
        protPick = self._pick('LoG for extraction')
        broken = MIC_NAMES[2]
        _break(os.path.join(self.micsDir, broken))

        print(magentaStr("\n==> Testing relion - extraction with a broken micrograph:"))
        prot = self.newProtocol(ProtRelionExtractParticles, boxSize=64,
                                doInvert=False, streamingBatchSize=1)
        prot.inputCoordinates.set(protPick.outputCoordinates)
        self.launchProtocol(prot)
        extracted = _countsBy(prot.outputParticles, '_micId')
        self.assertEqual(len(MIC_NAMES) - 1, len(extracted),
                         "Only the broken micrograph is left out.")

        self._restore(self.micsDir, self.micsSource, broken)
        print(magentaStr("\n==> Continue once the micrograph is readable:"))
        self._continue(prot)

        extractedAgain = _countsBy(prot.outputParticles, '_micId')
        self.assertEqual(len(MIC_NAMES), len(extractedAgain),
                         "Continue extracts the micrograph that failed.")

        for micId, count in extracted.items():
            self.assertEqual(count, extractedAgain[micId],
                             "A micrograph extracted before is not "
                             "extracted twice.")


class TestRelionMotioncorContinueAfterFailedMovies(TestRelionBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        ds = DataSet.getDataSet('relion30_tutorial')
        cls.moviesSource = ds.getFile('Movies/')
        cls.moviesDir = os.path.abspath(cls.proj.getTmpPath('movies'))
        os.makedirs(cls.moviesDir)

        # Linked, not copied: they are large, and a broken one is replaced
        # by a file of its own.
        for name in MOVIE_NAMES:
            os.symlink(os.path.join(cls.moviesSource, name),
                       os.path.join(cls.moviesDir, name))

        cls.protImport = cls.runImportMovies(filesPath=cls.moviesDir,
                                             filesPattern='*.tiff',
                                             samplingRate=0.885, voltage=200,
                                             sphericalAberration=1.4,
                                             dose=1.277,
                                             gain=ds.getFile("Movies/gain.mrc"))

    def testAFailedMovieIsAlignedOnContinue(self):
        broken = os.path.join(self.moviesDir, MOVIE_NAMES[1])
        _break(broken)

        print(magentaStr("\n==> Testing relion - motioncor with a broken movie:"))
        prot = self.newProtocol(ProtRelionMotioncor, objLabel='mc broken movie',
                                patchX=1, patchY=1, doDW=True,
                                numberOfThreads=CPUS)
        prot.inputMovies.set(self.protImport.outputMovies)
        self.launchProtocol(prot)

        self.assertEqual(len(MOVIE_NAMES) - 1, prot.outputMovies.getSize(),
                         "Only the broken movie is left out.")
        self.assertEqual(prot.outputMovies.getSize(),
                         prot.outputMicrographsDoseWeighted.getSize(),
                         "Every output has the same movies.")

        os.remove(broken)
        os.symlink(os.path.join(self.moviesSource, MOVIE_NAMES[1]), broken)
        print(magentaStr("\n==> Continue once the movie is readable:"))
        prot.runMode.set(MODE_RESUME)
        self.launchProtocol(prot)

        self.assertEqual(len(MOVIE_NAMES), prot.outputMovies.getSize(),
                         "Continue aligns the movie that failed.")
        self.assertEqual(len(MOVIE_NAMES),
                         prot.outputMicrographsDoseWeighted.getSize())
