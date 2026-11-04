#!/usr/bin/env python
# **************************************************************************
# *
# * Authors:     Mikel Iceta (miceta@cnb.csic.es)                     [1]
# *
# * [1] Unidad de  Biocomputacion, Centro Nacional de Biotecnologia, CSIC (CNB-CSIC)
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

""" Picking the preferential views to remove.

These tests drive the angle matching with particles built at known Euler
angles, and cover the two things that are easy to get wrong: which
particles fall within the limits, and that the selection survives being
read back by another protocol instance, as it is after a continue. They
need neither RELION binaries nor a GPU.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestOutput
from pwem.objects import Particle, SetOfParticles, Transform, Acquisition
import pwem.convert.transformations as tfs

from ..protocols import ProtRelionRemovePrefViews
from ..wizards import RelionWizAngDistDisplay


def transformFor(rot, tilt, psi):
    """ The transform a particle needs for the protocol to read back these
    Euler angles, that is, the inverse of what processAnglesStep does. """
    inverse = tfs.euler_matrix(*(-np.deg2rad([rot, tilt, psi])), axes='szyz')

    return Transform(np.linalg.inv(inverse))


class TestRelionRemovePrefViewsAngles(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _createSetOfParticles(self, angles):
        """ One particle per (rot, tilt, psi), with the object id of its
        position in the list, starting at 1. """
        fn = self.getOutputPath('%s.sqlite' % self._testMethodName)
        pwutils.cleanPath(fn)
        parts = SetOfParticles(filename=fn)
        parts.setSamplingRate(3.5)
        parts.setAcquisition(Acquisition(voltage=300, sphericalAberration=2.7,
                                         amplitudeContrast=0.1,
                                         magnification=50000))

        for i, (rot, tilt, psi) in enumerate(angles):
            part = Particle()
            part.setObjId(i + 1)
            part.setLocation(i + 1, 'particles.mrcs')
            part.setSamplingRate(3.5)
            part.setTransform(transformFor(rot, tilt, psi))
            parts.append(part)

        parts.setAlignmentProj()
        parts.write()

        return parts

    def _protocol(self, angles, **kwargs):
        """ A protocol with its own working dir, ready to run the step. """
        workingDir = self.getOutputPath(self._testMethodName)
        pwutils.makePath(os.path.join(workingDir, 'extra'))

        prot = ProtRelionRemovePrefViews()
        prot.setWorkingDir(workingDir)
        prot.inputParticles.set(self._createSetOfParticles(angles))
        for name, value in kwargs.items():
            getattr(prot, name).set(value)

        return prot

    def test_picksOnlyTheParticlesWithinTheLimits(self):
        """ Particles 2 and 3 sit inside the rot/tilt window. """
        prot = self._protocol([(-150., 30., 0.),    # rot below the limits
                               (10., 30., 0.),
                               (20., 40., 0.),
                               (10., 150., 0.)],    # tilt above the limits
                              numToRemove=4,
                              rotMin=0, rotMax=90, tiltMin=0, tiltMax=90)
        prot.processAnglesStep()

        self.assertEqual({2, 3}, prot._readRemovedIds())

    def test_limitsAreExclusive(self):
        """ withinLimits() is a strict comparison on both ends. """
        prot = self._protocol([(30., 45., 0.),
                               (30., 90., 0.)],   # exactly on the tilt limit
                              numToRemove=2,
                              rotMin=0, rotMax=90, tiltMin=0, tiltMax=90)
        prot.processAnglesStep()

        self.assertEqual({1}, prot._readRemovedIds())

    def test_psiIsOnlyCheckedWhenAsked(self):
        angles = [(30., 45., 0.), (30., 45., 170.)]
        limits = dict(rotMin=0, rotMax=90, tiltMin=0, tiltMax=90,
                      psiMin=-90, psiMax=90)

        prot = self._protocol(angles, numToRemove=2, removePsi=False, **limits)
        prot.processAnglesStep()
        self.assertEqual({1, 2}, prot._readRemovedIds())

        prot = self._protocol(angles, numToRemove=2, removePsi=True, **limits)
        prot.processAnglesStep()
        self.assertEqual({1}, prot._readRemovedIds())

    def test_removesOnlyTheRequestedNumber(self):
        prot = self._protocol([(30., 45., 0.), (31., 46., 0.),
                               (32., 47., 0.), (33., 48., 0.)],
                              numToRemove=2,
                              rotMin=0, rotMax=90, tiltMin=0, tiltMax=90)
        prot.processAnglesStep()

        removed = prot._readRemovedIds()
        self.assertEqual(2, len(removed))
        self.assertTrue(removed.issubset({1, 2, 3, 4}))

    def test_asksForMoreThanAvailable(self):
        """ Nothing to sample from: take every particle that matched. """
        prot = self._protocol([(30., 45., 0.), (-150., 30., 0.)],
                              numToRemove=10,
                              rotMin=0, rotMax=90, tiltMin=0, tiltMax=90)
        prot.processAnglesStep()

        self.assertEqual({1}, prot._readRemovedIds())

    def test_nothingMatchesTheLimits(self):
        prot = self._protocol([(30., 45., 0.)], numToRemove=1,
                              rotMin=100, rotMax=180, tiltMin=0, tiltMax=90)
        prot.processAnglesStep()

        self.assertEqual(set(), prot._readRemovedIds())

    def _runOutputStep(self, prot):
        """ The real output step, with the two calls that need a project
        stubbed out: what is under test is which particles it copies. """
        created = {}
        prot._defineOutputs = lambda **kwargs: created.update(kwargs)
        prot._defineTransformRelation = lambda *args: None
        prot.createOutputStep()

        return created['outputParticles']

    def test_theOutputLeavesOutTheSelectedParticles(self):
        prot = self._protocol([(30., 45., 0.),      # within
                               (-150., 30., 0.),
                               (20., 40., 0.),      # within
                               (10., 150., 0.)],
                              numToRemove=4,
                              rotMin=0, rotMax=90, tiltMin=0, tiltMax=90)
        prot.processAnglesStep()

        output = self._runOutputStep(prot)

        self.assertEqual(2, output.getSize())
        self.assertEqual([2, 4], [p.getObjId() for p in output.iterItems()])

    def test_theSelectionSurvivesAContinue(self):
        """ The output step runs on a protocol rebuilt from the database,
        so the selection cannot live in an attribute of the instance. """
        angles = [(30., 45., 0.), (-150., 30., 0.), (20., 40., 0.)]
        prot = self._protocol(angles, numToRemove=3,
                              rotMin=0, rotMax=90, tiltMin=0, tiltMax=90)
        prot.processAnglesStep()

        continued = ProtRelionRemovePrefViews()
        continued.setWorkingDir(prot.getWorkingDir())

        self.assertEqual({1, 3}, continued._readRemovedIds())


class TestRelionRemovePrefViewsSummary(BaseTest):
    """ The summary reads two different kinds of attribute: inputParticles
    is a Pointer and needs get(), outputParticles is the output Set and
    does not. Getting either one wrong raises instead of reporting.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def _set(self, name, size):
        fn = self.getOutputPath('%s_%s.sqlite' % (self._testMethodName, name))
        pwutils.cleanPath(fn)
        parts = SetOfParticles(filename=fn)
        parts.setSamplingRate(3.5)

        for i in range(size):
            part = Particle()
            part.setObjId(i + 1)
            part.setLocation(i + 1, 'particles.mrcs')
            part.setSamplingRate(3.5)
            parts.append(part)

        parts.write()

        return parts

    def test_summaryCountsBothSets(self):
        prot = ProtRelionRemovePrefViews()
        prot.inputParticles.set(self._set('input', 10))
        prot.outputParticles = self._set('output', 7)

        self.assertEqual(["Input particles: 10", "Output particles: 7"],
                         prot._summary())

    def test_quietBeforeThereIsAnOutput(self):
        prot = ProtRelionRemovePrefViews()
        prot.inputParticles.set(self._set('input', 10))

        self.assertEqual([], prot._summary())


class TestRelionAngDistWizard(BaseTest):
    """ The wizard asks pwem for the plot, and that API moved from a
    histogram flag to a plot type. Only drawing catches a stale call, so
    the backend is switched to a non-interactive one and put back.
    """
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
        cls._backend = plt.get_backend()
        plt.switch_backend('Agg')

    @classmethod
    def tearDownClass(cls):
        plt.switch_backend(cls._backend)

    def test_plotsTheAngularDistribution(self):
        fn = self.getOutputPath('%s.sqlite' % self._testMethodName)
        pwutils.cleanPath(fn)
        parts = SetOfParticles(filename=fn)
        parts.setSamplingRate(3.5)

        rng = np.random.default_rng(0)
        for i in range(50):
            quat = rng.normal(size=4)
            part = Particle()
            part.setObjId(i + 1)
            part.setLocation(i + 1, 'particles.mrcs')
            part.setSamplingRate(3.5)
            part.setTransform(
                Transform(tfs.quaternion_matrix(quat / np.linalg.norm(quat))))
            parts.append(part)

        parts.setAlignmentProj()
        parts.write()

        plotter = RelionWizAngDistDisplay.getPlotter(parts)

        self.assertIsNotNone(plotter)
        plotter.close()
