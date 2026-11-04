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

""" Which RELION version is installed.

Reading it off the name of the RELION_HOME folder is a guess, and when the
guess fails everything behind IS_GT50() disappears silently. These tests
pin down that the binary is asked first and that the folder name is still
there as a fallback. They need no RELION binaries: a stub script stands in
for relion_refine.
"""

import os

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestOutput

from relion import Plugin
from ..constants import RELION_HOME


class TestRelionActiveVersion(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)

    def setUp(self):
        self._home = Plugin._vars.get(RELION_HOME, None)
        self._cached = Plugin._activeVersion
        self._fromBinary = Plugin._getVersionFromBinary

    def tearDown(self):
        Plugin._vars[RELION_HOME] = self._home
        Plugin._activeVersion = self._cached
        Plugin._getVersionFromBinary = self._fromBinary

    def _setHome(self, name, binaryPrints=None):
        """ Point RELION_HOME at a folder with the given name, holding a
        bin/relion_refine that prints what it is told to (or no binary at
        all when nothing is given). """
        home = self.getOutputPath(self._testMethodName, name)
        pwutils.cleanPath(home)
        pwutils.makePath(os.path.join(home, 'bin'))

        if binaryPrints is not None:
            program = os.path.join(home, 'bin', 'relion_refine')
            with open(program, 'w') as f:
                f.write('#!/bin/bash\necho "%s"\n' % binaryPrints)
            os.chmod(program, 0o755)

        Plugin._vars[RELION_HOME] = home
        Plugin._activeVersion = None

        return home

    def test_readsTheVersionFromTheBinary(self):
        self._setHome('relion-5.1',
                      'RELION version: 5.1.0-commit-cb2683')

        self.assertEqual('5.1', Plugin._getVersionFromBinary())

    def test_theBinaryWinsOverTheFolderName(self):
        """ The folder name is only a hint, and it can be plain wrong: a
        5.1 build is often installed over a relion-5.0 directory. """
        self._setHome('relion-5.0',
                      'RELION version: 5.1.0-commit-cb2683')

        self.assertEqual('5.1', Plugin.getActiveVersion())

    def test_fallsBackToTheFolderNameWithoutABinary(self):
        self._setHome('relion-5.0')

        self.assertIsNone(Plugin._getVersionFromBinary())
        self.assertEqual('5.0', Plugin.getActiveVersion())

    def test_unnamedFolderIsNoLongerVersionless(self):
        """ The reason for all this: /opt/relion used to report nothing,
        and that turned off Blush, DynaMight and ModelAngelo. """
        self._setHome('relion', 'RELION version: 5.1.0-commit-cb2683')

        self.assertEqual('5.1', Plugin.getActiveVersion())
        self.assertTrue(Plugin.IS_GT50())

    def test_unnamedFolderWithoutABinaryStillReportsNothing(self):
        """ Nothing to go on: keep the old answer rather than invent one. """
        self._setHome('relion')

        self.assertEqual('', Plugin.getActiveVersion())
        self.assertFalse(Plugin.IS_GT50())

    def test_ignoresAVersionWeDoNotSupport(self):
        self._setHome('relion-5.1', 'RELION version: 9.9.0-commit-deadbee')

        self.assertIsNone(Plugin._getVersionFromBinary())
        self.assertEqual('5.1', Plugin.getActiveVersion())

    def test_survivesABinaryThatSaysSomethingElse(self):
        self._setHome('relion-5.1', 'Segmentation fault')

        self.assertIsNone(Plugin._getVersionFromBinary())
        self.assertEqual('5.1', Plugin.getActiveVersion())

    def test_theAnswerIsAskedForOnlyOnce(self):
        """ getActiveVersion() is called from isDisabled() at import time,
        so it must not spawn a process every time. """
        self._setHome('relion-5.0')
        calls = []

        def counted():
            calls.append(1)
            return '5.1'

        Plugin._getVersionFromBinary = counted

        self.assertEqual('5.1', Plugin.getActiveVersion())
        self.assertEqual('5.1', Plugin.getActiveVersion())
        self.assertEqual(1, len(calls))
