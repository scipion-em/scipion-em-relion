# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
import pyworkflow.utils as pwutils
import pwem

from .constants import RELION_HOME, RELION_CUDA_LIB, V3_0, V3_1


_logo = "relion_logo.png"
_references = ['Scheres2012a', 'Scheres2012b', 'Kimanius2016', 'Zivanov2018']


class Plugin(pwem.Plugin):
    _homeVar = RELION_HOME
    _supportedVersions = [V3_0, V3_1]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(RELION_HOME, 'relion-3.1')

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Relion. """
        environ = pwutils.Environ(os.environ)
        binPath = cls.getHome('bin')
        libPath = cls.getHome('lib') + ":" + cls.getHome('lib64')

        if binPath not in environ['PATH']:
            environ.update({'PATH': binPath,
                            'LD_LIBRARY_PATH': libPath,
                            'SCIPION_MPI_FLAGS': os.environ.get('RELION_MPI_FLAGS', ''),
                            }, position=pwutils.Environ.BEGIN)

        # Take Scipion CUDA library path
        cudaLib = environ.getFirst((RELION_CUDA_LIB, 'CUDA_LIB'))
        environ.addLibrary(cudaLib)

        if 'RELION_MPI_LIB' in os.environ:
            environ.addLibrary(os.environ['RELION_MPI_LIB'])

        if 'RELION_MPI_BIN' in os.environ:
            environ.set('PATH', os.environ['RELION_MPI_BIN'],
                        position=pwutils.Environ.BEGIN)

        return environ

    @classmethod
    def isVersion30Active(cls):
        return cls.getActiveVersion().startswith('3.0')

    @classmethod
    def isVersion31Active(cls):
        return cls.getActiveVersion().startswith('3.1')

    @classmethod
    def defineBinaries(cls, env):
        # Define FFTW3 path variables
        relion_vars = {'FFTW_LIB': env.getLibFolder(),
                       'FFTW_INCLUDE': env.getIncludeFolder()}

        relion_commands = [('cmake -DGUI=OFF -DCMAKE_INSTALL_PREFIX=./ .', []),
                           ('make -j %d' % env.getProcessors(),
                            ['bin/relion_refine'])]

        env.addPackage('relion', version='3.0',
                       url='https://github.com/3dem/relion/archive/3.0.tar.gz',
                       commands=relion_commands,
                       updateCuda=True,
                       vars=relion_vars)

        env.addPackage('relion', version='3.1',
                       url='https://github.com/3dem/relion/archive/3.1.tar.gz',
                       commands=relion_commands,
                       updateCuda=True,
                       vars=relion_vars,
                       default=True)


# Convenience constants to check what version of the binary is active
IS_30 = Plugin.isVersion30Active()
# Greater than 3.0?
IS_GT30 = not IS_30  # Now 2.x is deprecated
IS_31 = Plugin.isVersion31Active()
