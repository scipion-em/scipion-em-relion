# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
import pyworkflow.em
import pyworkflow.utils as pwutils

from .constants import RELION_HOME, V2_0, V2_1


_logo = "relion_logo.png"


# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented
class Plugin:
    __metaclass__ = pyworkflow.em.PluginMeta

    @classmethod
    def __getHome(cls, *paths):
        """ Return the binary home path and possible some subfolders. """
        return os.path.join(os.environ[RELION_HOME], *paths)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch Relion. """

        environ = pwutils.Environ(os.environ)
        binPath = cls.__getHome('bin')
        libPath = cls.__getHome('lib') + ":" + cls.__getHome('lib64')

        if not binPath in environ['PATH']:
            environ.update({'PATH': binPath,
                            'LD_LIBRARY_PATH': libPath,
                            'SCIPION_MPI_FLAGS': os.environ.get('RELION_MPI_FLAGS', ''),
                            }, position=pwutils.Environ.BEGIN)

        # Take Scipion CUDA library path
        cudaLib = environ.getFirst(('RELION_CUDA_LIB', 'CUDA_LIB'))
        environ.addLibrary(cudaLib)

        return environ

    @classmethod
    def getActiveVersion(cls):
        """ Return the version of the Relion binaries that is currently active.
        In the current implementation it will be inferred from the RELION_HOME
        variable, so it should contain the version number in it. """
        home = cls.__getHome()
        for v in cls.getSupportedVersions():
            if v in home:
                return v
        return ''

    @classmethod
    def isVersion2Active(cls):
        return cls.getActiveVersion().startswith("2.")

    @classmethod
    def getSupportedVersions(cls):
        """ Return the list of supported binary versions. """
        return [V2_0, V2_1]

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if RELION binaries are
        properly installed. """
        environ = cls.getEnviron()
        missingPaths = ["%s: %s" % (var, environ[var])
                        for var in [RELION_HOME]
                        if not os.path.exists(environ[var])]

        return (["Missing variables:"] + missingPaths) if missingPaths else []





