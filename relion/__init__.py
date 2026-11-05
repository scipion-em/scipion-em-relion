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
import re
import subprocess

import pyworkflow.utils as pwutils
from pyworkflow import Config, SPA, TOMO
import pwem

from .constants import *


__version__ = '5.1.1'
_logo = "relion_logo.jpg"
_references = ['Scheres2012a', 'Scheres2012b', 'Kimanius2016',
               'Zivanov2018', 'Kimanius2021', 'Kimanius2024']


class Plugin(pwem.Plugin):
    _homeVar = RELION_HOME
    _pathVars = [RELION_HOME]
    _supportedVersions = [V4_0, V5_0, V5_1]
    # getActiveVersion() runs a binary, so the answer is kept
    _activeVersion = None
    _url = "https://github.com/scipion-em/scipion-em-relion"
    _processingField = [SPA, TOMO]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(RELION_HOME, f'relion-{V5_1}')
        cls._defineVar(RELION_CUDA_LIB, pwem.Config.CUDA_LIB)
        cls._defineVar(RELION_CUDA_BIN, pwem.Config.CUDA_BIN)
        cls._defineVar(RELION_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD % V5_1)
        cls._defineVar(RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE,
                       os.getenv(RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE, None))
        cls._defineEmVar(TORCH_HOME_VAR, 'modelangelomodels-1.0')
        cls._defineVar(RELION_PREPEND, '',
                       description="Prepended to every RELION command of a "
                                   "job that uses the GPU, for loading "
                                   "modules and the like")
        cls._defineVar(RELION_PREPEND_CPU, '',
                       description="Prepended to every RELION command of a "
                                   "job that does not use the GPU")

    @classmethod
    def getEnviron(cls, useGpu=False):
        """ Setup the environment variables needed to launch Relion.

        The prepend is set here and not in ProtRelionBase so that it also
        reaches the protocols that do not inherit from it: Protocol
        ._getEnviron() falls back to this method for all of them.
        """
        environ = pwutils.Environ(os.environ)
        binPath = os.pathsep.join([cls.getHome('bin'),
                                   pwem.Config.MPI_BINDIR])
        libPath = os.pathsep.join([cls.getHome('lib'),
                                   pwem.Config.MPI_LIBDIR])

        if binPath not in environ['PATH']:
            environ.update({'PATH': binPath,
                            'LD_LIBRARY_PATH': libPath
                            }, position=pwutils.Environ.BEGIN)

        # Get Relion CUDA library path if defined
        cudaLib = cls.getVar(RELION_CUDA_LIB, pwem.Config.CUDA_LIB)
        environ.addLibrary(cudaLib)

        if 'RELION_MPI_LIB' in os.environ:
            environ.addLibrary(os.environ['RELION_MPI_LIB'])

        if 'RELION_MPI_BIN' in os.environ:
            environ.set('PATH', os.environ['RELION_MPI_BIN'],
                        position=pwutils.Environ.BEGIN)

        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']

        # Set SIDESPLITTER env var if possible
        if cls.getVar(SIDESPLITTER_HOME, None) is not None:
            environ.update({
                SIDESPLITTER:
                    os.path.join(cls.getVar(SIDESPLITTER_HOME), 'sidesplitter')
            })
            if cls.getVar(RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE) is None:
                environ.update({
                    RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE:
                        os.path.join(cls.getVar(SIDESPLITTER_HOME),
                                     'sidesplitter_wrapper.sh')
                })

        if TORCH_HOME_VAR not in environ:
            environ.set(TORCH_HOME_VAR, Plugin.getVar(TORCH_HOME_VAR))

        environ.setPrepend(cls.getVar(RELION_PREPEND if useGpu
                                      else RELION_PREPEND_CPU) or '')

        return environ

    @classmethod
    def getActiveVersion(cls, *args, **kwargs):
        """ Version of the binaries that are currently active.

        The inherited implementation reads it off the basename of
        RELION_HOME, so a perfectly good installation in /opt/relion or in
        relion-5.1-cuda12 reports no version at all. Everything behind
        IS_GT50() then disappears without a word: Blush, DynaMight and
        ModelAngelo. Ask the binary first and keep the folder name only as
        a fallback for when it cannot be run.
        """
        if cls._activeVersion is None:
            cls._activeVersion = (cls._getVersionFromBinary() or
                                  super().getActiveVersion(*args, **kwargs))

        return cls._activeVersion

    @classmethod
    def _getVersionFromBinary(cls):
        """ Ask relion_refine, which prints a line like
        'RELION version: 5.1.0-commit-cb2683'. Returns one of the supported
        versions, or None if the binary is missing or says something else.
        """
        program = cls.getHome('bin', 'relion_refine')

        if not os.path.exists(program):
            return None

        try:
            out = subprocess.run([program, '--version'], timeout=60,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT).stdout.decode()
        except Exception as e:
            pwutils.yellowStr("Could not ask %s for its version: %s"
                              % (program, e))
            return None

        match = re.search(r"RELION version:\s*(\d+\.\d+)", out)
        version = match.group(1) if match else None

        return version if version in cls._supportedVersions else None

    @classmethod
    def IS_GT50(cls):
        return cls.getActiveVersion().startswith('5')

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
                activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = ['git', 'gcc', 'cmake', 'make']
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def getRelionEnvActivation(cls):
        """ Remove the scipion home and activate the conda environment. """
        activation = cls.getVar(RELION_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep

        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getActivationCmd(cls):
        """ Return the activation command. """
        return '%s %s' % (cls.getCondaActivationCmd(),
                          cls.getRelionEnvActivation())

    @classmethod
    def defineBinaries(cls, env):
        for ver in cls._supportedVersions:
            torchHome = cls.getVar(TORCH_HOME_VAR)
            os.makedirs(torchHome, exist_ok=True)

            if ver == V4_0:
                envCmd = f"conda create -y -n relion-4.0 -c pytorch python=3.9 pytorch=1.10.0 numpy=1.20 &&"
            else:
                envCmd = "conda env create -f environment.yml -y &&"

            cmd = [
                f'cd .. && rmdir relion-{ver} &&',
                f'git clone https://github.com/3dem/relion.git relion-{ver} &&',
                f'cd relion-{ver} && git checkout ver{ver} &&',
                cls.getCondaActivationCmd(),
                envCmd,
                f'cmake -DCMAKE_INSTALL_PREFIX=./ -DTORCH_HOME_PATH={torchHome} .'
            ]

            installCmds = [
                (" ".join(cmd), ['Makefile']),
                (f'make -j {env.getProcessors()}', ['bin/relion_refine'])
            ]

            env.addPackage('relion', version=ver,
                           tar='void.tgz',
                           commands=installCmds,
                           neededProgs=cls.getDependencies(),
                           updateCuda=True,
                           default=ver == V5_1)
