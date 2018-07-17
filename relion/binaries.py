
import os

import pyworkflow.utils as pwutils
from .constants import RELION_HOME, V2_0, V2_1


def getHome(*paths):
    """ Return the binary home path and possible some subfolders. """
    return os.path.join(os.environ[RELION_HOME], *paths)


def getEnviron():
    """ Setup the environment variables needed to launch Relion. """

    environ = pwutils.Environ(os.environ)
    binPath = getHome('bin')
    libPath = getHome('lib') + ":" + getHome('lib64')

    if not binPath in environ['PATH']:
        environ.update({'PATH': binPath,
                        'LD_LIBRARY_PATH': libPath,
                        'SCIPION_MPI_FLAGS': os.environ.get('RELION_MPI_FLAGS', ''),
                        }, position=pwutils.Environ.BEGIN)

    # Take Scipion CUDA library path
    cudaLib = environ.getFirst(('RELION_CUDA_LIB', 'CUDA_LIB'))
    environ.addLibrary(cudaLib)

    return environ


def getActiveVersion():
    """ Return the version of the Relion binaries that is currently active.
    In the current implementation it will be inferred from the RELION_HOME
    variable, so it should contain the version number in it. """
    home = getHome()
    for v in getSupportedVersions():
        if v in home:
            return v
    return ''


def isVersion2Active():
    return getActiveVersion().startswith("2.")


def getSupportedVersions():
    """ Return the list of supported binary versions. """
    return [V2_0, V2_1]


def validateInstallation():
    """ This function will be used to check if RELION binaries are
    properly installed. """
    environ = getEnviron()
    missingPaths = ["%s: %s" % (var, environ[var])
                    for var in [RELION_HOME]
                    if not os.path.exists(environ[var])]

    return (["Missing variables:"] + missingPaths) if missingPaths else []
