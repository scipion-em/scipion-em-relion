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
"""Choices for building RELION that depend on the machine it is built on.

nvcc refuses host compilers newer than the ones its CUDA release supports
("unsupported GNU version!"), and a recent distribution can ship a default
gcc that the installed CUDA does not accept yet (Ubuntu 24.04: gcc 13 with
CUDA 11.8). RELION's GUI needs the Xft development files, without which
its bundled FLTK aborts the whole build.
"""

import os
import re
import shutil
import subprocess

# Newest gcc major each CUDA release accepts as host compiler, from the
# version check in its include/crt/host_config.h.
_MAX_GCC_BY_CUDA = (
    ((10, 1), 8),
    ((11, 0), 9),
    ((11, 1), 10),
    ((11, 4), 11),
    ((12, 0), 12),
    ((12, 4), 13),
    ((12, 8), 14),
)

# The oldest gcc worth trying: RELION needs C++17.
_OLDEST_GCC = 8


def maxGccForCuda(cudaVersion):
    """Newest gcc major ``cudaVersion`` (major, minor) accepts, or None if unknown."""
    if cudaVersion is None or cudaVersion < _MAX_GCC_BY_CUDA[0][0]:
        return None

    maxGcc = None

    for release, gcc in _MAX_GCC_BY_CUDA:
        if cudaVersion >= release:
            maxGcc = gcc

    return maxGcc


def _run(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True,
                              timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return ''


def cudaVersion(nvcc):
    """(major, minor) of the CUDA behind ``nvcc``, or None."""
    match = re.search(r'release (\d+)\.(\d+)', _run([nvcc, '--version']))

    return (int(match.group(1)), int(match.group(2))) if match else None


def gccMajor(compiler):
    """Major version of a gcc, or None when it is not there or not gcc."""
    match = re.match(r'(\d+)', _run([compiler, '-dumpversion']).strip())

    return int(match.group(1)) if match else None


def cudaHostCompiler(nvcc, defaultCompiler=None, which=shutil.which):
    """A gcc nvcc accepts, when the default one is too new for it.

    None when there is nothing to change: no CUDA, a default compiler
    CUDA accepts, an unknown CUDA release, or no suitable gcc installed
    (then the build reports CUDA's own error, as before).
    """
    maxGcc = maxGccForCuda(cudaVersion(nvcc)) if nvcc else None

    if maxGcc is None:
        return None

    defaultCompiler = defaultCompiler or os.environ.get('CC') or 'cc'
    defaultMajor = gccMajor(defaultCompiler)

    if defaultMajor is None or defaultMajor <= maxGcc:
        return None

    for major in range(maxGcc, _OLDEST_GCC - 1, -1):
        candidate = which('gcc-%d' % major)

        if candidate:
            return candidate

    return None


def cmakeToolchainFlags(cudaBin, which=shutil.which):
    """Extra cmake flags for the CUDA host compiler, as a string."""
    nvcc = os.path.join(cudaBin, 'nvcc') if cudaBin else None

    if not nvcc or not os.path.exists(nvcc):
        nvcc = which('nvcc')

    hostCompiler = cudaHostCompiler(nvcc, which=which)

    return '-DCUDA_HOST_COMPILER=%s' % hostCompiler if hostCompiler else ''


# A shell expression for the GUI flag: build it only when the Xft files
# RELION's FLTK needs are installed, and say so when they are not.
GUI_FLAG_SHELL = (
    '$(pkg-config --exists xft 2>/dev/null && echo -DGUI=ON || '
    '{ echo "RELION is built without its GUI: install the Xft, Xinerama, '
    'freetype and fontconfig development files to have it." >&2; '
    'echo -DGUI=OFF; })'
)
