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

from .convert_utils import *
from .convert_deprecated import *
from .dataimport import *
import relion

# Writing of star files will be handle by the Writer class
# We have a new implementation of it for Relion > 3.1 since
# the star file format has changed in 3.
from . import convert30
from . import convert31


def createReader(**kwargs):
    """ Create a new Reader instance.
    By default it will create the version (3.1 or older) based on the current
    plugin binary. It can also be forced to use old format by passing
    the format='30' argument.
    """
    is30 = kwargs.get('format', '') == '30' or relion.Plugin.IS_30()
    Reader = convert30.Reader if is30 else convert31.Reader

    return Reader(**kwargs)


def writeSetOfParticles(imgSet, starFile, outputDir, **kwargs):
    """ Convenience function to a SetOfImages as Relion metadata using a Writer.

    Params:
        imgSet: the SetOfImages instance.
        starFile: the filename where to write the meta
        outputDir: where binary files will be converted or linked.
        filesMapping: this dict will help when there is need to replace images names

    Keyword Arguments:
        blockName: The name of the data block (default particles)
        fillMagnification: If True set magnification values (default False)
        alignType:
        fillRandomSubset:
        extraLabels:
        postprocessImageRow:
        format: string value to specify STAR format, if '30' it will use
            Relion3.0 format, if not, it will depends on the binary version
    """
    is30 = kwargs.get('format', '') == '30' or relion.Plugin.IS_30()
    Writer = convert30.Writer if is30 else convert31.Writer

    writer = Writer(outputDir=outputDir)
    return writer.writeSetOfParticles(imgSet, starFile, **kwargs)


def readSetOfParticles(starFile, partsSet, **kwargs):
    """ Convert a star file into a set of particles.

    Params:
        starFile: the filename of the star file
        partsSet: output particles set

    Keyword Arguments:
        blockName: The name of the data block (default particles)
        alignType:
        removeDisabled:
        format: string value to specify STAR format, if '30' it will use
            Relion3.0 format, if not, it will depends on the binary version
    """
    return createReader(**kwargs).readSetOfParticles(starFile, partsSet, **kwargs)
