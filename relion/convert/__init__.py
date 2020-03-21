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
# the star file format has changed in 3.1
if relion.Plugin.IS_GT30():
    from .convert31 import Writer, Reader
else:
    from .convert30 import Writer, Reader


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
    """
    writer = Writer(outputDir=outputDir)
    return writer.writeSetOfParticles(imgSet, starFile, **kwargs)
