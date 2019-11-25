# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] MRC Laboratory of Molecular Biology, MRC-LMB
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


from relion.constants import *
from . import metadata as md


class Writer:
    """ Helper class to convert from Scipion SetOfImages subclasses
    with star file format previous to Relion>3.1, but providing the same
     interface as the new Writer class.
    """
    def __init__(self, **kwargs):
        """
        Create a new instance with some configuration parameters.

        Keyword Args:
            useBaseName: (bool) By default the writer will use the id to
                generate shorter names. If this option is True, then
                the images base name will be used instead. This option
                might be useful in export protocols.

        """
        self._optics = None

    def writeSetOfMovies(self, moviesIterable, starFile):
        with open(starFile, 'w') as f:
            table = md.Table(columns=['rlnMicrographMovieName'])
            for img in moviesIterable:
                table.addRow(os.path.basename(img.getFileName()))
            table.writeStar(f)

    def writeSetOfMicrographs(self, micsIterable, starFile):
        raise Exception("Not implemented")

