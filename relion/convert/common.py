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
"""
Common functions for conversions that will be used from both
newer Relion3.1 routines and old ones.
"""

import pyworkflow as pw
import pyworkflow.em.metadata as md


def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    if index != pw.em.NO_INDEX:
        return "%06d@%s" % (index, filename)

    return filename


def getImageLocation(location):
    return pw.em.ImageHandler.locationToXmipp(location)


def relionToLocation(filename):
    """ Return a location (index, filename) given
    a Relion filename with the index@filename structure. """
    if '@' in filename:
        indexStr, fn = filename.split('@')
        return int(indexStr), str(fn)
    else:
        return pw.em.NO_INDEX, str(filename)


def setRelionAttributes(obj, objRow, *labels):
    """ Set an attribute to obj from a label that is not
    basic ones. The new attribute will be named _rlnLabelName
    and the datatype will be set correctly.
    """
    for label in labels:
        setattr(obj, '_%s' % md.label2Str(label),
                objRow.getValueAsObject(label))
