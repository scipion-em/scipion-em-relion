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
from .convert_coordinates import *
from .dataimport import *
import relion
import math

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


def createWriter(**kwargs):
    """ Create a new Writer instance.
    By default it will create the version (3.1 or older) based on the current
    plugin binary. It can also be forced to use old format by passing
    the format='30' argument.
    """
    is30 = kwargs.get('format', '') == '30' or relion.Plugin.IS_30()
    Writer = convert30.Writer if is30 else convert31.Writer

    return Writer(**kwargs)


def writeSetOfParticles(imgSet, starFile, **kwargs):
    """ Convenience function to a SetOfImages as Relion metadata using a Writer.

    Params:
        imgSet: the SetOfImages instance.
        starFile: the filename where to write the meta

    Keyword Arguments:
        outputDir: where binary files will be converted or linked.
        blockName: The name of the data block (default particles)
        fillMagnification: If True set magnification values (default False)
        alignType:
        fillRandomSubset:
        extraLabels:
        postprocessImageRow:
        format: string value to specify STAR format, if '30' it will use
            Relion3.0 format, if not, it will depends on the binary version
    """
    return createWriter(**kwargs).writeSetOfParticles(imgSet, starFile, **kwargs)


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


class ClassesLoader:
    """ Helper class to read classes information from star files produced
    by Relion classification runs (2D or 3D).
    """
    def __init__(self, protocol, alignType):
        self._protocol = protocol
        self._alignType = alignType

    def _loadClassesInfo(self, iteration):
        """ Read some information about the produced Relion 3D classes
        from the *model.star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        modelFn = self._protocol._getFileName('model', iter=iteration)
        modelIter = Table.iterRows('model_classes@' + modelFn)

        for classNumber, row in enumerate(modelIter):
            index, fn = relionToLocation(row.rlnReferenceImage)
            # Store info indexed by id
            self._classesInfo[classNumber + 1] = (index, fn, row)

    def fillClassesFromIter(self, clsSet, iteration):
        """ Create the SetOfClasses3D from a given iteration. """
        prot = self._protocol  # shortcut
        self._loadClassesInfo(iteration)

        tableName = 'particles@' if Plugin.IS_GT30() else ''
        dataStar = prot._getFileName('data', iter=iteration)

        pixelSize = prot.inputParticles.get().getSamplingRate()
        self._reader = createReader(alignType=self._alignType,
                                    pixelSize=pixelSize)

        mdIter = Table.iterRows(tableName + dataStar, key='rlnImageId')
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=mdIter,
                             doClone=False)

    def _updateParticle(self, item, row):
        item.setClassId(row.rlnClassNumber)
        self._reader.setParticleTransform(item, row)

        # Try to create extra objects only once if item is reused
        if not hasattr(item, '_rlnNormCorrection'):
            item._rlnNormCorrection = Float()
            item._rlnLogLikeliContribution = Float()
            item._rlnMaxValueProbDistribution = Float()

        item._rlnNormCorrection.set(row.rlnNormCorrection)
        item._rlnLogLikeliContribution.set(row.rlnLogLikeliContribution)
        item._rlnMaxValueProbDistribution.set(row.rlnMaxValueProbDistribution)

        if hasattr(item, '_rlnGroupName'):
            item._rlnGroupName.set(row.rlnGroupName)
        elif hasattr(row, 'rlnGroupName'):
            item._rlnGroupName = String(row.rlnGroupName)

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            item.setAlignment(self._alignType)
            if self._alignType == pwem.ALIGN_PROJ:
                fn += ':mrc'  # mark reference as a MRC volume
            item.getRepresentative().setLocation(index, fn)
            item._rlnClassDistribution = Float(row.rlnClassDistribution)
            item._rlnAccuracyRotations = Float(row.rlnAccuracyRotations)
            if Plugin.IS_GT30():
                item._rlnAccuracyTranslationsAngst = Float(row.rlnAccuracyTranslationsAngst)
            else:
                item._rlnAccuracyTranslations = Float(row.rlnAccuracyTranslations)


class DefocusGroups:
    """ Helper class to create defocus groups for particles. """
    class Group:
        """ Single CTF group. """
        def __init__(self, id):
            self.id = id
            self.count = 0
            self.minDefocus = math.inf
            self.maxDefocus = -math.inf

        def addDefocus(self, defocus):
            self.count += 1
            if defocus < self.minDefocus:
                self.minDefocus = defocus

            if defocus > self.maxDefocus:
                self.maxDefocus = defocus

    def __initGroups(self):
        self._groups = []

    def __addGroup(self):
        group = self.Group(len(self._groups) + 1)
        self._groups.append(group)
        return group

    def __init__(self):
        self._groups = []

    def __len__(self):
        return len(self._groups)

    def __iter__(self):
        return iter(self._groups)

    def __str__(self):
        s = ">>> Defocus groups: %d\n" % len(self)
        row_format = u"{:>15}{:>15}{:>10}\n"
        s += row_format.format("Min (A)", "Max (A)", "Count")

        for group in self._groups:
            s += row_format.format("%0.3f" % group.minDefocus,
                                   "%0.3f" % group.maxDefocus,
                                   group.count)
        return s

    def splitByDiff(self, inputParts, defocusDiff=1000, minGroupSize=10):
        self.__initGroups()
        group = self.__addGroup()

        for part in inputParts.iterItems(orderBy=['_ctfModel._defocusU']):
            defocus = part.getCTF().getDefocusU()
            # Only when we reach the min number of particles
            # and the defocus difference, we create a new group
            if (group.count >= minGroupSize
                    and (defocus - group.minDefocus > defocusDiff)):
                group = self.__addGroup()

            group.addDefocus(defocus)

    def getGroup(self, defocus):
        """ Return the group that this defocus belong. """
        if (defocus < self._groups[0].minDefocus
                or defocus > self._groups[-1].maxDefocus):
            return None

        for group in self._groups:
            if defocus <= group.maxDefocus:
                return group
