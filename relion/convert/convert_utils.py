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
"""
Utility functions for conversions that will be used from both
newer Relion3.1 routines and old ones.
"""

import os

import pyworkflow.utils as pwutils
import pwem
from pwem.emlib.image import ImageHandler


def locationToRelion(index, filename):
    """ Convert an index and filename location
    to a string with @ as expected in Relion.
    """
    if index != pwem.NO_INDEX:
        return "%06d@%s" % (index, filename)

    return filename


def getImageLocation(location):
    return ImageHandler.locationToXmipp(location)


def relionToLocation(filename):
    """ Return a location (index, filename) given
    a Relion filename with the index@filename structure. """
    if '@' in filename:
        indexStr, fn = filename.split('@')
        return int(indexStr), str(fn)
    else:
        return pwem.NO_INDEX, str(filename)


def setRelionAttributes(obj, objRow, *labels):
    """ Set an attribute to obj from a label that is not
    basic ones. The new attribute will be named _rlnLabelName
    and the datatype will be set correctly.
    """
    # FIXME: Remove dependency from md (and thus from Xmipp)
    import pwem.emlib.metadata as md
    for label in labels:
        setattr(obj, '_%s' % md.label2Str(label),
                objRow.getValueAsObject(label))


def convertBinaryFiles(imgSet, outputDir, extension='mrcs', forceConvert=False):
    """ Convert binary images files to a format read by Relion.
    Or create links if there is no need to convert the binary files.

    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
        extension: extension accepted by the program
        forceConvert: if True, the files will be converted and no root will be used
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesDict = {}
    ih = ImageHandler()
    outputRoot = outputDir if forceConvert else os.path.join(outputDir, 'input')
    # Get the extension without the dot
    stackFiles = imgSet.getFiles()
    ext = pwutils.getExt(next(iter(stackFiles)))[1:]
    rootDir = pwutils.commonPath(list(stackFiles))

    def getUniqueFileName(fn, extension):
        """ Get an unique file for either link or convert files.
        It is possible that the base name overlap if they come
        from different runs. (like particles.mrcs after relion preprocess)
        """
        newFn = os.path.join(outputRoot, pwutils.replaceBaseExt(fn, extension))
        newRoot = pwutils.removeExt(newFn)

        values = filesDict.values()
        counter = 1

        while newFn in values:
            counter += 1
            newFn = '%s_%05d.%s' % (newRoot, counter, extension)

        return newFn

    def createBinaryLink(fn):
        """ Just create a link named .mrcs to Relion understand
        that it is a binary stack file and not a volume.
        """
        newFn = getUniqueFileName(fn, extension)
        if not os.path.exists(newFn):
            pwutils.createLink(fn, newFn)
            print("   %s -> %s" % (newFn, fn))
        return newFn

    def convertStack(fn):
        """ Convert from a format that is not read by Relion
        to an spider stack.
        """
        newFn = getUniqueFileName(fn, 'mrcs')
        ih.convertStack(fn, newFn)
        print("   %s -> %s" % (newFn, fn))
        return newFn

    def replaceRoot(fn):
        """ Link create to the root folder, so just replace that
        in the name, no need to do anything else.
        """
        return fn.replace(rootDir, outputRoot)

    if forceConvert:
        mapFunc = convertStack
    elif ext == extension:
        print("convertBinaryFiles: creating soft links.")
        print("   Root: %s -> %s" % (outputRoot, rootDir))
        mapFunc = replaceRoot
        # FIXME: There is a bug in pwutils.createLink when input is a single folder
        # pwutils.createLink(rootDir, outputRoot)
        # relativeOutput = os.path.join(os.path.relpath(rootDir, outputRoot), rootDir)
        # If the rootDir is a prefix in the outputRoot (usually Runs)
        # we need to prepend that basename to make the link works
        if rootDir in outputRoot:
            relativeOutput = os.path.join(os.path.relpath(rootDir, outputRoot),
                                          os.path.basename(rootDir))
        else:
            relativeOutput = os.path.relpath(rootDir,
                                             os.path.dirname(outputRoot))
        if not os.path.exists(outputRoot):
            os.symlink(relativeOutput, outputRoot)
    elif ext == 'mrc' and extension == 'mrcs':
        print("convertBinaryFiles: creating soft links (mrcs -> mrc).")
        mapFunc = createBinaryLink
    elif ext.endswith('hdf'):  # assume eman .hdf format
        print("convertBinaryFiles: converting stacks. (%s -> %s)"
              % (extension, ext))
        mapFunc = convertStack
    else:
        mapFunc = None

    if mapFunc is not None:
        pwutils.makePath(outputRoot)
        for fn in stackFiles:
            newFn = mapFunc(fn)  # convert or link
            filesDict[fn] = newFn  # map new filename

    return filesDict


def relativeFromFileName(imgRow, prefixPath):
    """ Remove some prefix from filename in row. """
    index, imgPath = relionToLocation(imgRow['rlnImageName'])
    newImgPath = os.path.relpath(imgPath, prefixPath)
    imgRow['rlnImageName'] = locationToRelion(index, newImgPath)
