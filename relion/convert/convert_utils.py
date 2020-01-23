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

import pyworkflow.utils as pwutils
import pwem


def convertBinaryFiles(imgSet, outputDir, extension='mrcs'):
    """ Convert binary images files to a format read by Relion.
    Or create links if there is no need to convert the binary files.

    Params:
        imgSet: input image set to be converted.
        outputDir: where to put the converted file(s)
    Return:
        A dictionary with old-file as key and new-file as value
        If empty, not conversion was done.
    """
    filesDict = {}
    ih = pwem.convert.ImageHandler()
    outputRoot = os.path.join(outputDir, 'input')
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
        newFn = getUniqueFileName(fn, 'stk')
        ih.convertStack(fn, newFn)
        print("   %s -> %s" % (newFn, fn))
        return newFn

    def replaceRoot(fn):
        """ Link create to the root folder, so just replace that
        in the name, no need to do anything else.
        """
        return fn.replace(rootDir, outputRoot)

    if ext == extension:
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