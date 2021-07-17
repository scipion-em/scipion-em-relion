# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se)
# *
# * SciLifeLab, Stockholm University
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
from collections import OrderedDict
from emtable import Table

from pyworkflow.object import Float
from pwem.constants import ALIGN_PROJ, ALIGN_2D, ALIGN_NONE
from pwem.objects import Micrograph
import pwem.emlib.metadata as md
import pyworkflow.utils as pwutils

from .convert31 import OpticsGroups
from .convert_utils import relionToLocation


class FileTransform:
    """ Simple classes to handle input files.
    This classes will create a destination folder from which the input files
    will be accessible. The destination folder can be a symlink if the files
    are not required to be copied.
    """
    def __init__(self, srcPath, dstPath, copyFiles=False):
        self._src = srcPath
        self._dst = dstPath
        self._processed = set()

        if copyFiles:
            pwutils.makePath(dstPath)
            self.transform = self._copyFile
        else:
            pwutils.makePath(os.path.dirname(srcPath))
            pwutils.createLink(srcPath, dstPath)
            self.transform = self._addPrefix

    def _copyFile(self, inputFile):
        """ Copy file to the destination folder and return new filename. """
        if inputFile not in self._processed:
            pwutils.copyFile(os.path.join(self._src, inputFile), self._dst)
            self._processed.add(inputFile)

        return os.path.join(self._dst, os.path.basename(inputFile))

    def _addPrefix(self, inputFile):
        """ Just prepend the new symlink folder to the filename. """
        return os.path.join(self._dst, inputFile)


class RelionImport:
    """ Protocol to import existing Relion runs. """
    def __init__(self, protocol, starFile):
        self.protocol = protocol
        self._starFile = starFile
        self.copyOrLink = protocol.getCopyOrLink()
        self.version30 = False

    def importParticles(self):
        """ Import particles from 'run_data.star' """
        self.ignoreIds = self.protocol.ignoreIdColumn.get()
        self._imgDict = {}  # store which images stack have been linked/copied and the new path
        self._findImagesPath('rlnImageName')
        if self._micIdOrName:
            # If rlnMicrographName or rlnMicrographId then
            # create a set to link from particles
            self.micSet = self.protocol._createSetOfMicrographs()
            self.protocol.setSamplingRate(self.micSet)
            self.micSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
            self.protocol.fillAcquisition(self.micSet.getAcquisition())

        partSet = self.protocol._createSetOfParticles()
        partSet.setObjComment('Particles imported from Relion star file:\n%s' % self._starFile)

        # Update both samplingRate and acquisition with parameters
        # selected in the protocol form
        self.protocol.setSamplingRate(partSet)
        self._pixelSize = self.protocol.samplingRate.get()
        partSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
        self.protocol.fillAcquisition(partSet.getAcquisition())
        # Read the micrographs from the 'self._starFile' metadata
        # but fixing the filenames with new ones (linked or copy to extraDir)

        if self.version30:
            from .convert_deprecated import readSetOfParticles
            readSetOfParticles(
                self._starFile, partSet,
                preprocessImageRow=self._preprocessImageRow30,
                postprocessImageRow=self._postprocessImageRow30,
                readAcquisition=False, alignType=self.alignType)
        else:
            from relion.convert import readSetOfParticles
            readSetOfParticles(
                self._starFile, partSet,
                preprocessImageRow=None,
                postprocessImageRow=self._postprocessImageRow,
                readAcquisition=False, alignType=self.alignType,
                pixelSize=self._pixelSize)

        if self._micIdOrName:
            self.protocol._defineOutputs(outputMicrographs=self.micSet)
        self.protocol._defineOutputs(outputParticles=partSet)

        if self._classesFunc is not None:
            self._createClasses(partSet)

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesDict:
            index, fn, row = self._classesDict[classId]
            if fn.endswith('.mrc'):
                fn += ':mrc'  # Specify that are volumes to read them properly in xmipp
            item.getRepresentative().setLocation(index, fn)
            item._rlnclassDistribution = Float(row.get('rlnClassDistribution'))
            item._rlnAccuracyRotations = Float(row.get('rlnAccuracyRotations'))
            if self.version30:
                accInAngst = row.get('rlnAccuracyTranslations') * self._pixelSize
                item._rlnAccuracyTranslationsAngst = Float(accInAngst)
            else:
                item._rlnAccuracyTranslationsAngst = Float(row.get('rlnAccuracyTranslationsAngst'))

    def _createClasses(self, partSet):
        self._classesDict = {}  # store classes info, indexed by class id
        pathDict = {}

        self.protocol.info('Loading classes info from: %s' % self._modelStarFile)
        table = Table(fileName=self._modelStarFile, tableName='model_classes')
        for classNumber, row in enumerate(table):
            index, fn = relionToLocation(row.rlnReferenceImage)

            if fn in pathDict:
                newFn = pathDict.get(fn)
            else:
                clsPath = pwutils.findRootFrom(self._modelStarFile, fn)
                if clsPath is None:
                    newFn = fn
                else:
                    newFn = self.protocol._getExtraPath(os.path.basename(fn))
                    self.copyOrLink(os.path.join(clsPath, fn), newFn)
                pathDict[fn] = newFn

            self._classesDict[classNumber+1] = (index, newFn, row)

        clsSet = self._classesFunc(partSet)
        clsSet.classifyItems(updateClassCallback=self._updateClass)

        self.protocol._defineOutputs(outputClasses=clsSet)
        self.protocol._defineSourceRelation(partSet, clsSet)

    # -------------------------- INFO functions -------------------------------
    def validateParticles(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        self._findImagesPath("rlnImageName", warnings=False)

    def summaryParticles(self):
        """ Should be overwritten in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []

    def _getModelFile(self, dataStar):
        """ Retrieve the model star file from a given
        _data.star file.
        """
        modelStarFile = dataStar.replace('_data.star', '_model.star')

        if os.path.exists(modelStarFile):
            result = modelStarFile
        else:
            modelHalfStarFile = self._starFile.replace('_data.star',
                                                       '_half1_model.star')
            if os.path.exists(modelHalfStarFile):
                result = modelHalfStarFile
            else:
                result = None
                
        return result

    def _findImagesPath(self, label, warnings=True):
        # read the first table
        table = Table(fileName=self._starFile)
        acqRow = row = table[0]

        if row is None:
            raise Exception("Cannot import from empty metadata: %s"
                            % self._starFile)

        if not row.get('rlnOpticsGroup', False):
            self.version30 = True
            self.protocol.warning("Import from Relion version < 3.1 ...")
        else:
            acqRow = OpticsGroups.fromStar(self._starFile).first()
            # read particles table
            table = Table(fileName=self._starFile, tableName='particles')
            row = table[0]

        if not row.get(label, False):
            raise Exception("Label *%s* is missing in metadata: %s"
                            % (label, self._starFile))

        index, fn = relionToLocation(row.get(label))
        # Relion does not allow abs paths
        if fn.startswith("/"):
            raise Exception("ERROR: %s cannot be an absolute path: %s\n"
                            "Please create a symlink to an abs path instead."
                            % (label, fn))
        self._imgPath = pwutils.findRootFrom(self._starFile, fn)

        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: %s"
                                  % self._starFile)

        if (self._starFile.endswith('_data.star') and
                self._getModelFile(self._starFile)):
            self._modelStarFile = self._getModelFile(self._starFile)
            modelRow = Table(fileName=self._modelStarFile,
                             tableName='model_general')[0]
            classDimensionality = int(modelRow.rlnReferenceDimensionality)
            self._optimiserFile = self._starFile.replace('_data.star',
                                                         '_optimiser.star')
            if not os.path.exists(self._optimiserFile):
                autoRefine = int(modelRow.rlnNrClasses) == 1
            else:
                optimiserRow = Table(fileName=self._optimiserFile,
                                     tableName='optimiser_general')[0]
                autoRefine = optimiserRow.get('rlnModelStarFile2', False)

            self.alignType = ALIGN_PROJ

            if not autoRefine:
                if classDimensionality == 3:
                    self._classesFunc = self.protocol._createSetOfClasses3D
                else:
                    self._classesFunc = self.protocol._createSetOfClasses2D
                    self.alignType = ALIGN_2D
            else:
                self._classesFunc = None
        else:
            self._classesFunc = None
            self._modelStarFile = None
            modelRow = None
                       
            # Check if we have rot angle -> ALIGN_PROJ,
            # if only psi angle -> ALIGN_2D
            if (row.get('rlnAngleRot', False) and
                    float(row.rlnAngleRot) != 0.0):
                self.alignType = ALIGN_PROJ
            elif (row.get('rlnAnglePsi', False) and
                  float(row.rlnAnglePsi) != 0.0):
                self.alignType = ALIGN_2D
            else:
                self.alignType = ALIGN_NONE

        print("alignType: ", self.alignType)
            
        # Check if the MetaData contains either rlnMicrographName
        # or rlnMicrographId, this will be used when imported
        # particles to keep track of the particle's micrograph
        self._micIdOrName = (row.get('rlnMicrographName', False) or
                             row.get('rlnMicrographId', False))
        # init dictionary. It will be used in the preprocessing
        self.micDict = {}
        self._stackTrans = None
        self._micTrans = None

        return row, modelRow, acqRow

    def _preprocessImageRow30(self, img, imgRow):
        from .convert_deprecated import setupCTF, copyOrLinkFileName
        if self._imgPath is not None:
            copyOrLinkFileName(imgRow, self._imgPath, self.protocol._getExtraPath())
        setupCTF(imgRow, self.protocol.samplingRate.get())

        if self._micIdOrName:
            micId = imgRow.get('rlnMicrographId', None)
            micName = imgRow.get('rlnMicrographName', None)

            # Check which is the key to identify micrographs (id or name)
            if micId is not None:
                micKey = micId
            else:
                micKey = micName

            mic = self.micDict.get(micKey, None)

            # First time I found this micrograph (either by id or name)
            if mic is None:
                mic = Micrograph()
                mic.setObjId(micId)
                if micName is None:
                    micName = self.protocol._getExtraPath('fake_micrograph%6d' % micId)
                mic.setFileName(micName)
                mic.setMicName(os.path.basename(micName))
                self.micSet.append(mic)
                # Update dict with new Micrograph
                self.micDict[micKey] = mic

            # Update the row to set a MDL_MICROGRAPH_ID
            imgRow['rlnMicrographId'] = int(mic.getObjId())

    def _postprocessImageRow30(self, img, imgRow):
        if self.ignoreIds:
            img.setObjId(None)  # Force to generate a new id in Set

        if self._micIdOrName:
            micId = imgRow.get('rlnMicrographId', None)
            micName = imgRow.get('rlnMicrographName', None)
            if img.hasCoordinate():
                coord = img.getCoordinate()
                coord.setMicId(micId)
                coord.setMicName(os.path.basename(micName))

    def _postprocessImageRow(self, img, imgRow):
        # shortcut notation
        prot = self.protocol
        imgPath = self._imgPath

        if self.ignoreIds:
            img.setObjId(None)  # Force to generate a new id in Set

        if imgPath is not None:
            if self._stackTrans is None:
                self._stackTrans = FileTransform(imgPath,
                                                 prot._getExtraPath('Particles'),
                                                 prot.copyFiles)
            img.setFileName(self._stackTrans.transform(img.getFileName()))

        if self._micIdOrName:
            micId = imgRow.get('rlnMicrographId', None)
            micName = imgRow.get('rlnMicrographName', None)

            # Check which is the key to identify micrographs (id or name)
            if micId is not None:
                micKey = micId
            else:
                micKey = micName

            mic = self.micDict.get(micKey, None)

            # First time I found this micrograph (either by id or name)
            if mic is None:
                mic = Micrograph()
                mic.setObjId(micId)
                if micName is None:
                    micName = prot._getExtraPath('fake_micrograph%6d' % micId)
                else:
                    if not len(self.micDict):  # first time
                        if os.path.exists(os.path.join(imgPath, micName)):
                            micRoot = imgPath
                        else:
                            micRoot = pwutils.findRootFrom(self._starFile,
                                                           micName)
                        if micRoot is not None:
                            self._micTrans = FileTransform(
                                micRoot,
                                prot._getExtraPath('Micrographs'),
                                prot.copyFiles)
                    if self._micTrans is not None:
                        micName = self._micTrans.transform(micName)
                mic.setFileName(micName)
                mic.setMicName(os.path.basename(micName))
                mic.setAcquisition(img.getAcquisition())
                self.micSet.append(mic)
                # Update dict with new Micrograph
                self.micDict[micKey] = mic

            img.setMicId(mic.getObjId())

            if img.hasCoordinate():
                coord = img.getCoordinate()
                coord.setMicId(mic.getObjId())
                coord.setMicName(os.path.basename(micName))
    
    def loadAcquisitionInfo(self):
        """ Return a dictionary with acquisition values and 
        the sampling rate information.
        In the case of Relion, they are stored in the optics table
        """
        acquisitionDict = OrderedDict()

        try:
            _, modelRow, acqRow = self._findImagesPath('rlnImageName', warnings=False)

            if acqRow.get('rlnVoltage', False):
                acquisitionDict['voltage'] = acqRow.rlnVoltage

            if acqRow.get('rlnAmplitudeContrast', False):
                acquisitionDict['amplitudeContrast'] = acqRow.rlnAmplitudeContrast

            if acqRow.get('rlnSphericalAberration', False):
                acquisitionDict['sphericalAberration'] = acqRow.rlnSphericalAberration

            if modelRow is not None and modelRow.get('rlnPixelSize', False):
                acquisitionDict['samplingRate'] = modelRow.rlnPixelSize

        except Exception as ex:
            print("Error loading acquisition: ", str(ex))

        return acquisitionDict

    def importCoordinates(self, fileName, addCoordinate):
        from .convert_deprecated import rowToCoordinate
        for row in md.iterRows(fileName):
            coord = rowToCoordinate(row)
            addCoordinate(coord)
