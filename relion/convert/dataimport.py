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
from os.path import exists
from collections import OrderedDict

from pyworkflow.object import Float
from pwem.constants import ALIGN_PROJ, ALIGN_2D, ALIGN_NONE
from pwem.objects import Micrograph
import pwem.emlib.metadata as md
from pyworkflow.utils.path import findRootFrom

from .convert_utils import relionToLocation
from relion.convert.metadata import Table
#from .convert_deprecated import readSetOfParticles, rowToCoordinate


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
        self._findImagesPath(label='rlnImageName')
        if self._micIdOrName:
            # If MDL_MICROGRAPH_ID or MDL_MICROGRAPH then
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
                postprocessImageRow=None,
                readAcquisition=False, alignType=self.alignType)

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
            item._rlnclassDistribution = Float(row.getValue('rlnClassDistribution'))
            item._rlnAccuracyRotations = Float(row.getValue('rlnAccuracyRotations'))
            item._rlnAccuracyTranslations = Float(row.getValue('rlnAccuracyTranslations'))

    def _createClasses(self, partSet):
        self._classesDict = {}  # store classes info, indexed by class id
        pathDict = {}

        self.protocol.info('Loading classes info from: %s' % self._modelStarFile)
        modelMd = md.MetaData('model_classes@' + self._modelStarFile)
        for classNumber, objId in enumerate(modelMd):
            row = md.Row()
            row.readFromMd(modelMd, objId)
            index, fn = relionToLocation(row.getValue('rlnReferenceImage'))

            if fn in pathDict:
                newFn = pathDict.get(fn)
            else:
                clsPath = findRootFrom(self._modelStarFile, fn)
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
        errors = []
        try:
            self._findImagesPath(label="rlnImageName", warnings=False)
        except Exception as ex:
            errors.append(str(ex))

        return errors

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

        if exists(modelStarFile):
            result = modelStarFile
        else:
            modelHalfStarFile = self._starFile.replace('_data.star', '_half1_model.star')
            if exists(modelHalfStarFile):
                result = modelHalfStarFile
            else:
                result = None
                
        return result

    def _findImagesPath(self, label, warnings=True):
        # read the first table
        table = Table(fileName=self._starFile)
        row = table[0]
        acqRow = row

        if row is None:
            raise Exception("Cannot import from empty metadata: %s" % self._starFile)

        if not getattr(row, 'rlnOpticsGroup', False):
            self.version30 = True
            self.protocol.warning("Import from Relion version < 3.1 ...")
        else:
            acqRow = Table(fileName=self._starFile, tableName='optics')[0]
            # read particles table
            table = Table(fileName=self._starFile, tableName='particles')
            row = table[0]

        if not getattr(row, label, False):
            raise Exception("Label *%s* is missing in metadata: %s" % (label,
                                                                       self._starFile))

        index, fn = relionToLocation(getattr(row, label))
        self._imgPath = findRootFrom(self._starFile, fn)

        if warnings and self._imgPath is None:
            self.protocol.warning("Binary data was not found from metadata: %s" % self._starFile)

        if (self._starFile.endswith('_data.star') and
                self._getModelFile(self._starFile)):
            self._modelStarFile = self._getModelFile(self._starFile)
            modelRow = Table(fileName=self._modelStarFile,
                             tableName='model_general')[0]
            classDimensionality = int(modelRow.rlnReferenceDimensionality)
            self._optimiserFile = self._starFile.replace('_data.star',
                                                         '_optimiser.star')
            if not exists(self._optimiserFile):
                autoRefine = int(modelRow.rlnNrClasses) == 1
            else:
                optimiserRow = Table(fileName=self._optimiserFile,
                                     tableName='optimiser_general')[0]
                autoRefine = getattr(optimiserRow, 'rlnModelStarFile2', False)

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
            if (getattr(row, 'rlnAngleRot', False) and
                float(row.rlnAngleRot) != 0.0):
                self.alignType = ALIGN_PROJ
            elif (getattr(row, 'rlnAnglePsi', False) and
                  float(row.rlnAnglePsi) != 0.0):
                self.alignType = ALIGN_2D
            else:
                self.alignType = ALIGN_NONE
            
        # Check if the MetaData contains either MDL_MICROGRAPH_ID
        # or MDL_MICROGRAPH, this will be used when imported
        # particles to keep track of the particle's micrograph
        self._micIdOrName = (getattr(row, 'rlnMicrographName', False) or
                             getattr(row, 'rlnMicrographId', False))
        # init dictionary. It will be used in the preprocessing
        self.micDict = {}

        return row, modelRow, acqRow

    def _preprocessImageRow30(self, img, imgRow):
        from .convert_deprecated import setupCTF, copyOrLinkFileName
        if self._imgPath is not None:
            copyOrLinkFileName(imgRow, self._imgPath, self.protocol._getExtraPath())
        setupCTF(imgRow, self.protocol.samplingRate.get())

        if self._micIdOrName:
            micId = imgRow.getValue('rlnMicrographId', None)
            micName = imgRow.getValue('rlnMicrographName', None)

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
            imgRow.setValue('rlnMicrographId', int(mic.getObjId()))

    def _postprocessImageRow30(self, img, imgRow):
        if self.ignoreIds:
            img.setObjId(None)  # Force to generate a new id in Set

        if self._micIdOrName:
            micId = imgRow.getValue('rlnMicrographId', None)
            micName = imgRow.getValue('rlnMicrographName', None)
            if img.hasCoordinate():
                coord = img.getCoordinate()
                coord.setMicId(micId)
                coord.setMicName(os.path.basename(micName))
    
    def loadAcquisitionInfo(self):
        """ Return a dictionary with acquisition values and 
        the sampling rate information.
        In the case of Xmipp, they are stored in files:
        acquisition_info.xmd and microscope.xmd 
        """
        acquisitionDict = OrderedDict()

        try:
            _, modelRow, acqRow = self._findImagesPath(label='rlnImageName', warnings=False)

            if getattr(acqRow, 'rlnVoltage', False):
                acquisitionDict['voltage'] = acqRow.rlnVoltage

            if getattr(acqRow, 'rlnAmplitudeContrast', False):
                acquisitionDict['amplitudeContrast'] = acqRow.rlnAmplitudeContrast

            if getattr(acqRow, 'rlnSphericalAberration', False):
                acquisitionDict['sphericalAberration'] = acqRow.rlnSphericalAberration

            if (modelRow is not None and
                    getattr(modelRow, 'rlnPixelSize', False)):
                acquisitionDict['samplingRate'] = modelRow.rlnPixelSize

        except Exception as ex:
            print("Error loading acquisition: ", str(ex))

        return acquisitionDict

    def importCoordinates(self, fileName, addCoordinate):
        from .convert_deprecated import rowToCoordinate
        for row in md.iterRows(fileName):
            coord = rowToCoordinate(row)
            addCoordinate(coord)
