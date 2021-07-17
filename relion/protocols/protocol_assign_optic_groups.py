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
import emtable

from pyworkflow.object import Integer
import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pwem.objects import SetOfMovies, SetOfParticles, SetOfMicrographs

from relion.convert.convert31 import OpticsGroups, getPixelSizeLabel
from .protocol_base import ProtRelionBase


class ProtRelionAssignOpticsGroup(ProtRelionBase):
    """ Assign Optics Group name and related parameters to an input set.
     Input set can be: movies, micrographs or particles.
    """
    _label = 'assign optics groups'
    _devStatus = PROD
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam,
                      pointerClass='SetOfMovies,SetOfMicrographs,SetOfParticles',
                      label="Input set", important=True,
                      help='Select the input set (Movies, Micrographs or '
                           'Particles) to assign Optics Group parameters.')

        form.addParam('inputType', params.EnumParam, default=0,
                      choices=['single group params', 'groups from star'],
                      display=params.EnumParam.DISPLAY_LIST,
                      label='Input Optics Groups from',
                      help='Select how to provide information about the optics '
                           'groups. In the case of single group provide '
                           'the parameters below, otherwise you need to provide '
                           'a star file.')

        group = form.addGroup("Optics Group params",
                              condition="inputType == 0")

        group.addParam('opticsGroupName', params.StringParam,
                       default='opticsGroup1',
                       label='Optics group name',
                       help='Relion-specific option. Name of this optics group. '
                            'Each group of movies with different '
                            'optics characteristics for CTF refinement '
                            'should have a unique name.')

        group.addParam('mtfFile', params.FileParam, allowsNull=True,
                       label='MTF-curve file',
                       help='User-provided STAR-file with the MTF-curve '
                            'of the detector. Use the wizard to load one '
                            'of the predefined ones provided at:\n'
                            '- [[https://www3.mrc-lmb.cam.ac.uk/relion/index.php/'
                            'FAQs#Where_can_I_find_MTF_curves_for_typical_detectors.3F]'
                            '[Relion\'s Wiki FAQs]]\n'
                            ' - [[https://www.gatan.com/techniques/cryo-em#MTF][Gatan\'s website]]\n\n'
                            'Relion param: *--mtf*')

        line = group.addLine('Beam tilt (mrad)',
                             help='Known beam tilt in the X/Y-direction (in mrad). '
                                  'Set to zero if unknown.')
        line.addParam('beamTiltX', params.FloatParam, default=0.,
                      label='X')
        line.addParam('beamTiltY', params.FloatParam, default=0.,
                      label='Y')

        group.addParam('gainFile', params.FileParam,
                       label='Gain reference',
                       help='A gain reference file is required for gain '
                            'correction')

        group.addParam('gainRot', params.EnumParam, default=0,
                       choices=['No rotation (0)',
                                ' 90 degrees (1)',
                                '180 degrees (2)',
                                '270 degrees (3)'],
                       label='Gain rotation',
                       help="Rotate the gain reference by this number times 90 "
                            "degrees clockwise in relion_display. This is the "
                            "same as -RotGain in MotionCor2. \n"
                            "Note that MotionCor2 uses a different convention "
                            "for rotation so it says 'counter-clockwise'.")

        group.addParam('gainFlip', params.EnumParam, default=0,
                       choices=['No flipping        (0)',
                                'Flip upside down   (1)',
                                'Flip left to right (2)'],
                       label='Gain flip',
                       help="Flip the gain reference after rotation. "
                            "This is the same as -FlipGain in MotionCor2. "
                            "0 means do nothing, 1 means flip Y (upside down) "
                            "and 2 means flip X (left to right).")

        group.addParam('defectFile', params.FileParam, allowsNull=True,
                       label='Defects file',
                       help='Location of a UCSF MotionCor2-style '
                            'defect text file or a defect map that '
                            'describe the defect pixels on the detector. '
                            'Each line of a defect text file should contain '
                            'four numbers specifying x, y, width and height '
                            'of a defect region. A defect map is an image '
                            '(MRC or TIFF), where 0 means good and 1 means '
                            'bad pixels. The coordinate system is the same '
                            'as the input movie before application of '
                            'binning, rotation and/or flipping.\n\n'
                            '_Note that the format of the defect text is '
                            'DIFFERENT from the defect text produced '
                            'by SerialEM!_\n One can convert a SerialEM-style '
                            'defect file into a defect map using IMOD '
                            'utilities e.g.:\n'
                            '*clip defect -D defect.txt -f tif movie.tif defect_map.tif*\n'
                            'See explanations in the SerialEM manual.\n'
                            'Leave empty if you do not have any defects, '
                            'or do not want to correct for defects on your detector.')

        form.addParam('inputStar', params.FileParam,
                      condition='inputType == 1',
                      label='Input Star file with optics groups',
                      help='Provide input star file with Optics groups '
                           'information. The input Star file should contain: \n'
                           '- *data_optics* table with values for each group.\n'
                           '- *data_micrographs* table with two columns: \n\n'
                           '\trlnMicrographName with the micName associated to '
                           'the input set.\n'
                           '\trlnOpticsGroup with the group number associated '
                           'to this micrograph.\n\nIf you provide rlnMicrographGainName '
                           'in the optics table, it has to point to a transformed '
                           'gain reference (rotated and flipped if necessary).')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep',
                                 self.inputSet.get().getObjId())

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self, inputId):
        inputSet = self.inputSet.get()

        getMicName = lambda item: item.getMicName()

        if isinstance(inputSet, SetOfMovies):
            outputSet = self._createSetOfMovies()
            outputName = 'outputMovies'
        elif isinstance(inputSet, SetOfMicrographs):
            outputSet = self._createSetOfMicrographs()
            outputName = 'outputMicrographs'
        elif isinstance(inputSet, SetOfParticles):
            outputSet = self._createSetOfParticles()
            outputName = 'outputParticles'
            getMicName = lambda item: item.getCoordinate().getMicName()
        else:
            raise Exception("Invalid input of type %s, expecting:\n"
                            "SetOfMovies, SetOfMicrographs or SetOfParticles"
                            % inputSet.getClassName())

        # Copy general info from input set
        outputSet.copyInfo(inputSet)

        if self.inputType == 0:  # single group params
            og = OpticsGroups.fromImages(inputSet)
            outputSet.copyItems(inputSet)

            if len(og) > 1:
                raise Exception("Multiple optics groups detected in the input!\n"
                                "Assigning single optics group params "
                                "is valid only when the input set contains "
                                "one optics group.")

            acq = inputSet.getAcquisition()
            params = {
                'rlnVoltage': acq.getVoltage(),
                'rlnSphericalAberration': acq.getSphericalAberration(),
                'rlnAmplitudeContrast': acq.getAmplitudeContrast(),
                getPixelSizeLabel(inputSet): inputSet.getSamplingRate(),
                'rlnImageSize': inputSet.getXDim(),
                'rlnOpticsGroupName': self.opticsGroupName.get(),
                'rlnBeamTiltX': self.beamTiltX.get(),
                'rlnBeamTiltY': self.beamTiltY.get()
            }
            og.updateAll(**params)

            if self.mtfFile.hasValue():
                inputMtf = self.mtfFile.get()
                outputMtf = self._getPath(os.path.basename(inputMtf))
                pwutils.copyFile(inputMtf, outputMtf)
                og.addColumns(rlnMtfFileName=outputMtf)
            if self.gainFile.hasValue():
                gainFn = self._convertGain()
                og.addColumns(rlnMicrographGainName=gainFn)
            if self.defectFile.hasValue():
                inputDef = self.defectFile.get()
                outputDef = self._getPath(os.path.basename(inputDef))
                pwutils.copyFile(inputDef, outputDef)
                og.addColumns(rlnMicrographDefectFile=outputDef)
        else:
            inputStar = self.inputStar.get()
            og = OpticsGroups.fromStar(inputStar)
            micTable = emtable.Table(fileName=inputStar,
                                     tableName='micrographs')
            micDict = {row.rlnMicrographName: row.rlnOpticsGroup
                       for row in micTable}

            # check if MTF file exists
            if og.hasColumn('rlnMtfFileName'):
                for i in og:
                    if not os.path.exists(i.rlnMtfFileName):
                        self.warning("MTF file %s not found for %s" % (
                            i.rlnMtfFileName, i.rlnOpticsGroupName
                        ))

            # check if gain file exists
            if og.hasColumn('rlnMicrographGainName'):
                for i in og:
                    if not pwutils.exists(i.rlnMicrographGainName):
                        self.warning("Gain reference file %s not found for %s" % (
                            i.rlnMicrographGainName, i.rlnOpticsGroupName
                        ))

            def updateItem(item, row):
                micName = getMicName(item)

                if micName not in micDict:
                    raise Exception("Micrograph name (aka micName) '%s' was "
                                    "not found in the 'data_micrographs' table of "
                                    "the input star file: %s"
                                    % (micName, inputStar))

                ogNumber = micDict[micName]

                if not hasattr(item, '_rlnOpticsGroup'):
                    item._rlnOpticsGroup = Integer()

                item._rlnOpticsGroup.set(ogNumber)

            outputSet.copyItems(inputSet,
                                updateItemCallback=updateItem,
                                doClone=False)

        print(og)

        og.toImages(outputSet)

        self._defineOutputs(**{outputName: outputSet})
        self._defineTransformRelation(inputSet, outputSet)
    
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validateMsgs = []

        if self.mtfFile.hasValue() and not os.path.exists(self.mtfFile.get()):
            validateMsgs.append("MTF file %s does not exist!" % self.mtfFile.get())

        if self.defectFile.hasValue() and not os.path.exists(self.defectFile.get()):
            validateMsgs.append("Defect file %s does not exist!" % self.defectFile.get())

        if self.gainRot.get() or self.gainFlip.get():
            try:
                from pwem import Domain
                eman2 = Domain.importFromPlugin('eman2', doRaise=True)
            except:
                validateMsgs.append("EMAN2 plugin not found!\nTransforming gain image "
                                    "requires EMAN2 plugin and binaries installed.")

        return validateMsgs
    
    def _summary(self):
        summary = []
        return summary
    
    def _methods(self):
        return []

    # -------------------------- UTILS functions ------------------------------
    def _convertGain(self):
        """ We need to transform gain file for a possible polishing job. """
        rotation = self.gainRot.get()
        flip = self.gainFlip.get()
        gainFn = self.gainFile.get()

        if rotation or flip:
            args = "%s %s " % (gainFn, self._getPath(os.path.basename(gainFn)))

            if flip:
                # flip axis Y - left to right
                args += "--process xform.flip:axis=%s " % ("y" if flip == 2 else "x")

            if rotation:
                args += "--rotate %d " % (rotation * 90)

            from pwem import Domain
            eman2 = Domain.importFromPlugin('eman2')
            pwutils.runJob(self._log, eman2.Plugin.getProgram('e2proc2d.py'), args,
                           env=eman2.Plugin.getEnviron())

            return self._getPath(os.path.basename(gainFn))
        else:
            outputGain = self._getPath(os.path.basename(gainFn))
            pwutils.createAbsLink(gainFn, outputGain)
            return outputGain
