# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              Laura del Cano (ldelcano@cnb.csic.es) [2]
# *              Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [3]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [3] MRC Laboratory of Molecular Biology, MRC-LMB
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

from os.path import basename
import numpy as np

from pyworkflow.object import ObjectWrap, String, Integer
from pwem.constants import ALIGN_2D, ALIGN_3D, ALIGN_PROJ, ALIGN_NONE
import pwem.convert.transformations as tfs
import pwem.objects as pwobj

from relion.constants import *
from .convert_utils import *


def objectToRow(obj, row, attrDict, extraLabels=[]):
    """ This function will convert an EMObject into a XmippMdRow.
    Params:
        obj: the EMObject instance (input)
        row: the XmippMdRow instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).        
        extraLabels: a list with extra labels that could be included
            as _xmipp_labelName
    """
    row.set(md.RLN_IMAGE_ENABLED, obj.isEnabled())

    for attr, label in attrDict.items():
        if hasattr(obj, attr):
            valueType = md.label2Python(label)
            row.set(label, valueType(getattr(obj, attr).get()))

    attrLabels = attrDict.values()

    for label in extraLabels:
        attrName = '_' + md.label2Str(label)
        if label not in attrLabels and hasattr(obj, attrName):
            value = obj.getAttributeValue(attrName)
            row.set(label, value)


def rowToObject(row, obj, attrDict, extraLabels=[]):
    """ This function will convert from a XmippMdRow to an EMObject.
    Params:
        row: the XmippMdRow instance (input)
        obj: the EMObject instance (output)
        attrDict: dictionary with the map between obj attributes(keys) and 
            row MDLabels in Xmipp (values).
        extraLabels: a list with extra labels that could be included
            as properties with the label name such as: _rlnSomeThing
    """
    obj.setEnabled(row.getValue(md.RLN_IMAGE_ENABLED, 1) > 0)

    for attr, label in attrDict.items():
        value = row.getValue(label)
        if not hasattr(obj, attr):
            setattr(obj, attr, ObjectWrap(value))
        else:
            getattr(obj, attr).set(value)

    attrLabels = attrDict.values()

    for label in extraLabels:
        if label not in attrLabels and row.hasLabel(label):
            labelStr = md.label2Str(label)
            setattr(obj, '_' + labelStr, row.getValueAsObject(label))


def setObjId(obj, mdRow, label=md.RLN_IMAGE_ID):
    obj.setObjId(mdRow.getValue(label, None))


def setRowId(mdRow, obj, label=md.RLN_IMAGE_ID):
    mdRow.set(label, int(obj.getObjId()))


def acquisitionToRow(acquisition, ctfRow):
    """ Set labels values from acquisition to md row. """
    objectToRow(acquisition, ctfRow, ACQUISITION_DICT)


def rowToAcquisition(acquisitionRow):
    """ Create an acquisition from a row of a meta """
    if acquisitionRow.containsAll(ACQUISITION_DICT):
        acquisition = pwobj.Acquisition()
        rowToObject(acquisitionRow, acquisition, ACQUISITION_DICT)
    else:
        acquisition = None

    return acquisition


def setPsdFiles(ctfModel, ctfRow):
    """ Set the PSD files of CTF estimation related
    to this ctfModel. The values will be read from
    the ctfRow if present.
    """
    for attr, label in CTF_PSD_DICT.items():
        if ctfRow.containsLabel(label):
            setattr(ctfModel, attr, String(ctfRow.getValue(label)))


def ctfModelToRow(ctfModel, ctfRow):
    """ Set labels values from ctfModel to md row. """
    # Refresh phase shift!
    phaseShift = ctfModel.getPhaseShift()

    if phaseShift is not None:
        ctfRow.set(md.RLN_CTF_PHASESHIFT, phaseShift)

    objectToRow(ctfModel, ctfRow, CTF_DICT, extraLabels=CTF_EXTRA_LABELS)


def rowToCtfModel(ctfRow):
    """ Create a CTFModel from a row of a meta """
    if ctfRow.containsAll(CTF_DICT):
        ctfModel = pwobj.CTFModel()

        rowToObject(ctfRow, ctfModel, CTF_DICT, extraLabels=CTF_EXTRA_LABELS)
        if ctfRow.hasLabel(md.RLN_CTF_PHASESHIFT):
            ctfModel.setPhaseShift(ctfRow.getValue(md.RLN_CTF_PHASESHIFT, 0))
        ctfModel.standardize()
        setPsdFiles(ctfModel, ctfRow)
    else:
        ctfModel = None

    return ctfModel


def geometryFromMatrix(matrix, inverseTransform):
    from pwem.convert.transformations import translation_from_matrix, euler_from_matrix

    if inverseTransform:
        matrix = np.linalg.inv(matrix)
        shifts = -translation_from_matrix(matrix)
    else:
        shifts = translation_from_matrix(matrix)
    angles = -np.rad2deg(euler_from_matrix(matrix, axes='szyz'))
    return shifts, angles


def matrixFromGeometry(shifts, angles, inverseTransform):
    """ Create the transformation matrix from a given
    2D shifts in X and Y...and the 3 euler angles.
    """
    radAngles = -np.deg2rad(angles)
    M = tfs.euler_matrix(radAngles[0], radAngles[1], radAngles[2], 'szyz')
    if inverseTransform:
        M[:3, 3] = -shifts[:3]
        M = np.linalg.inv(M)
    else:
        M[:3, 3] = shifts[:3]

    return M


def alignmentToRow(alignment, alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
                          -> for xmipp implies alignment
    """
    is2D = alignType == ALIGN_2D
    is3D = alignType == ALIGN_3D
    inverseTransform = alignType == ALIGN_PROJ
    matrix = alignment.getMatrix()
    shifts, angles = geometryFromMatrix(matrix, inverseTransform)

    alignmentRow.set(md.RLN_ORIENT_ORIGIN_X, shifts[0])
    alignmentRow.set(md.RLN_ORIENT_ORIGIN_Y, shifts[1])

    if is2D:
        angle = angles[0] + angles[2]
        alignmentRow.set(md.RLN_ORIENT_PSI, -angle)

        flip = bool(np.linalg.det(matrix[0:2, 0:2]) < 0)
        if flip:
            print("FLIP in 2D not implemented")
    elif is3D:
        raise Exception("3D alignment conversion for Relion not implemented. "
                        "It seems the particles were generated with an "
                        "incorrect alignment type. You may either re-launch "
                        "the protocol that generates the particles "
                        "with angles or set 'Consider previous alignment?' "
                        "to No")
    else:
        alignmentRow.set(md.RLN_ORIENT_ORIGIN_Z, shifts[2])
        alignmentRow.set(md.RLN_ORIENT_ROT, angles[0])
        alignmentRow.set(md.RLN_ORIENT_TILT, angles[1])
        alignmentRow.set(md.RLN_ORIENT_PSI, angles[2])


def rowToAlignment(alignmentRow, alignType):
    """
    is2D == True-> matrix is 2D (2D images alignment)
            otherwise matrix is 3D (3D volume alignment or projection)
    invTransform == True  -> for xmipp implies projection
    """
    if alignType == ALIGN_3D:
        raise Exception("3D alignment conversion for Relion not implemented.")

    is2D = alignType == ALIGN_2D
    inverseTransform = alignType == ALIGN_PROJ
    if alignmentRow.containsAny(ALIGNMENT_DICT):
        alignment = pwobj.Transform()
        angles = np.zeros(3)
        shifts = np.zeros(3)
        shifts[0] = alignmentRow.get(md.RLN_ORIENT_ORIGIN_X, 0.)
        shifts[1] = alignmentRow.get(md.RLN_ORIENT_ORIGIN_Y, 0.)
        if not is2D:
            angles[0] = alignmentRow.get(md.RLN_ORIENT_ROT, 0.)
            angles[1] = alignmentRow.get(md.RLN_ORIENT_TILT, 0.)
            angles[2] = alignmentRow.get(md.RLN_ORIENT_PSI, 0.)
            shifts[2] = alignmentRow.get(md.RLN_ORIENT_ORIGIN_Z, 0.)
        else:
            angles[2] = - alignmentRow.get(md.RLN_ORIENT_PSI, 0.)
        M = matrixFromGeometry(shifts, angles, inverseTransform)
        alignment.setMatrix(M)
    else:
        alignment = None

    return alignment


def coordinateToRow(coord, coordRow, copyId=True):
    """ Set labels values from Coordinate coord to md row. """
    if copyId:
        setRowId(coordRow, coord)
    objectToRow(coord, coordRow, COOR_DICT, extraLabels=COOR_EXTRA_LABELS)
    if coord.getMicName():
        micName = coord.getMicName()
        coordRow.set(md.RLN_MICROGRAPH_NAME, str(micName.replace(" ", "")))
    else:
        if coord.getMicId():
            coordRow.set(md.RLN_MICROGRAPH_NAME, str(coord.getMicId()))


def rowToCoordinate(coordRow):
    """ Create a Coordinate from a row of a meta """
    # Check that all required labels are present in the row
    if coordRow.containsAll(COOR_DICT):
        coord = pwobj.Coordinate()
        rowToObject(coordRow, coord, COOR_DICT, extraLabels=COOR_EXTRA_LABELS)

        micName = None

        if coordRow.hasLabel(md.RLN_MICROGRAPH_ID):
            micId = int(coordRow.get(md.RLN_MICROGRAPH_ID))
            coord.setMicId(micId)
            # If RLN_MICROGRAPH_NAME is not present, use the id as a name
            micName = micId

        if coordRow.hasLabel(md.RLN_MICROGRAPH_NAME):
            micName = coordRow.get(md.RLN_MICROGRAPH_NAME)

        coord.setMicName(micName)

    else:
        coord = None

    return coord


def readSetOfCoordinates(coordSet, coordFiles, micList=None):
    """ Read a set of coordinates from given coordinate files
    associated to some SetOfMicrographs.
    Params:
        micSet and coordFiles should have same length and same order.
        coordSet: empty SetOfCoordinates to be populated.
    """
    if micList is None:
        micList = coordSet.getMicrographs()

    for mic, coordFn in zip(micList, coordFiles):

        if not os.path.exists(coordFn):
            print("WARNING: Missing coordinates star file: ", coordFn)

        try:
            readCoordinates(mic, coordFn, coordSet)
        except Exception:
            print("WARNING: Error reading coordinates star file: ", coordFn)


def readCoordinates(mic, fileName, coordsSet):
    for row in md.iterRows(fileName):
        coord = rowToCoordinate(row)
        coord.setX(coord.getX())
        coord.setY(coord.getY())
        coord.setMicrograph(mic)
        coordsSet.append(coord)


def imageToRow(img, imgRow, imgLabel=md.RLN_IMAGE_NAME, **kwargs):
    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, imgRow)

    setRowId(imgRow, img)  # Set the id in the metadata as MDL_ITEM_ID
    index, fn = img.getLocation()
    # check if the is a file mapping
    filesDict = kwargs.get('filesDict', {})
    filename = filesDict.get(fn, fn)

    imgRow.set(imgLabel, locationToRelion(index, filename))

    if kwargs.get('writeCtf', True) and img.hasCTF():
        ctfModelToRow(img.getCTF(), imgRow)

    # alignment is mandatory at this point, it should be check
    # and detected defaults if not passed at readSetOf.. level
    alignType = kwargs.get('alignType')

    if alignType != ALIGN_NONE and img.hasTransform():
        alignmentToRow(img.getTransform(), imgRow, alignType)

    if kwargs.get('writeAcquisition', True) and img.hasAcquisition():
        acquisitionToRow(img.getAcquisition(), imgRow)

    # Write all extra labels to the row    
    objectToRow(img, imgRow, {},
                extraLabels=IMAGE_EXTRA_LABELS + kwargs.get('extraLabels', []))

    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, imgRow)


def volumeToRow(vol, volRow, **kwargs):
    """ Set labels values from Micrograph mic to md row. """
    imageToRow(vol, volRow, writeAcquisition=False, **kwargs)


def rowToVolume(volRow, **kwargs):
    """ Create a Volume object from a row of metadata. """
    return rowToParticle(volRow, particleClass=pwobj.Volume, **kwargs)


def particleToRow(part, partRow, **kwargs):
    """ Set labels values from Particle to md row. """
    coord = part.getCoordinate()
    if coord is not None:
        coordinateToRow(coord, partRow, copyId=False)
    if part.hasMicId():
        partRow.set(md.RLN_MICROGRAPH_ID, int(part.getMicId()))
        # If the row does not contains the micrograph name
        # use a fake micrograph name using id to relion
        # could at least group for CTF using that
        if not partRow.hasLabel(md.RLN_MICROGRAPH_NAME):
            partRow.set(md.RLN_MICROGRAPH_NAME,
                        'fake_micrograph_%06d.mrc' % part.getMicId())
    if part.hasAttribute('_rlnParticleId'):
        partRow.set(md.RLN_PARTICLE_ID, int(part._rlnParticleId.get()))

    if part.hasAttribute('_rlnRandomSubset'):
        partRow.set(md.RLN_PARTICLE_RANDOM_SUBSET,
                    int(part._rlnRandomSubset.get()))
        if part.hasAttribute('_rlnBeamTiltX'):
            partRow.set('rlnBeamTiltX',
                        float(part._rlnBeamTiltX.get()))
            partRow.set('rlnBeamTiltY',
                        float(part._rlnBeamTiltY.get()))

    imageToRow(part, partRow, md.RLN_IMAGE_NAME, **kwargs)


def rowToParticle(partRow, particleClass=pwobj.Particle, **kwargs):
    """ Create a Particle from a row of a meta """
    img = particleClass()

    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    preprocessImageRow = kwargs.get('preprocessImageRow', None)
    if preprocessImageRow:
        preprocessImageRow(img, partRow)

    # Decompose Relion filename
    index, filename = relionToLocation(partRow.get(md.RLN_IMAGE_NAME))
    img.setLocation(index, filename)

    if partRow.containsLabel(md.RLN_PARTICLE_CLASS):
        img.setClassId(partRow.get(md.RLN_PARTICLE_CLASS))

    if kwargs.get('readCtf', True):
        img.setCTF(rowToCtfModel(partRow))

    # alignment is mandatory at this point, it should be check
    # and detected defaults if not passed at readSetOf.. level
    alignType = kwargs.get('alignType')

    if alignType != ALIGN_NONE:
        img.setTransform(rowToAlignment(partRow, alignType))

    if kwargs.get('readAcquisition', True):
        img.setAcquisition(rowToAcquisition(partRow))

    if kwargs.get('magnification', None):
        img.getAcquisition().setMagnification(kwargs.get("magnification"))

    setObjId(img, partRow)
    # Read some extra labels
    rowToObject(partRow, img, {},
                extraLabels=IMAGE_EXTRA_LABELS + kwargs.get('extraLabels', []))

    img.setCoordinate(rowToCoordinate(partRow))

    # copy micId if available from row to particle
    if partRow.hasLabel(md.RLN_MICROGRAPH_ID):
        img.setMicId(partRow.get(md.RLN_MICROGRAPH_ID))

    # copy particleId if available from row to particle
    if partRow.hasLabel(md.RLN_PARTICLE_ID):
        img._rlnParticleId = Integer(partRow.get(md.RLN_PARTICLE_ID))

    # Provide a hook to be used if something is needed to be 
    # done for special cases before converting image to row
    postprocessImageRow = kwargs.get('postprocessImageRow', None)
    if postprocessImageRow:
        postprocessImageRow(img, partRow)
    return img


def readSetOfParticles(filename, imgSet, rowToFunc=rowToParticle, **kwargs):
    """read from Relion image meta
        filename: The metadata filename where the image are.
        imgSet: the SetOfParticles that will be populated.
        rowToParticle: this function will be used to convert the row to Object
    """
    imgMd = md.MetaData(filename)
    # By default remove disabled items from metadata
    # be careful if you need to preserve the original number of items
    if kwargs.get('removeDisabled', True):
        imgMd.removeDisabled()

    for imgRow in md.iterRows(imgMd):
        img = rowToFunc(imgRow, **kwargs)
        imgSet.append(img)

    imgSet.setHasCTF(img.hasCTF())
    imgSet.setAlignment(kwargs['alignType'])


def setOfImagesToMd(imgSet, imgMd, imgToFunc, **kwargs):
    """ This function will fill Relion metadata from a SetOfMicrographs
    Params:
        imgSet: the set of images to be converted to metadata
        md: metadata to be filled
        rowFunc: this function can be used to setup the row before 
            adding to meta
    """

    if 'alignType' not in kwargs:
        kwargs['alignType'] = imgSet.getAlignment()

    for img in imgSet:
        objId = imgMd.addObject()
        imgRow = md.Row()
        imgToFunc(img, imgRow, **kwargs)
        imgRow.writeToMd(imgMd, objId)


def writeSetOfImages(imgSet, filename, imgToFunc,
                     blockName='Images', **kwargs):
    """ This function will write a SetOfImages as metadata.
    Params:
        imgSet: the set of images to be written (particles,
        micrographs or volumes)
        filename: the filename where to write the metadata.`
        rowFunc: this function can be used to setup the row before
            adding to metadata.
    """
    mdFn = md.MetaData()

    setOfImagesToMd(imgSet, mdFn, imgToFunc, **kwargs)
    mdFn.write('%s@%s' % (blockName, filename))


def writeSetOfVolumes(volSet, filename, blockName='Volumes', **kwargs):
    writeSetOfImages(volSet, filename, volumeToRow, blockName, **kwargs)


def writeReferences(inputSet, outputRoot, useBasename=False, **kwargs):
    """
    Write references star and stack files from SetOfAverages or SetOfClasses2D/3D.
    Params:
        inputSet: the input SetOfParticles to be converted
        outputRoot: where to write the output files.
        basename: If True, use the basename of the stack for setting path.
    """
    refsMd = md.MetaData()
    stackFile = outputRoot + '.stk'
    stackName = basename(stackFile) if useBasename else stackFile
    starFile = outputRoot + '.star'
    ih = ImageHandler()
    row = md.Row()

    def _convert(item, i, convertFunc):
        index = i + 1
        ih.convert(item, (index, stackFile))
        item.setLocation(index, stackName)
        convertFunc(item, row, **kwargs)
        row.writeToMd(refsMd, refsMd.addObject())

    if isinstance(inputSet, pwobj.SetOfAverages):
        for i, img in enumerate(inputSet):
            _convert(img, i, particleToRow)

    elif isinstance(inputSet, pwobj.SetOfClasses2D):
        for i, rep in enumerate(inputSet.iterRepresentatives()):
            _convert(rep, i, particleToRow)

    elif isinstance(inputSet, pwobj.SetOfClasses3D):
        for i, rep in enumerate(inputSet.iterRepresentatives()):
            _convert(rep, i, imageToRow)

    elif isinstance(inputSet, pwobj.SetOfVolumes):
        for i, vol in enumerate(inputSet):
            _convert(vol, i, imageToRow)

    elif isinstance(inputSet, pwobj.Volume):
        _convert(inputSet, 0, imageToRow)

    else:
        raise Exception('Invalid object type: %s' % type(inputSet))

    refsMd.write(starFile)


def copyOrLinkFileName(imgRow, prefixDir, outputDir, copyFiles=False):
    index, imgPath = relionToLocation(imgRow.get(md.RLN_IMAGE_NAME))
    baseName = os.path.basename(imgPath)
    newName = os.path.join(outputDir, baseName)
    if not os.path.exists(newName):
        if copyFiles:
            pwutils.copyFile(os.path.join(prefixDir, imgPath), newName)
        else:
            pwutils.createLink(os.path.join(prefixDir, imgPath), newName)

    imgRow.set(md.RLN_IMAGE_NAME, locationToRelion(index, newName))


def setupCTF(imgRow, sampling):
    """ Do some validations and set some values
    for Relion import.
    """
    imgRow.set(md.MDL_SAMPLINGRATE, sampling)
    # TODO: check if we want to move this behaviour to setup CTFModel by default
    hasDefocusU = imgRow.containsLabel(md.MDL_CTF_DEFOCUSU)
    hasDefocusV = imgRow.containsLabel(md.MDL_CTF_DEFOCUSV)
    hasDefocusAngle = imgRow.containsLabel(md.MDL_CTF_DEFOCUS_ANGLE)

    if hasDefocusU or hasDefocusV:
        if not hasDefocusU:
            imgRow.set(md.MDL_CTF_DEFOCUSU,
                       imgRow.get(md.MDL_CTF_DEFOCUSV))
        if not hasDefocusV:
            imgRow.set(md.MDL_CTF_DEFOCUSV,
                       imgRow.get(md.MDL_CTF_DEFOCUSU))
        if not hasDefocusAngle:
            imgRow.set(md.MDL_CTF_DEFOCUS_ANGLE, 0.)
