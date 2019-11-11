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
New conversion functions dealing with Relion3.1 new star files format.
"""

import os
from os.path import join, basename
import numpy
from itertools import izip
from collections import OrderedDict

import pyworkflow as pw
from pyworkflow.object import ObjectWrap, String, Integer
import pyworkflow.em.metadata as md
from pyworkflow.em.constants import NO_INDEX

from relion.constants import *
from .metadata import Table, Column


class SetOfImagesWriter:
    """ Helper class to convert from Scipion SetOfImages subclasses
    into Relion>3.1 star files (and binaries if conversion needed).
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
        pass

    def writeSetOfMovies(self, movieSet, starFile):
        pass  # This is used later in the pipeline in extract particles movies

    def writeSetOfMicrographs(self, micsIterable, starFile):
        # Process the first item and create the table based
        # on the generated columns
        self._optics = OrderedDict()
        micRow = OrderedDict()
        iterMics = iter(micsIterable)
        mic = next(iterMics)
        self._micToRow(mic, micRow)

        opticsTable = self._createTableFromDict(self._optics.values()[0])
        micsTable = self._createTableFromDict(micRow)

        while mic is not None:
            self._micToRow(mic, micRow)
            micsTable.addRow(**micRow)
            mic = next(iterMics, None)

        for opticsDict in self._optics.values():
            opticsTable.addRow(**opticsDict)

        with open(starFile, 'w') as f:
            f.write("# Star file generated with Scipion\n")
            f.write("# version 30001\n")
            opticsTable.writeStar(f, tableName='optics')
            f.write("# version 30001\n")
            micsTable.writeStar(f, tableName='micrographs')

    def _createTableFromDict(self, rowDict):
        """ Helper function to create a Table instance from
        an input dict with keys as columns names and type
        the type of the values in the dict.
        """
        from pyworkflow.utils import prettyDict
        prettyDict(rowDict)

        return Table(columns=[
            Column(k, type=type(v)) for k, v in rowDict.iteritems()])

    def _micToRow(self, mic, row):
        row['rlnMicrographName'] = mic.getFileName()

        acq = mic.getAcquisition()
        ogName = acq.opticsGroupName.get()

        if not ogName in self._optics:
            ogNumber = len(self._optics) + 1
            self._optics[ogName] = {
                'rlnOpticsGroupName': ogName,
                'rlnOpticsGroup': ogNumber,
                'rlnMtfFileName': acq.mtfFile.get(),
                # FIXME: Check when we need to update the following
                'rlnMicrographOriginalPixelSize': mic.getSamplingRate(),
                'rlnVoltage': acq.getVoltage(),
                'rlnSphericalAberration': acq.getSphericalAberration(),
                'rlnAmplitudeContrast': acq.getAmplitudeContrast(),
                'rlnMicrographPixelSize': mic.getSamplingRate()
            }
        else:
            ogNumber = self._optics[ogName]['rlnOpticsGroup']

        row['rlnOpticsGroup'] = ogNumber

        if mic.hasCTF():
            self._ctfToRow(mic.getCTF(), row)

    def _ctfToRow(self, ctf, row):
        row['rlnCtfImage'] = ctf.getPsdFile()
        dU, dV, dAngle = ctf.getDefocus()
        row['rlnDefocusU'] = dU
        row['rlnDefocusV'] = dV
        # FIXME Check how astigmatism is defined in Relion
        row['rlnCtfAstigmatism'] = dU / dV
        row['rlnDefocusAngle'] = dAngle
        row['rlnCtfFigureOfMerit'] = ctf.getFitQuality()
        row['rlnCtfMaxResolution'] = ctf.getResolution()

        phaseShift = ctf.getPhaseShift()

        if phaseShift is not None:
            row['rlnCtfPhaseShift'] = phaseShift

"""

# version 30001

data_optics

loop_
_rlnOpticsGroupName #1
_rlnOpticsGroup #2
_rlnMtfFileName #3
_rlnMicrographOriginalPixelSize #4
_rlnVoltage #5
_rlnSphericalAberration #6
_rlnAmplitudeContrast #7
_rlnMicrographPixelSize #8
opticsGroup1            1 mtf_k2_200kV.star     0.885000   200.000000     1.400000     0.100000     0.885000


# version 30001

data_micrographs

loop_
_rlnMicrographName #1
_rlnOpticsGroup #2
_rlnCtfImage #3
_rlnDefocusU #4
_rlnDefocusV #5
_rlnCtfAstigmatism #6
_rlnDefocusAngle #7
_rlnCtfFigureOfMerit #8
_rlnCtfMaxResolution #9
MotionCorr/job002/Movies/20170629_00021_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00021_frameImage_PS.ctf:mrc 10863.857422 10575.721680   288.135742    77.967194     0.131144     4.809192
MotionCorr/job002/Movies/20170629_00022_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00022_frameImage_PS.ctf:mrc  9836.475586  9586.718750   249.756836    70.291290     0.168926     3.619038
MotionCorr/job002/Movies/20170629_00023_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00023_frameImage_PS.ctf:mrc 10678.056641 10365.234375   312.822266    71.958046     0.166332     3.412236
MotionCorr/job002/Movies/20170629_00024_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00024_frameImage_PS.ctf:mrc 11693.039062 11353.885742   339.153320    74.367035     0.153686     3.495461
MotionCorr/job002/Movies/20170629_00025_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00025_frameImage_PS.ctf:mrc 10656.021484 10381.144531   274.876953    73.184746     0.156609     3.478493
MotionCorr/job002/Movies/20170629_00026_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00026_frameImage_PS.ctf:mrc  8346.364258  8101.361816   245.002441    72.656509     0.187263     2.998199
MotionCorr/job002/Movies/20170629_00027_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00027_frameImage_PS.ctf:mrc  9145.218750  8826.969727   318.249023    78.355713     0.173906     3.213317
MotionCorr/job002/Movies/20170629_00028_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00028_frameImage_PS.ctf:mrc 10159.675781  9805.832031   353.843750    75.627907     0.164042     3.075406
MotionCorr/job002/Movies/20170629_00029_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00029_frameImage_PS.ctf:mrc  8303.106445  7963.178711   339.927734    79.366615     0.193207     3.272007
MotionCorr/job002/Movies/20170629_00030_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00030_frameImage_PS.ctf:mrc  9116.507812  8828.420898   288.086914    76.925766     0.202328     3.348456
MotionCorr/job002/Movies/20170629_00031_frameImage.mrc            1 CtfFind/job003/Movies/20170629_00031_frameImage_PS.ctf:mrc 10123.724609  9808.489258   315.235352    75.411964     0.198348     3.198971
"""