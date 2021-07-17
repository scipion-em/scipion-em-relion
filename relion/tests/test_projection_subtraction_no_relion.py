# ***************************************************************************
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
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
# ***************************************************************************/

import os
from tempfile import NamedTemporaryFile

from pwem.protocols import ProtImportParticles
from pyworkflow.tests import BaseTest, setupTestProject

from ..protocols import ProtRelionSubtract


class TestSubtractionProjection(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.xmippAvailable = True
        cls.sampling = 1.30
        from pwem import Domain
        try:
            cls.xmipp3 = \
                Domain.importFromPlugin('xmipp3', doRaise=True)
        except:
            cls.xmippAvailable = False

    def _create3Dmap(self, fileName, label="feat mask"):
        from xmipp3.protocols import XmippProtCreateMask3D
        protMask = self.newProtocol(XmippProtCreateMask3D,
                                    featureFilePath=fileName,
                                    source=2,
                                    samplingRate=self.sampling)
        protMask.setObjLabel(label)
        self.launchProtocol(protMask)
        return protMask.outputMask, protMask

    def createVolume(self):
        # volume
        f = NamedTemporaryFile(delete=False, suffix=".feat")
        command = """# XMIPP_STAR_1 *
# Type of feature (sph, blo, gau, Cyl, dcy, cub, ell, con)(Required)	
data_block1
 _dimensions3D  '752 752 752'
 _phantomBGDensity  0.
 _scale  0.1702127659574468
data_block2
loop_
 _featureType
 _featureOperation
 _featureDensity
 _featureCenter
 _featureSpecificVector	
blo	+	1	'-165.51	0.00	267.78'	'60 1.4 2'
blo	+	1	'165.51	0.00	267.78'		'60 1.4 2'
blo	+	1	'0.00	267.78	165.51'		'60 1.4 2'
blo	+	1	'-267.78	165.51	0.00'		'60 1.4 2'
blo	+	1	'-267.78	-165.51	0.00'		'60 1.4 2'
blo	+	1	'0.00	-267.78	165.51'		'60 1.4 2'
blo	+	1	'165.51	0.00	-267.78'		'60 1.4 2'
blo	+	1	'-165.51	0.00	-267.78'		'60 1.4 2'
blo	+	1	'0.00	267.78	-165.51'		'60 1.4 2'
blo	+	1	'267.78	165.51	0.00'		'60 1.4 2'
blo	+	1	'267.78	-165.51	0.00'		'60 1.4 2'
blo	+	1	'0.00	-267.78	-165.51'		'60 1.4 2'
"""
        f.write(command.encode('utf8'))
        f.close()
        return self._create3Dmap(f.name, "create volume")

    def createMask(self):
        # volume
        f = NamedTemporaryFile(delete=False, suffix=".feat")
        command = """# XMIPP_STAR_1 *
# Type of feature (sph, blo, gau, Cyl, dcy, cub, ell, con)(Required)	
data_block1
 _dimensions3D  '752 752 752'
 _phantomBGDensity  0.
 _scale  0.1702127659574468
data_block2
loop_
 _featureType
 _featureOperation
 _featureDensity
 _featureCenter
 _featureSpecificVector
sph	+	1.	'-165.51 0.00 267.78'	'65'
"""
        f.write(command.encode('utf8'))
        f.close()
        return self._create3Dmap(f.name, "create mask")

    def runImportParticles(self, pattern, samplingRate):
        """ Run an Import particles protocol. """
        protImport = self.newProtocol(ProtImportParticles,
                                      objLabel='from Xmipp createprojectionlibrary',
                                      importFrom=ProtImportParticles.IMPORT_FROM_XMIPP3,
                                      mdFile=pattern,
                                      magnification=65000,
                                      samplingRate=samplingRate,
                                      haveDataBeenPhaseFlipped=True
                                      )
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputParticles,
                             "There was a problem with import particles.")
        return protImport

    def testProjectionSubtraction(self):
        if not self.xmippAvailable:
            return
        # create volume
        volume, protocolVolume = self.createVolume()
        # create mask
        mask, protocolMask = self.createMask()
        # create projection
        # unfortunately there is no available protocol
        # so we will run xmipp_angular_project_library
        # -i kk.mrc  -o kk.mrcs  --sampling_rate 15 --method real_space
        docFileParticles = os.path.abspath(protocolVolume._getExtraPath("proj"))
        args = " -i %s" % volume.getFileName()
        args += " -o %s" % docFileParticles + ".mrcs"
        args += " --sampling_rate 15"
        args += " --method real_space"
        progname = "xmipp_angular_project_library"
        self.xmipp3.Plugin.runXmippProgram(progname, args)

        # rename doc file to xmp since import assumes the
        # xmipp metadata files extension is xmd.
        os.rename(docFileParticles + ".doc", docFileParticles + ".xmd")
        # import projection
        protImport = self.runImportParticles(docFileParticles + ".xmd",
                                             1.3)
        # subtract
        protSubtract = self.newProtocol(ProtRelionSubtract,
                                        relionInput=False,
                                        refMask=mask,
                                        doCTF=False,
                                        numberOfMpi=1)
        protSubtract.inputParticlesAll.set(protImport.outputParticles)
        protSubtract.inputVolume.set(volume)
        self.launchProtocol(protSubtract)
        self.assertIsNotNone(protSubtract.outputParticles,
                             "There was a problem with subtract projection")
