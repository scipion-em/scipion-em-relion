# ***************************************************************************
# * Authors:   Roberto Marabini (roberto@cnb.csic.es)
# *
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

from os.path import exists

from pyworkflow.tests import BaseTest, setupTestProject
import pwem.protocols as emprot
from ..protocols import ProtRelionModelAngelo


class TestModelAngelo(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_ProtModelAngeloBuild(self):
        # import vol
        EMDBiD = 406
        args = {
                'importFrom': emprot.ProtImportVolumes.IMPORT_FROM_EMDB,
                'emdbId': EMDBiD,
                }

        prot = self.newProtocol(emprot.ProtImportVolumes, **args)
        prot.setObjLabel('import vol')
        self.launchProtocol(prot)
        vol = prot.outputVolume

        # import sequence
        PDBiD = '6nbb'
        args = {'inputSequenceName': PDBiD,
                'inputProteinSequence': emprot.ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'pdbId': PDBiD,
                'inputStructureChain': '{"model": 0, "chain": "A", "residues": 119}',
                }

        prot1 = self.newProtocol(emprot.ProtImportSequence, **args)
        prot1.setObjLabel('import aminoacid seq,\n from uniprot id')
        self.launchProtocol(prot1)
        sequence = prot1.outputSequence
        # run protocol
        listSeq = [sequence, sequence]
        args = {'inputVolume': vol,
                'inputSequenceS': listSeq,
                'useGpu': True,
                'gpuList': "0",
                }
        prot2 = self.newProtocol(ProtRelionModelAngelo, **args)
        prot2.setObjLabel('model angelo')
        self.launchProtocol(prot2)
        # check results
        self.assertTrue(exists(prot2.pruned.getFileName()))
        self.assertTrue(exists(prot2.raw.getFileName()))

    def test_ProtModelAngeloBuildNoSeq(self):
        # import vol
        EMDBiD = 406
        args = {
                'importFrom': emprot.ProtImportVolumes.IMPORT_FROM_EMDB,
                'emdbId': EMDBiD,
                }

        prot = self.newProtocol(emprot.ProtImportVolumes, **args)
        prot.setObjLabel('import vol')
        self.launchProtocol(prot)
        vol = prot.outputVolume

        # run protocol
        args = {'inputVolume': vol,
                'useGpu': True,
                'gpuList': "0",
                }
        prot2 = self.newProtocol(ProtRelionModelAngelo, **args)
        prot2.setObjLabel('model angelo\n No sequence')
        self.launchProtocol(prot2)
        # check results
        self.assertTrue(exists(prot2.pruned.getFileName()))
