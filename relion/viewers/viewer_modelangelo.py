# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es)
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
# *
# **************************************************************************

import os
from pyworkflow.viewer import Viewer
from pwem.viewers.viewer_chimera import Chimera
from pwem.emlib.image import ImageHandler
from ..protocols.protocol_modelangelo import ProtRelionModelAngelo


class ProtModelAngeloViewer(Viewer):
    _label = 'viewer model angelo'
    _targets = [ProtRelionModelAngelo]

    def visInputVolume(self, f, vol, counter):
        inputVolFileName = ImageHandler.removeFileType(vol.getFileName())
        f.write("open %s\n" % inputVolFileName)
        # model_angelo does not honor origin in 3D maps
        # if vol.hasOrigin():
        #    x, y, z = vol.getOrigin().getShifts()
        # else:
        #    x, y, z = vol.getOrigin(force=True).getShifts()
        f.write("volume #%d style surface voxelSize %f\n"
                "volume #%d origin 0,0,0\n"
                % (counter, vol.getSamplingRate(), counter))

    def visualize(self, obj, **args):
        # create axis file
        models = 1
        _inputVol = self.protocol.inputVolume.get()
        if _inputVol is not None:
            dim = _inputVol.getDim()[0]
            sampling = _inputVol.getSamplingRate()

        bildFileName = os.path.abspath(self.protocol._getExtraPath("axis_input.bild"))
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=bildFileName,
                                         sampling=sampling)

        fnCmd = self.protocol._getExtraPath("model_angelo_viewer.cxc")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % bildFileName)
        models += 1
        f.write("cofr 0,0,0\n")  # set center of coordinates
        # change to workingDir
        # If we do not use cd and the project name has a space
        # the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())

        # get path to atomstructs
        for output in self.protocol._outputs:
            # if the file is an atomic struct show it in chimera
            fileName = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
            if fileName.endswith(".cif") or fileName.endswith(".pdb"):
                f.write("open %s\n" % fileName)
                models += 1
        # set alphafold colormap
        f.write("color bfactor palette alphafold\n")
        f.write("key red:low orange: yellow: cornflowerblue: blue:high\n")
        # open 3D map
        models += 1
        self.visInputVolume(f, _inputVol, models)
        f.write("color #%d #bf404071 models\n" % models)
        f.write("show cartoons\n")
        f.write("hide atoms\n")
        f.close()
        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")
