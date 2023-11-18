# ******************************************************************************
# *
# * Authors:     Grigory Sharov     (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology, MRC-LMB
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
# ******************************************************************************
import logging
logger = logging.getLogger(__name__)

from pyworkflow.protocol.constants import STATUS_FINISHED

from .viewer_base import *
from ..protocols import ProtRelionMultiBody


class MultibodyViewer(RelionViewer):
    """ Visualization of Relion 3D multi-body results. """
    _targets = [ProtRelionMultiBody]
    _environments = [DESKTOP_TKINTER]

    _label = 'viewer multi-body'

    def _defineParams(self, form):
        showAnimation = self.protocol.runFlexAnalysis.get()

        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam,
                      choices=['last', 'selection'], default=ITER_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize",
                      help="""
        *last*: only the last iteration will be visualized.
        *selection*: you may specify a range of iterations.
        Examples:
        "1,5-8,10" -> [1,5,6,7,8,10]
        "2,6,9-11" -> [2,6,9,10,11]
        "2 5, 6-8" -> [2,5,6,7,8]                      
                                   """)
        form.addParam('iterSelection', params.NumericRangeParam,
                      condition='viewIter==%d' % ITER_SELECTION,
                      label="Iterations list",
                      help="Write the iteration list to visualize.")

        group = form.addGroup('3D Bodies')
        group.addParam('showHalves', params.EnumParam, default=3,
                       choices=['half1', 'half2', 'both', 'final'],
                       label='Volume to visualize',
                       help='Select which half do you want to visualize.')
        group.addParam('displayVol', params.EnumParam,
                       choices=['slices', 'chimera'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       default=VOLUME_SLICES,
                       label='Display volumes with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')

        if showAnimation:
            group = form.addGroup('Flexibility analysis')
            group.addParam('component', params.IntParam, default=1,
                           label="Principal component to show")
            group.addParam('showMovie', params.LabelParam,
                           label="Show Chimera animation")

    def _getVisualizeDict(self):
        self._load()
        return {'displayVol': self._showVolumes,
                'showMovie': self._showMovie
                }

    def _createVolumesSqlite(self):
        """ Write a sqlite with all volumes selected for visualization. """
        path = self.protocol._getExtraPath('relion_viewer_volumes.sqlite')
        if self.protocol.doContinue:
            protRef = self.protocol.continueRun.get().protRefine.get()
        else:
            protRef = self.protocol.protRefine.get()

        samplingRate = protRef._getInputParticles().getSamplingRate()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            if not os.path.exists(volFn):
                raise FileNotFoundError(f"Missing volume file: {volFn}\n Please select "
                                        "a valid class or iteration number.")
            logger.debug(f"Adding vol: {volFn}")
            if not volFn.endswith(":mrc"):
                files.append(volFn + ":mrc")

        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

    def _getVolumeNames(self):
        vols = []
        prefixes = self._getVolumePrefixes()
        num = self.protocol._getNumberOfBodies()
        logger.debug(f"self._iterations: {self._iterations}")
        for it in self._iterations:
            for ref3d in range(1, num+1):
                for prefix in prefixes:
                    volFn = self.protocol._getFileName(prefix + 'volume_mbody',
                                                       iter=it, ref3d=ref3d)
                    if os.path.exists(volFn):
                        volFn = volFn.replace(":mrc", "")
                        vols.append(volFn)
                    else:
                        raise FileNotFoundError(f"Volume {volFn} does not exists.\n"
                                                "Please select a valid iteration number.")
        return vols

    @protected_show
    def _showMovie(self, paramName=None):
        """ Create a chimera script to visualize animation. """
        prot = self.protocol
        if prot.getStatus() != STATUS_FINISHED:
            raise ValueError("Protocol has not finished yet, results are not ready!")

        compToShow = self.component.get()
        totalComp = prot.numberOfEigenvectors.get()
        if compToShow > totalComp:
            raise ValueError("You only have %d components!" % totalComp)

        vols = "analyse_component%03d_bin*.mrc" % compToShow
        cmdFile = prot._getExtraPath('chimera_animation_comp%d.cxc' % compToShow)

        with open(cmdFile, 'w') as f:
            f.write("open %s vseries true\n" % vols)
            f.write('vseries play #1 loop true direction oscillate\n')

        views = ChimeraView(cmdFile)

        return [views]
