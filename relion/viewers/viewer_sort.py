# ******************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es) [1]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pyworkflow.viewer import Viewer

from .viewer_base import *
from ..protocols import ProtRelionSortParticles


# TODO: deprecate this class
class RelionSortViewer(Viewer):
    """ Visualization of Relion sorting results."""

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtRelionSortParticles]

    def _visualize(self, obj, **kwargs):
        views = []

        if obj.hasAttribute('outputParticles'):  # Protocol finished
            particles = obj.outputParticles
            labels = ('id enabled _index _filename _rlnSelectParticlesZscore '
                      '_coordinate._rlnAutopickFigureOfMerit _sampling '
                      '_ctfModel._defocusU _ctfModel._defocusV '
                      '_ctfModel._defocusAngle _transform._matrix')
            sortBy = '_rlnSelectParticlesZscore asc'
            strId = particles.strId()
            fn = particles.getFileName()
            views.append(ObjectView(self._project, strId, fn,
                                    viewParams={showj.ORDER: labels,
                                                showj.VISIBLE: labels,
                                                showj.SORT_BY: sortBy,
                                                showj.RENDER: '_filename'}))

            fn = obj._getExtraPath('input_particles_sorted.star')
            mdFn = md.MetaData(fn)
            # If Zscore in output images plot Zscore particle sorting
            if mdFn.containsLabel(md.RLN_SELECT_PARTICLES_ZSCORE):
                # sort output by Z-score
                mdFn.sort(md.RLN_SELECT_PARTICLES_ZSCORE)
                xplotter = RelionPlotter(windowTitle="Zscore particles sorting")
                xplotter.createSubPlot("Particle sorting", "Particle number",
                                       "Zscore")
                xplotter.plotMd(mdFn, False,
                                mdLabelY=md.RLN_SELECT_PARTICLES_ZSCORE)
                views.append(xplotter)

        return views
