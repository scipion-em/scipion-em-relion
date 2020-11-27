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
import relion

from .protocol_create_mask3d import ProtRelionCreateMask3D
from .protocol_classify2d import ProtRelionClassify2D
from .protocol_classify3d import ProtRelionClassify3D
from .protocol_refine3d import ProtRelionRefine3D
from .protocol_reconstruct import ProtRelionReconstruct
from .protocol_postprocess import ProtRelionPostprocess
from .protocol_sort import ProtRelionSortParticles
from .protocol_expand_symmetry import ProtRelionExpandSymmetry
from .protocol_initialmodel import ProtRelionInitialModel
from .protocol_localres import ProtRelionLocalRes
from .protocol_export_coords import ProtRelionExportCoordinates
from .protocol_extract_particles import ProtRelionExtractParticles
from .protocol_export_ctf import ProtRelionExportCtf
from .protocol_center_averages import ProtRelionCenterAverages
from .protocol_export_particles import ProtRelionExportParticles
from .protocol_autopick_ref import ProtRelion2Autopick
from .protocol_autopick_log import ProtRelionAutopickLoG
from .protocol_remove_views import ProtRelionRemovePrefViews
from .protocol_crop_resize import ProtRelionResizeVolume

from .protocol_motioncor import ProtRelionMotioncor
from .protocol_multibody import ProtRelionMultiBody
from .protocol_symmetrize_volume import ProtRelionSymmetrizeVolume

if relion.Plugin.IS_30():
    from ._legacy.protocol30_ctf_refinement import ProtRelionCtfRefinement
    from ._legacy.protocol30_preprocess import ProtRelionPreprocessParticles
    from ._legacy.protocol30_subtract import ProtRelionSubtract
    from ._legacy.protocol30_bayesian_polishing import ProtRelionBayesianPolishing
else:
    from .protocol_ctf_refinement import ProtRelionCtfRefinement
    from .protocol_preprocess import ProtRelionPreprocessParticles
    from .protocol_subtract import ProtRelionSubtract
    from .protocol_bayesian_polishing import ProtRelionBayesianPolishing

# New protocol from Relion v3.1:
# FIXME: this has to be imported for tests / wizards to work
from .protocol_assign_optic_groups import ProtRelionAssignOpticsGroup
from .protocol_compress_movies import ProtRelionCompressMovies
from .protocol_compress_estimate_gain import ProtRelionCompressEstimateGain
