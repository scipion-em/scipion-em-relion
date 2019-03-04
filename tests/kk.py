from __future__ import print_function
"""
edit file postprocess.star and add to

data_general

_rlnFinalResolution                           5.541053
_rlnBfactorUsedForSharpening               -439.707536
_rlnFittedSlopeGuinierPlot                 -109.926884
_rlnFittedInterceptGuinierPlot              -13.944270
_rlnCorrelationFitGuinierPlot                 0.962517


_rlnUnfilteredMapHalf1                    Runs/008287_ProtRelionRefine3D/extra/relion_half1_class001_unfil.mrc
_rlnUnfilteredMapHalf2                    Runs/008287_ProtRelionRefine3D/extra/relion_half1_class002_unfil.mrc
_rlnMaskName                              Runs/009075_ProtRelionPostprocess/extra/cess/input/postprocess_automask.mrc
_rlnRandomiseFrom                             7.91579

Randon noise is in the stdout of postprocess

"""

"""
file input_particles.start
does not have column rlnRandomSubset
put this available in star file in
"""

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.em.constants import ALIGN_PROJ

import relion
from relion.convert.metadata import Table

# fnStar = '/home/roberto/ScipionUserData/projects/20170529_CARMENSM_delta7_crio_tomo_4_control_wt/Runs/008287_ProtRelionRefine3D/extra/relion_data.star'
#fnStar = '/home/roberto/tmp/relion_data.star'
fnStar = '/home/roberto/ScipionUserData/projects/20170529_CARMENSM_delta7_crio_tomo_4_control_wt/Runs/010387_ProtRelionCtfRefinement/extra/particles_ctf_refine.star'
mdAll = md.MetaData(fnStar)
# extraLabels = [md.RLN_PARTICLE_RANDOM_SUBSET]
extraLabels = [md.RLN_PARTICLE_RANDOM_SUBSET,
               md.RLN_IMAGE_BEAMTILT_X,
               md.RLN_IMAGE_BEAMTILT_Y]
imgSet = em.SetOfParticles(filename="particles_v2.sqlite")
alignType = ALIGN_PROJ
relion.convert.readSetOfParticles(fnStar,imgSet,
                                  alignType=alignType,
                                  extraLabels=extraLabels)  # _rlnRandomSubset
imgSet.write()

# writeSqliteIterData
#relion.convert.writeSetOfParticles(inputParts, imgStar, inputFolder,
#                                    alignType=em.ALIGN_PROJ,
#                                    fillMagnification=True,
#                                    fillRandomSubset=True)


