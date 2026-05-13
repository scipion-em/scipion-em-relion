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

import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
from pwem.objects import Volume
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtAlignVolume


class ProtRelionSymmetrizeVolume(ProtAlignVolume):
    """
    Symmetrizes a three-dimensional cryo-EM volume according to a user-defined
    point-group symmetry in order to improve structural consistency and enhance
    interpretable signal across equivalent regions of the map.

    AI Generated:

    Symmetrize Volume (ProtRelionSymmetrizeVolume) — User Manual
        Overview

        The Symmetrize Volume protocol applies symmetry operations to a 3D
        cryo-EM reconstruction using Relion tools. Its primary objective is to
        transform a reconstructed density map into a symmetrized representation
        consistent with the biological symmetry of the macromolecular assembly.
        This process is commonly used for particles that naturally adopt cyclic,
        dihedral, tetrahedral, octahedral, or icosahedral symmetry.

        In cryo-EM workflows, symmetry application is often essential for
        improving map quality and interpretability. By averaging equivalent
        structural regions, symmetry enhances the signal-to-noise ratio and can
        reveal features that are difficult to observe in asymmetric maps.
        Symmetrization is especially important in studies of viruses, molecular
        cages, membrane channels, and highly ordered oligomeric assemblies.

        Biological Meaning of Symmetry

        Biological symmetry reflects the repeated organization of identical or
        near-identical subunits within a macromolecular complex. Correctly
        applying symmetry allows the reconstruction to better represent the
        expected structural arrangement and often improves the apparent
        resolution of the final map.

        However, symmetry should only be imposed when biologically justified.
        Applying incorrect symmetry can introduce artificial features, obscure
        meaningful conformational variability, and distort biologically relevant
        asymmetry. For this reason, users should carefully evaluate the
        structural evidence supporting the selected symmetry group before
        running the protocol.

        Input Volume Requirements

        The protocol requires a single input volume representing the
        reconstruction to be symmetrized. The volume should already be properly
        reconstructed and approximately centered. Maps that are severely
        misaligned or contain strong artifacts may produce unstable or
        biologically misleading results after symmetry application.

        The voxel size associated with the input map is preserved throughout the
        process. Consistent sampling information is important because the
        protocol uses the physical scale of the reconstruction during alignment
        and symmetry operations.

        Symmetry Selection

        The user specifies the symmetry group to apply to the reconstruction.
        Common examples include c1 for asymmetric maps, c2 or higher cyclic
        symmetries for ring-like complexes, dn symmetries for dihedral
        assemblies, and higher-order point groups for highly symmetric
        particles.

        Choosing the correct symmetry is one of the most critical decisions in
        structural analysis. Even when a particle is expected to be symmetric,
        partial flexibility or compositional heterogeneity may reduce the
        validity of strict symmetry averaging. In such situations, users often
        compare asymmetric and symmetrized reconstructions to determine whether
        important biological information is being lost.

        Alignment Before Symmetrization

        Before symmetry is imposed, the protocol aligns the input map relative
        to the selected symmetry axes. This step ensures that the symmetry
        operators are applied consistently and that equivalent regions overlap
        correctly in three-dimensional space.

        Proper alignment is especially important for rotational symmetries,
        where even small angular deviations can reduce the quality of the final
        symmetrized map. For highly symmetric particles, accurate orientation
        relative to the symmetry axis is essential for achieving optimal signal
        enhancement.

        Outputs and Interpretation

        The protocol produces an aligned volume and a symmetrized volume. The
        aligned map represents the reconstruction positioned according to the
        selected symmetry convention, while the symmetrized map contains the
        final averaged density after symmetry application.

        From a biological perspective, the symmetrized reconstruction often
        appears cleaner and more interpretable than the original map because
        repeated structural information is reinforced. This can facilitate model
        building, map interpretation, visualization, and downstream refinement.

        Nevertheless, users should remain cautious when interpreting subtle
        structural features. Flexible domains, asymmetric ligands, or partial
        occupancies may become blurred or disappear entirely after symmetry
        averaging. Comparing the symmetrized and unsymmetrized maps is often a
        valuable strategy for identifying biologically meaningful deviations
        from perfect symmetry.

        Practical Recommendations

        In routine cryo-EM analysis, it is generally advisable to first inspect
        the asymmetric reconstruction before imposing symmetry. If the assembly
        clearly exhibits repeated organization and no major asymmetric features
        are present, symmetrization can substantially improve map quality.

        For oligomeric complexes with uncertain symmetry, exploratory runs using
        different symmetry groups may help determine the most biologically
        appropriate choice. However, users should avoid selecting higher-order
        symmetries solely to improve apparent resolution, since artificially
        imposed symmetry can generate misleading structural interpretations.

        When preparing maps for publication or atomic modeling, it is often good
        practice to document whether symmetry was applied and to preserve access
        to the original asymmetric reconstruction for comparison and validation.

        Final Perspective

        Symmetry application is one of the most powerful tools in single-particle
        cryo-EM because it can dramatically improve map clarity and structural
        detail when used appropriately. At the same time, it requires careful
        biological judgment. Reliable results depend on selecting the correct
        symmetry group, validating the presence of genuine structural symmetry,
        and understanding how averaging may influence the interpretation of
        flexible or asymmetric regions.
    """
    _label = 'symmetrize volume'
    _devStatus = PROD
    _possibleOutputs = {
        'outputVolumeAligned': Volume,
        'outputVolumeSymmetrized': Volume
    }
    
    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume',
                      label="Input volume", important=True,
                      help='Select the input volume to be symmetrized. ')

        form.addParam('symmetryGroup', params.StringParam, default='c1',
                      label="Symmetry",
                      help='Select which symmetry do you want to apply. ')
        
        form.addParallelSection(threads=0, mpi=0)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep,
                                 self.inputVolume.getObjId(),
                                 needsGPU=False)

    # --------------------------- STEPS functions ------------------------------
    def createOutputStep(self, volId):
        sym = self.symmetryGroup.get()

        inFn = self._getPath('input_volume.mrc')
        alignedFn = self._getPath('volume_aligned_sym%s.mrc' % sym)
        symFn = self._getPath('volume_sym%s.mrc' % sym)
        pixSize = self.inputVolume.get().getSamplingRate()

        ImageHandler().convert(self.inputVolume.get(), inFn)

        self.runJob("relion_align_symmetry",
                    "--i %s --o %s --sym %s --angpix %0.5f" % (
                        inFn, alignedFn, sym, pixSize))

        self.runJob("relion_image_handler",
                    "--i %s --o %s --sym %s" % (alignedFn, symFn, sym))

        def _defineOutputVol(name, fn):
            vol = Volume()
            vol.copyInfo(self.inputVolume.get())
            vol.setLocation(fn)
            self._defineOutputs(**{name: vol})
            self._defineTransformRelation(self.inputVolume, vol)

        _defineOutputVol('outputVolumeAligned', alignedFn)
        _defineOutputVol('outputVolumeSymmetrized', symFn)

    # -------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumeSymmetrized'):
            summary.append("Output is not ready yet.")
        else:
            summary.append("Symmetry used: *%s*" % self.symmetryGroup.get())
        return summary
