# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)     [2]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
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


import os.path
from pathlib import Path
from typing import Optional
import torch
import numpy as np
import mrcfile

from torch.utils.data import DataLoader

from dynamight.data.handlers.particle_image_preprocessor import ParticleImagePreprocessor
from dynamight.data.dataloaders.relion import RelionDataset
from dynamight.evaluation.utils import compute_latent_space_and_colors


class HeterogeneityProgramInterface:
    def __init__(self, _path_template: str, _program_loading_params: dict):
        self.model = self.prepare_heterogeneity_program(**_program_loading_params)
        self.path_template = _path_template

    def prepare_heterogeneity_program(self, **kwargs) -> object:
        gpu_id = kwargs.pop("gpu_id", None)
        checkpoint_file = kwargs.pop("checkpoint_file", None)

        self.device = "cpu" if gpu_id is None else 'cuda:' + str(int(gpu_id))
        if checkpoint_file is None:
            checkpoint_file = output_directory / \
                              'forward_deformations/checkpoints/checkpoint_final.pth'
        cp = torch.load(checkpoint_file, map_location=self.device , weights_only=False)

        refinement_star_file = cp['refinement_directory']
        if refinement_star_file.suffix == '.star':
            pass
        else:
            refinement_star_file = refinement_star_file / 'run_data.star'

        encoder_h1 = cp['encoder_half1']
        decoder_h1 = cp['decoder_half1']
        encoder_h2 = cp['encoder_half2']
        decoder_h2 = cp['decoder_half2']

        poses = cp['poses']

        relion_dataset = RelionDataset(
            path=refinement_star_file.resolve(),
            circular_mask_thickness=20,
            particle_diameter=None,
        )
        dataset = relion_dataset.make_particle_dataset()

        encoder_h1.load_state_dict(
            cp['encoder_half1_state_dict'])
        decoder_h1.load_state_dict(
            cp['decoder_half1_state_dict'])
        encoder_h2.load_state_dict(
            cp['encoder_half2_state_dict'])
        decoder_h2.load_state_dict(
            cp['decoder_half2_state_dict'])

        poses.load_state_dict(cp['poses_state_dict'])

        '''Computing indices for the second half set'''
        self.indices_h1 = cp['indices_half1'].cpu().numpy()
        self.indices_h2 = np.asarray(list(set(range(len(dataset))) - set(list(self.indices_h1))))

        decoder_h1.p2i.device = self.device
        decoder_h1.projector.device = self.device
        decoder_h1.image_smoother.device = self.device
        decoder_h1.p2v.device = self.device
        decoder_h1.device = self.device
        decoder_h1.to(self.device )

        decoder_h2.p2i.device = self.device
        decoder_h2.projector.device = self.device
        decoder_h2.image_smoother.device = self.device
        decoder_h2.p2v.device = self.device
        decoder_h2.device = self.device
        decoder_h2.to(self.device )

        return [decoder_h1, decoder_h2]

    def decode_state_from_latent(self, latent: np.array) -> None:
        latent = torch.from_numpy(latent.astype(np.float32)).to(self.device)
        r = torch.zeros([2, 3]).to(self.device)
        t = torch.zeros([2, 2]).to(self.device)
        idx = 1
        for l in latent:
            l = l[None, ...]
            vol_h1 = self.model[0].generate_volume(l, r, t).detach().float().cpu().numpy()
            vol_h2 = self.model[1].generate_volume(l, r, t).detach().float().cpu().numpy()
            vol = 0.5 * (vol_h1 + vol_h2)
            with mrcfile.new(self.path_template.format(idx)) as mrc:
                mrc.set_data(vol)
            idx += 1
