#!/usr/bin/env python
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

from torch.utils.data import DataLoader

from dynamight.data.handlers.particle_image_preprocessor import ParticleImagePreprocessor
from dynamight.data.dataloaders.relion import RelionDataset
from dynamight.evaluation.utils import compute_latent_space_and_colors


def encode_latent_space(
    output_directory: Path,
    checkpoint_file: Optional[Path] = None,
    batch_size: int = 100,
    gpu_id: Optional[int] = 0,
    n_workers: int = 8,
    reduce_by_deformation: bool = False,
):
    device = "cpu" if gpu_id is None else 'cuda:' + str(gpu_id)
    if checkpoint_file is None:
        checkpoint_file = output_directory / \
            'forward_deformations/checkpoints/checkpoint_final.pth'
    cp = torch.load(checkpoint_file, map_location=device, weights_only=False)

    refinement_star_file = cp['refinement_directory']
    if refinement_star_file.suffix == '.star':
        pass
    else:
        refinement_star_file = refinement_star_file / 'run_data.star'

    circular_mask_thickness = 20
    encoder_h1 = cp['encoder_half1']
    decoder_h1 = cp['decoder_half1']
    encoder_h2 = cp['encoder_half2']
    decoder_h2 = cp['decoder_half2']

    poses = cp['poses']

    relion_dataset = RelionDataset(
        path=refinement_star_file.resolve(),
        circular_mask_thickness=circular_mask_thickness,
        particle_diameter=None,
    )
    dataset = relion_dataset.make_particle_dataset()
    diameter_ang = relion_dataset.particle_diameter
    ang_pix = relion_dataset.pixel_spacing_angstroms

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
    indices_h1 = cp['indices_half1'].cpu().numpy()
    indices_h2 = np.asarray(list(set(range(len(dataset))) - set(list(indices_h1))))

    decoder_h1.p2i.device = device
    decoder_h1.projector.device = device
    decoder_h1.image_smoother.device = device
    decoder_h1.p2v.device = device
    decoder_h1.device = device
    decoder_h1.to(device)

    decoder_h2.p2i.device = device
    decoder_h2.projector.device = device
    decoder_h2.image_smoother.device = device
    decoder_h2.p2v.device = device
    decoder_h2.device = device
    decoder_h2.to(device)

    dataset_h1 = torch.utils.data.Subset(dataset, indices_h1)
    dataset_h2 = torch.utils.data.Subset(dataset, indices_h2)
    dataloader_h1 = DataLoader(
        dataset=dataset_h1,
        batch_size=batch_size,
        num_workers=n_workers,
        shuffle=False,
        pin_memory=True
    )
    dataloader_h2 = DataLoader(
        dataset=dataset_h2,
        batch_size=batch_size,
        num_workers=n_workers,
        shuffle=False,
        pin_memory=True
    )

    batch_h1 = next(iter(dataloader_h1))
    batch_h2 = next(iter(dataloader_h2))

    data_preprocessor_h1 = ParticleImagePreprocessor()
    data_preprocessor_h2 = ParticleImagePreprocessor()
    data_preprocessor_h1.initialize_from_stack(
        stack=batch_h1['image'],
        circular_mask_radius=diameter_ang / (2 * ang_pix),
        circular_mask_thickness=circular_mask_thickness / ang_pix)
    data_preprocessor_h2.initialize_from_stack(
        stack=batch_h2['image'],
        circular_mask_radius=diameter_ang / (2 * ang_pix),
        circular_mask_thickness=circular_mask_thickness / ang_pix)

    with torch.no_grad():
        latent_space_h1, _, _, _ = compute_latent_space_and_colors(
            encoder_h1, decoder_h1, dataloader_h1, poses, data_preprocessor_h1, indices_h1, reduce_by_deformation
        )
        latent_space_h2, _, _, _ = compute_latent_space_and_colors(
            encoder_h2, decoder_h2, dataloader_h2, poses, data_preprocessor_h2, indices_h2, reduce_by_deformation
        )

    # Merge latent spaces
    latent_space = torch.vstack((latent_space_h1, latent_space_h2)).cpu().numpy()
    indices = np.hstack((indices_h1, indices_h2))
    sorting_inds = np.argsort(indices)

    # Save latent space
    np.save(os.path.join(output_directory, "latent_vectors.npy"), latent_space[sorting_inds])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_directory", required=True, type=str)
    parser.add_argument("--checkpoint_file", required=True, type=str)
    parser.add_argument("--gpu_id", required=True, type=int, default=0)
    args = parser.parse_args()

    encode_latent_space(**args.__dict__)
