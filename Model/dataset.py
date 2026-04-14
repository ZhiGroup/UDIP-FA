# please use affine registered MRI. Instructions in /training/readme.md

import os

import nibabel as nib
import numpy as np
import pandas as pd
import torch
from monai import transforms

transforms_monai = transforms.Compose(
    [transforms.AddChannel(), transforms.ToTensor()]
)


def load_and_normalize_nifti(img_path):
    """
    Load a NIfTI image and z-score non-zero voxels only.

    Returns:
        tuple[np.ndarray, np.ndarray]: normalized image data and non-zero mask.
    """
    if not os.path.exists(img_path):
        raise FileNotFoundError(f"Image file not found: {img_path}")

    img = nib.load(img_path).get_fdata()
    mask = img != 0

    if not mask.any():
        return np.zeros_like(img, dtype=np.float32), mask

    nonzero = img[mask]
    std = nonzero.std()
    if std == 0:
        normalized = np.zeros_like(img, dtype=np.float32)
        normalized[mask] = nonzero - nonzero.mean()
        return normalized, mask

    normalized = (img - nonzero.mean()) / std
    return normalized.astype(np.float32), mask


class aedataset(torch.utils.data.Dataset):
    def __init__(self, datafile, modality, transforms=None):
        """
        Args:
            datafile (str): CSV file containing image paths.
            modality (str): Column containing the image path for the modality of interest.
            transforms (callable): Transforms to add a channel and convert to tensor.
        Returns:
            img (torch.Tensor): Normalized image tensor.
            mask (torch.Tensor): Boolean mask excluding background.
        """
        self.datafile = pd.read_csv(datafile)
        if modality not in self.datafile.columns:
            raise ValueError(f"Column '{modality}' not found in {datafile}")

        self.unbiased_brain = self.datafile[modality].tolist()
        self.transforms = transforms if transforms is not None else transforms_monai

    def __len__(self):
        return len(self.unbiased_brain)

    def __getitem__(self, idxx=int):
        img, mask = load_and_normalize_nifti(self.unbiased_brain[idxx])
        img = self.transforms(img)

        img = img.type(torch.float)
        mask = torch.tensor(mask, dtype=torch.bool)

        return img, mask
