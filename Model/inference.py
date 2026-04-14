
"""
Inference Script for UDIP-FA Model.

This script loads a trained AutoEncoder model and generates latent embeddings (endophenotypes)
for a given dataset of MRI images.

Usage:
    python inference.py --input_csv /path/to/data.csv --checkpoint /path/to/model.ckpt --output_dir /path/to/output

Author: Antigravity (Refactored)
"""

import argparse
import os
import pickle
import sys

import pandas as pd
import torch
from monai import transforms
from tqdm import tqdm

# Import model definition
# Ensure the directory containing model.py is in the python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from model import engine_AE
    from dataset import load_and_normalize_nifti
except ImportError:
    # If running from parent directory
    from Model.model import engine_AE
    from Model.dataset import load_and_normalize_nifti

# Define default transforms using MONAI
# These transforms add a channel dimension and convert the numpy array to a torch tensor.
transforms_monai = transforms.Compose(
    [
        transforms.AddChannel(),
        transforms.ToTensor(),
    ]
)


class InferenceDataset(torch.utils.data.Dataset):
    """
    Dataset class for inference.
    
    Reads a CSV file containing paths to MRI images, loads them, 
    normalizes them, and applies transforms.
    
    Returns:
        img (Tensor): Preprocessed image tensor.
        name (str): ID or filename of the image, used for tracking results.
    """

    def __init__(self, datafile, modality_column, transforms=None):
        """
        Args:
            datafile (str): Path to the CSV file.
            modality_column (str): Column name in CSV containing the image paths.
            transforms (callable, optional): Transforms to apply to the image.
        """
        self.dataframe = pd.read_csv(datafile)
        
        if modality_column not in self.dataframe.columns:
            raise ValueError(f"Column '{modality_column}' not found in CSV file.")
            
        self.image_paths = self.dataframe[modality_column].tolist()
        self.transforms = transforms

    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        img_path = self.image_paths[idx]
        
        img_data, _ = load_and_normalize_nifti(img_path)

        # Apply transforms (Add channel, ToTensor)
        if self.transforms:
            img_data = self.transforms(img_data)

        # Ensure correct type
        img_data = img_data.type(torch.float)

        return img_data, img_path


def main(args):
    """
    Main inference function.
    """
    # 1. Setup Device
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # 2. Setup Output Directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 3. Prepare Data
    print(f"Loading data from {args.input_csv}...")
    dataset = InferenceDataset(
        datafile=args.input_csv,
        modality_column=args.modality_col,
        transforms=transforms_monai,
    )
    
    dataloader = torch.utils.data.DataLoader(
        dataset, 
        batch_size=args.batch_size, 
        pin_memory=True if device.type == 'cuda' else False,
        num_workers=args.num_workers, 
        shuffle=False
    )
    print(f"Dataset size: {len(dataset)}")

    # 4. Load Model
    print("Initializing model...")
    model = engine_AE()
    model = model.to(device)

    print(f"Loading checkpoint from {args.checkpoint}...")
    if not os.path.exists(args.checkpoint):
        raise FileNotFoundError(f"Checkpoint not found at {args.checkpoint}")
        
    checkpoint = torch.load(args.checkpoint, map_location=device)
    
    # Handle state dict compatibility
    # The saved checkpoint might contain Batch Normalization stats (running_mean, running_var)
    # or other tracking stats that are not needed or compatible if we strictly match the model architecture.
    # Specifically, if the model architecture in code was changed (e.g., from BN to InstanceNorm or customized layers),
    # we need to filter out these keys.
    state_dict = checkpoint["state_dict"] if "state_dict" in checkpoint else checkpoint
    
    new_state_dict = {
        k: v for k, v in state_dict.items() 
        if "running_mean" not in k 
        and "running_var" not in k 
        and "num_batches_tracked" not in k
    }
    
    # Load weights
    # strict=False allows loading even if some keys are missing or extra keys exist, 
    # which is useful when dealing with these specific normalization layer changes.
    model.load_state_dict(new_state_dict, strict=False)
    model.eval()
    print("Model loaded successfully.")

    # 5. Run Inference
    embeddings = []
    image_names = []
    
    print("Starting inference...")
    with torch.no_grad():
        for batch_imgs, batch_names in tqdm(dataloader, desc="Processing Batches"):
            batch_imgs = batch_imgs.to(device)
            
            # Forward pass
            # model returns (reconstruction, latent_code)
            _, latent_code = model(batch_imgs)
            
            # Move to CPU and store
            embeddings.extend(latent_code.cpu().numpy())
            image_names.extend(batch_names)

    # 6. Save Results
    output_pkl_features = os.path.join(args.output_dir, f"{args.prefix}_features.pkl")
    output_pkl_names = os.path.join(args.output_dir, f"{args.prefix}_names.pkl")
    
    print(f"Saving features to {output_pkl_features}...")
    with open(output_pkl_features, 'wb') as f:
        pickle.dump(embeddings, f)
        
    print(f"Saving names to {output_pkl_names}...")
    with open(output_pkl_names, 'wb') as f:
        pickle.dump(image_names, f)

    print("Inference completed successfully!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run inference for UDIP-FA AutoEncoder.")
    
    # Data arguments
    parser.add_argument("--input_csv", type=str, required=True, 
                        help="Path to the input CSV file containing image paths.")
    parser.add_argument("--modality_col", type=str, default="mri_names", 
                        help="Column name in the CSV file that contains the file paths.")
    
    # Model arguments
    parser.add_argument("--checkpoint", type=str, required=True, 
                        help="Path to the model checkpoint (.ckpt or .pth).")
    
    # Output arguments
    parser.add_argument("--output_dir", type=str, default="./results", 
                        help="Directory to save the output pickle files.")
    parser.add_argument("--prefix", type=str, default="output", 
                        help="Prefix for the output filenames (e.g., 'T1_128').")
    
    # Runtime arguments
    parser.add_argument("--batch_size", type=int, default=1, 
                        help="Batch size for inference.")
    parser.add_argument("--num_workers", type=int, default=4, 
                        help="Number of data loading workers.")
    parser.add_argument("--device", type=str, default="cuda:0", 
                        help="Device to run inference on (e.g., 'cuda:0', 'cpu').")

    args = parser.parse_args()
    
    main(args)
