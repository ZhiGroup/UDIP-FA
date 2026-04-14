# imports

import argparse
import os

import pytorch_lightning as pl
import torch
from pytorch_lightning.callbacks import (
    LearningRateMonitor,
    ModelCheckpoint,
    ProgressBar,
)
from pytorch_lightning.loggers import CSVLogger, TensorBoardLogger
from torch.optim.lr_scheduler import ReduceLROnPlateau

from dataset import aedataset, transforms_monai
from model import engine_AE


class LitAutoEncoder(pl.LightningModule):
    def __init__(self, lr):
        super().__init__()
        self.save_hyperparameters()
        self.model = engine_AE()

        self.train_loss_function = torch.nn.MSELoss(reduction="none")
        self.valid_loss_function = torch.nn.MSELoss(reduction="none")

    def forward(self, x):
        return self.model(x)

    def _masked_reconstruction_loss(self, x, recon, mask, loss_fn):
        loss = loss_fn(x, recon)
        loss = loss.squeeze(1) * mask
        return loss.sum() / mask.sum()

    def training_step(self, batch, batch_idx):
        x, mask = batch
        recon, _ = self(x)
        loss = self._masked_reconstruction_loss(x, recon, mask, self.train_loss_function)
        self.log("train_loss", loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, mask = batch
        recon, _ = self(x)
        loss = self._masked_reconstruction_loss(x, recon, mask, self.valid_loss_function)
        self.log("val_loss", loss, prog_bar=True, sync_dist=True)
        return loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.hparams["lr"])
        lr_scheduler_config = {
            "scheduler": ReduceLROnPlateau(
                optimizer,
                "min",
                patience=4,
                min_lr=self.hparams["lr"] / 1000,
                factor=0.5,
            ),
            "interval": "epoch",
            "frequency": 1,
            "monitor": "val_loss",
            "strict": True,
            "name": None,
        }
        return {"optimizer": optimizer, "lr_scheduler": lr_scheduler_config}


def parse_args():
    parser = argparse.ArgumentParser(description="Train the UDIP-FA autoencoder.")
    parser.add_argument(
        "--train_csv",
        type=str,
        default="train_mixed_ethnicity.csv",
        help="CSV file containing training image paths.",
    )
    parser.add_argument(
        "--val_csv",
        type=str,
        default="val_mixed_ethnicity.csv",
        help="CSV file containing validation image paths.",
    )
    parser.add_argument(
        "--modality_col",
        type=str,
        default="T1_unbiased_linear",
        help="CSV column containing NIfTI paths.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="T1_128_7_gpus",
        help="Directory used to store checkpoints and logs.",
    )
    parser.add_argument(
        "--learning_rate",
        type=float,
        default=0.0005248074602497723,
        help="Learning rate for Adam.",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=9,
        help="Batch size for training and validation.",
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=3,
        help="Number of DataLoader workers.",
    )
    parser.add_argument(
        "--max_epochs",
        type=int,
        default=60,
        help="Maximum number of epochs.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--gpus",
        type=int,
        nargs="*",
        default=[0, 1, 2, 3, 4, 5, 6],
        help="GPU ids to use. Pass no values to run on CPU.",
    )
    return parser.parse_args()


def build_dataloader(datafile, modality, batch_size, num_workers, shuffle):
    dataset = aedataset(
        datafile=datafile,
        modality=modality,
        transforms=transforms_monai,
    )
    return torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        pin_memory=torch.cuda.is_available(),
        num_workers=num_workers,
        shuffle=shuffle,
    )


def validate_paths(args):
    for path in [args.train_csv, args.val_csv]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"CSV file not found: {path}")


def build_trainer(args):
    callbacks = [
        LearningRateMonitor(logging_interval="epoch"),
        ModelCheckpoint(
            dirpath=args.output_dir,
            monitor="val_loss",
            save_last=True,
            filename="{epoch}-{train_loss:.6f}-{val_loss:.6f}",
            save_top_k=5,
        ),
        ProgressBar(refresh_rate=2),
    ]
    loggers = [
        TensorBoardLogger(save_dir=os.path.join(args.output_dir, "tb_logs")),
        CSVLogger(save_dir=os.path.join(args.output_dir, "csv_logs")),
    ]

    trainer_kwargs = {
        "logger": loggers,
        "callbacks": callbacks,
        "log_every_n_steps": 10,
        "benchmark": True,
        "max_epochs": args.max_epochs,
    }

    if args.gpus:
        trainer_kwargs["gpus"] = args.gpus
        trainer_kwargs["sync_batchnorm"] = len(args.gpus) > 1
        if len(args.gpus) > 1:
            trainer_kwargs["accelerator"] = "ddp"

    return pl.Trainer(**trainer_kwargs)


def main():
    args = parse_args()
    validate_paths(args)
    os.makedirs(args.output_dir, exist_ok=True)
    pl.seed_everything(args.seed, workers=True)

    train_dataloader = build_dataloader(
        datafile=args.train_csv,
        modality=args.modality_col,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        shuffle=True,
    )
    val_dataloader = build_dataloader(
        datafile=args.val_csv,
        modality=args.modality_col,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        shuffle=False,
    )

    model = LitAutoEncoder(args.learning_rate)
    trainer = build_trainer(args)
    trainer.fit(
        model,
        train_dataloaders=train_dataloader,
        val_dataloaders=val_dataloader,
    )


if __name__ == "__main__":
    main()
