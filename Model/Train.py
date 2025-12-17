# imports

# PyTorch
import torch
from torch.nn import functional as F
from torch import nn
from torch.optim.lr_scheduler import ReduceLROnPlateau

# PyTorch Lightning

import pytorch_lightning as pl
from pytorch_lightning.callbacks import (
    LearningRateMonitor,
    ModelCheckpoint,
    ProgressBar,
)
from pytorch_lightning.loggers import TensorBoardLogger, CSVLogger

# Custom imports
from dataset import *
from model import engine_AE

# Model architecture and forward pass to Pytorch lightning module.


class LitAutoEncoder(pl.LightningModule):
    def __init__(self, lr):
        super().__init__()
        self.save_hyperparameters()
        self.model = engine_AE()
        
        # loss function to be used in training loop
        self.train_loss_function1 = torch.nn.MSELoss(
            size_average=None, reduce=None, reduction="none"
        )
        # loss function to be used in validation loop
        self.valid_loss_function = torch.nn.MSELoss(
            size_average=None, reduce=None, reduction="none"
        )

    # Forward function
    def forward(self, x):
        return self.model(x)

    # pytorch lightning training step
    def training_step(self, batch, batch_idx):
        # x, reg_input = batch
        x, mask = batch
        recon, _ = self(x)

        loss1 = self.train_loss_function1(x, recon)
        loss1 = loss1.squeeze(1) * mask
        loss1 = loss1.sum()
        loss1 = loss1 / mask.sum()
        # loss2 = self.train_loss_function(reg_input, reg)
        # loss = loss1 + loss2
        self.log("train_loss", loss1, prog_bar=True)
        return loss1

    # pytorch lightning validation step
    def validation_step(self, batch, batch_idx):
        x, mask = batch
        recon, _ = self(x)
        loss1 = self.valid_loss_function(x, recon)
        loss1 = loss1.squeeze(1) * mask
        loss1 = loss1.sum()
        loss1 = loss1 / mask.sum()
        # loss2 = self.train_loss_function(reg_input, reg)
        # loss = loss1 + loss2
        self.log("val_loss", loss1, prog_bar=True, sync_dist=True)
        return loss1

    # pytorch lightning optimizer configuration
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
        return {
            "optimizer": optimizer,
            "lr_scheduler": lr_scheduler_config,
        }


# defining train dataset
train_dataset = aedataset(
    datafile="train_mixed_ethnicity.csv", modality="T1_unbiased_linear", transforms=transforms_monai
    #datafile="sub2.csv", modality="T1_unbiased_linear"   
)

valid_dataset = aedataset(
    datafile="val_mixed_ethnicity.csv", modality="T1_unbiased_linear", transforms=transforms_monai
    #datafile="sub2.csv", modality="T1_unbiased_linear"
)

# defining train dataloader
train_dataloader = torch.utils.data.DataLoader(
    train_dataset, batch_size=9, pin_memory=True, num_workers=3, shuffle=True,
)

# defining validation dataloader
val_dataloader = torch.utils.data.DataLoader(
    valid_dataset, batch_size=9, pin_memory=True, num_workers=3, shuffle=False
)

# directory name to save checkpoints and metrics
dir_name = "T1_128_7_gpus"

# initiaing the model
# AE_model = LitAutoEncoder(0.0010964781961431851)
AE_model = LitAutoEncoder(0.0005248074602497723)
# learning rate monitor as using scheduler
lr_monitor = LearningRateMonitor(logging_interval="epoch")

# saving checkpoints monitoring validation loss
model_checkpoint = ModelCheckpoint(
    dirpath=dir_name,
    monitor="val_loss",
    save_last=True,
    filename="{epoch}-{train_loss:.6f}-{val_loss:.6f}",
    save_top_k=5,
)

# Loggers
tb_logger = TensorBoardLogger(save_dir=dir_name + "/tb_logs")
csv_logger = CSVLogger(save_dir=dir_name + "/csv_logs")
pb = ProgressBar(refresh_rate=2)

# main training
if __name__ == "__main__":
    trainer = pl.Trainer(
        logger=[tb_logger, csv_logger],
        # Change the number of GPUs here
        gpus=[0,1,2,3,4,5,6],
        callbacks=[lr_monitor, model_checkpoint, pb],
        sync_batchnorm=True,
        log_every_n_steps=10,
        #accelerator="dp",
        accelerator="ddp",
        benchmark=True,
        max_epochs=60,
    )

    trainer.fit(
        AE_model, train_dataloaders=train_dataloader, val_dataloaders=val_dataloader
    )
