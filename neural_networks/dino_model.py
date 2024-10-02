"""
From https://github.com/RobvanGastel/dinov2-finetune
"""

from dino_finetune import LinearClassifier, FPNDecoder
import math
from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


class DINOV2FeatureExtractor(nn.Module):
    def __init__(
        self,
        encoder,
        decoder,
        r: int = 3,
        alpha=16,
        use_lora: bool = False,
    ):
        """The DINOv2 encoder-decoder model for finetuning to downstream tasks.

        Args:
            encoder (nn.Module): The ViT encoder model loaded with the DINOv2 model weights.
            r (int, optional): The rank parameter of the LoRA weights. Defaults to 3.
            emb_dim (int, optional): The embedding dimension of the encoder. Defaults to 1024.
            n_classes (int, optional): The number of classes to output. Defaults to 1000.
            use_lora (bool, optional): Determines whether to use LoRA. Defaults to False.
        """
        super().__init__()
        assert r > 0

        self.use_lora = use_lora

        # Number of previous layers to use as input
        self.inter_layers = 4

        self.encoder = encoder
        for param in self.encoder.parameters():
            param.requires_grad = False

        # Decoder
        # Patch size is given by (490/14)**2 = 35 * 35

        self.decoder = decoder
        # Add LoRA layers to the encoder
        if self.use_lora:
            self.lora_layers = list(range(len(self.encoder.blocks)))
            self.w_a = []
            self.w_b = []

            for i, block in enumerate(self.encoder.blocks):
                if i not in self.lora_layers:
                    continue
                w_qkv_linear = block.attn.qkv
                dim = w_qkv_linear.in_features

                w_a_linear_q, w_b_linear_q = self._create_lora_layer(dim, r)
                w_a_linear_v, w_b_linear_v = self._create_lora_layer(dim, r)

                self.w_a.extend([w_a_linear_q, w_a_linear_v])
                self.w_b.extend([w_b_linear_q, w_b_linear_v])

                block.attn.qkv = LoRA(
                    w_qkv_linear,
                    w_a_linear_q,
                    w_b_linear_q,
                    w_a_linear_v,
                    w_b_linear_v,
                    alpha=alpha,
                )
            self._reset_lora_parameters()

    def _create_lora_layer(self, dim: int, r: int):
        w_a = nn.Linear(dim, r, bias=False)
        w_b = nn.Linear(r, dim, bias=False)
        return w_a, w_b

    def _reset_lora_parameters(self) -> None:
        for w_a in self.w_a:
            nn.init.kaiming_uniform_(w_a.weight, a=math.sqrt(5))
        for w_b in self.w_b:
            nn.init.zeros_(w_b.weight)

    def forward(self, x: torch.Tensor) -> torch.Tensor:

        feature = self.encoder.forward(x)

        # get the patch embeddings - so we exclude the CLS token
        # patch_embeddings = feature["x_norm_patchtokens"]
        logits = self.decoder(feature)

        # logits = F.interpolate(
        #     logits,
        #     size=x.shape[2:],
        #     mode="bilinear",
        #     align_corners=False,
        # )
        return logits

    def save_parameters(self, filename: str) -> None:
        """Save the LoRA weights and decoder weights to a .pt file

        Args:
            filename (str): Filename of the weights
        """
        w_a, w_b = {}, {}
        if self.use_lora:
            w_a = {f"w_a_{i:03d}": self.w_a[i].weight for i in range(len(self.w_a))}
            w_b = {f"w_b_{i:03d}": self.w_b[i].weight for i in range(len(self.w_a))}

        decoder_weights = self.decoder.state_dict()
        torch.save({**w_a, **w_b, **decoder_weights}, filename)

    def load_parameters(self, filename: str) -> None:
        """Load the LoRA and decoder weights from a file

        Args:
            filename (str): File name of the weights
        """
        state_dict = torch.load(filename)

        # Load the LoRA parameters
        if self.use_lora:
            for i, w_A_linear in enumerate(self.w_a):
                saved_key = f"w_a_{i:03d}"
                saved_tensor = state_dict[saved_key]
                w_A_linear.weight = nn.Parameter(saved_tensor)

            for i, w_B_linear in enumerate(self.w_b):
                saved_key = f"w_b_{i:03d}"
                saved_tensor = state_dict[saved_key]
                w_B_linear.weight = nn.Parameter(saved_tensor)

        # Load decoder parameters
        decoder_head_dict = self.decoder.state_dict()
        decoder_head_keys = [k for k in decoder_head_dict.keys()]
        decoder_state_dict = {k: state_dict[k] for k in decoder_head_keys}

        self.decoder.load_state_dict(decoder_state_dict)


# class LoRA(nn.Module):
#     """Low-Rank Adaptation for the for Query (Q), Key (Q), Value (V) matrices"""

#     def __init__(
#         self,
#         qkv: nn.Module,
#         linear_a_q: nn.Module,
#         linear_b_q: nn.Module,
#         linear_a_v: nn.Module,
#         linear_b_v: nn.Module,
#         alpha: float = 1.0,
#     ):
#         super().__init__()
#         self.qkv = qkv
#         self.linear_a_q = linear_a_q
#         self.linear_b_q = linear_b_q
#         self.linear_a_v = linear_a_v
#         self.linear_b_v = linear_b_v
#         self.dim = qkv.in_features
#         self.w_identity = torch.eye(self.dim)
#         self.alpha = alpha

#     def forward(self, x) -> torch.Tensor:
#         # Compute the original qkv
#         qkv = self.qkv(x)  # Shape: (B, N, 3 * org_C)

#         # Compute the new q and v components
#         new_q = self.linear_b_q(self.linear_a_q(x))
#         new_v = self.linear_b_v(self.linear_a_v(x))

#         # Add new q and v components to the original qkv tensor
#         qkv[:, :, : self.dim] += new_q
#         qkv[:, :, -self.dim :] += new_v

#         return qkv


class LoRA(nn.Module):
    """Low-Rank Adaptation for the Query (Q), Key (K), and Value (V) matrices"""

    def __init__(
        self,
        qkv: nn.Module,
        linear_a_q: nn.Module,
        linear_b_q: nn.Module,
        linear_a_v: nn.Module,
        linear_b_v: nn.Module,
        alpha: float = 16,  # Scaling factor for the adaptation
    ):
        super().__init__()
        self.qkv = qkv
        self.linear_a_q = linear_a_q
        self.linear_b_q = linear_b_q
        self.linear_a_v = linear_a_v
        self.linear_b_v = linear_b_v
        self.dim = qkv.in_features
        self.alpha = alpha  # Store alpha parameter

    def forward(self, x) -> torch.Tensor:
        # Compute the original qkv
        qkv = self.qkv(x)  # Shape: (B, N, 3 * org_C)

        # Compute the new q and v components using low-rank adaptation
        new_q = self.linear_b_q(self.linear_a_q(x))
        new_v = self.linear_b_v(self.linear_a_v(x))

        # Scale the new q and v components by alpha and rank
        new_q = (self.alpha / self.linear_a_q.out_features) * new_q
        new_v = (self.alpha / self.linear_a_v.out_features) * new_v

        # Instead of clone, use torch.add_ to modify the original qkv safely
        qkv[:, :, : self.dim].add_(new_q)
        qkv[:, :, -self.dim :].add_(new_v)

        return qkv
