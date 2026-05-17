#!/usr/bin/env python3
"""
CONGASp CNV Inference – DV90 scATAC-seq
=========================================
Repo:    https://github.com/caravagnalab/CONGASp
Paper:   Patruno et al., PLOScompbio 2023

Input: count matrix già costruita da build_matrix.py (PCR-corrected)
       /sharedFolder/Results/AtaCNV/DV90/count.rds
       celle x bin, frammenti unici (colonna 5 ignorata)
"""

import os
import sys
import logging
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
import scipy.sparse as sp
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pyro.infer import SVI, Trace_ELBO
from pyro.optim import ClippedAdam

from congas.Interface import Interface
from congas.models.LatentCategorical import LatentCategorical

warnings.filterwarnings("ignore", category=UserWarning)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURAZIONE
# ══════════════════════════════════════════════════════════════════════════════
COUNT_RDS  = "/sharedFolder/Results/AtaCNV/DV90/count.rds"
BARCODES_F = "/sharedFolder/Data/filtered_feature_bc_matrix/barcodes.tsv.gz"
OUTDIR     = Path("/sharedFolder/Results/CONGAS_python/DV90")

CHROMS       = ["chr{}".format(i) for i in range(1, 23)]
K_RANGE      = [2, 3, 4]
N_STEPS      = 500
PATIENCE     = 10
LR           = 0.05
SEED         = 42
MIN_FRAGS    = 500
MIN_BIN_OCC  = 0.05
HIDDEN_DIM   = 5
LAMBDA       = 0.0

# ══════════════════════════════════════════════════════════════════════════════
# LOGGING
# ══════════════════════════════════════════════════════════════════════════════
OUTDIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
    handlers=[
        logging.FileHandler(str(OUTDIR / "run.log"), mode="w"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("congas_dv90")
log.info("=" * 65)
log.info("  CONGASp – DV90 scATAC-seq CNV inference (v2, PCR-corrected)")
log.info("=" * 65)
log.info("Input count matrix: {}".format(COUNT_RDS))
log.info("Output dir: {}".format(OUTDIR))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 – CARICA COUNT MATRIX DA AtaCNV (già PCR-corrected)
# ══════════════════════════════════════════════════════════════════════════════
log.info("[1/5] Caricamento count matrix da AtaCNV ...")

cr       = pyreadr.read_r(COUNT_RDS)
count_df = list(cr.values())[0]   # DataFrame celle x bin

# Barcodes dai rownames
barcodes_ordered = count_df.index.tolist()
n_cells          = len(barcodes_ordered)
log.info("      Celle: {}".format(n_cells))

# Trasponi a bin x celle (formato CONGASp)
count_mat = count_df.values.T.astype(np.float32)  # bin x celle

# Ricostruisci bins_df dai nomi colonna tipo "chr1_0_1000000"
bin_records = []
for col in count_df.columns:
    parts = col.split("_")
    bin_records.append({
        "chr":   parts[0],
        "start": int(float(parts[1])),
        "end":   int(float(parts[2])),
    })
bins_df = pd.DataFrame(bin_records)
n_bins  = len(bins_df)
log.info("      Bin: {}".format(n_bins))

count_csr    = sp.csr_matrix(count_mat)
barcodes_arr = np.array(barcodes_ordered, dtype=object)

log.info("      Matrice caricata: {} bin x {} celle".format(n_bins, n_cells))
log.info("      Somma totale conteggi: {:,}".format(int(count_mat.sum())))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 – QC FILTRO CELLE E BIN
# ══════════════════════════════════════════════════════════════════════════════
log.info("[2/5] Quality filtering celle e bin ...")

cell_totals = np.asarray(count_csr.sum(axis=0)).flatten()

# Filtro celle
cell_mask = cell_totals >= MIN_FRAGS
n_pass    = int(cell_mask.sum())
log.info("      Celle con >= {} frammenti: {} / {}".format(MIN_FRAGS, n_pass, n_cells))

# Filtro bin
count_cells_pass = count_csr[:, cell_mask]
bin_occ  = np.asarray((count_cells_pass > 0).sum(axis=1)).flatten()
min_occ  = int(n_pass * MIN_BIN_OCC)
bin_mask = bin_occ >= min_occ
log.info("      Bin con >= {:.0%} occupancy: {} / {}".format(
    MIN_BIN_OCC, bin_mask.sum(), n_bins))

# Applica maschere
cell_idx = np.where(cell_mask)[0]
bin_idx  = np.where(bin_mask)[0]

count_filt    = count_csr[bin_idx, :][:, cell_idx].toarray().astype(np.float32)
barcodes_filt = np.array(barcodes_ordered, dtype=object)[cell_idx]
bins_filt     = bins_df.iloc[bin_idx].reset_index(drop=True)

n_bins_filt  = count_filt.shape[0]
n_cells_filt = count_filt.shape[1]
log.info("      Matrice finale: {} bin x {} celle".format(n_bins_filt, n_cells_filt))

# Salva matrice filtrata
np.save(str(OUTDIR / "count_matrix_filtered.npy"), count_filt)
bins_filt.to_csv(str(OUTDIR / "bins_filtered.csv"), index=False)
pd.DataFrame({"barcode": barcodes_filt}).to_csv(
    str(OUTDIR / "barcodes_filtered.csv"), index=False)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 – DIZIONARIO INPUT CONGASp
# ══════════════════════════════════════════════════════════════════════════════
log.info("[3/5] Costruzione dizionario input CONGASp ...")

cell_totals_filt = count_filt.sum(axis=0)
mean_total       = float(cell_totals_filt.mean())
norm_factor      = (cell_totals_filt / mean_total).astype(np.float32)

data_atac_t   = torch.tensor(count_filt,  dtype=torch.float32)
norm_factor_t = torch.tensor(norm_factor, dtype=torch.float32)
pld_t         = torch.full((n_bins_filt,), 2.0, dtype=torch.float32)

theta_shape_val  = 10.0
theta_rate_val   = 0.2
theta_shape_atac = torch.full((n_bins_filt,), theta_shape_val, dtype=torch.float32)
theta_rate_atac  = torch.full((n_bins_filt,), theta_rate_val,  dtype=torch.float32)

data_dict = {
    "data_atac"        : data_atac_t,
    "norm_factor_atac" : norm_factor_t,
    "pld"              : pld_t,
    "segments"         : bins_filt,
}

log.info("      data_atac shape : {}".format(tuple(data_atac_t.shape)))
log.info("      norm_factor range: [{:.3f}, {:.3f}]".format(
    float(norm_factor.min()), float(norm_factor.max())))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 – FIT LatentCategorical per K = 2, 3, 4
# ══════════════════════════════════════════════════════════════════════════════
log.info("[4/5] Fit CONGASp (LatentCategorical) per K = {} ...".format(K_RANGE))

results = {}

for K in K_RANGE:
    log.info("")
    log.info("  ── K = {} ──────────────────────────────────────────".format(K))
    torch.manual_seed(SEED)

    model_params = {
        "K"                  : K,
        "hidden_dim"         : HIDDEN_DIM,
        "probs"              : torch.tensor([0.05, 0.25, 0.4, 0.2, 0.1]),
        "likelihood_atac"    : "NB",
        "likelihood_rna"     : "NB",
        "lambda"             : float(LAMBDA),
        "theta_shape_atac"   : theta_shape_atac,
        "theta_rate_atac"    : theta_rate_atac,
        "binom_prior_limits" : [10, 10_000],
        "equal_sizes_sd"     : False,
        "nb_size_init_atac"  : torch.full((K, n_bins_filt), 50.0, dtype=torch.float32),
        "init_probs"         : 0.7,
        "mixture"            : torch.full((K,), 1.0 / K, dtype=torch.float32),
        "CUDA"               : False,
        "multiome"           : False,
        "normal_cells"       : False,
    }

    interface = Interface(
        model     = LatentCategorical,
        optimizer = ClippedAdam,
        loss      = Trace_ELBO,
    )
    interface.initialize_model(data_dict)
    interface.set_model_params(model_params)

    log.info("      Avvio SVI ...")
    loss_trace, n_obs = interface.run(
        steps           = N_STEPS,
        param_optimizer = {"lr": LR},
        e               = 1e-3,
        patience        = PATIENCE,
        seed            = SEED,
    )

    params = interface.learned_parameters()
    ICs    = interface.calculate_ICs()

    log.info("      BIC = {:.4f}   ICL = {:.4f}   NLL = {:.4f}".format(
        float(ICs["BIC"]), float(ICs["ICL"]), float(ICs["NLL"])))

    results[K] = {"loss": loss_trace, "params": params, "ICs": ICs}

    with open(str(OUTDIR / "fit_K{}.pkl".format(K)), "wb") as fh:
        pickle.dump({
            "params"  : params,
            "ICs"     : ICs,
            "loss"    : loss_trace,
            "barcodes": barcodes_filt,
            "bins"    : bins_filt,
        }, fh)
    log.info("      Salvato -> fit_K{}.pkl".format(K))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 – MODEL SELECTION + OUTPUT
# ══════════════════════════════════════════════════════════════════════════════
log.info("")
log.info("[5/5] Model selection ...")

bic_rows = []
for K in K_RANGE:
    bic_rows.append({
        "K"  : K,
        "BIC": float(results[K]["ICs"]["BIC"]),
        "ICL": float(results[K]["ICs"]["ICL"]),
        "NLL": float(results[K]["ICs"]["NLL"]),
    })
bic_table = pd.DataFrame(bic_rows)
bic_table.to_csv(str(OUTDIR / "model_selection_BIC.csv"), index=False)
log.info("\n" + bic_table.to_string(index=False))

best_K      = int(bic_table.loc[bic_table["BIC"].idxmin(), "K"])
best_params = results[best_K]["params"]
log.info("\n  ★  Miglior K per BIC = {}".format(best_K))

# CNV calls
if "CNA" in best_params:
    cna = best_params["CNA"]
    for k in range(best_K):
        bins_filt["clone_{}".format(k)] = cna[k]
    bins_filt.to_csv(str(OUTDIR / "cnv_calls_K{}.csv".format(best_K)), index=False)
    log.info("  CNV calls -> cnv_calls_K{}.csv".format(best_K))

# Cluster assignments
if "assignment_atac" in best_params:
    cell_cluster = best_params["assignment_atac"]
    cell_df = pd.DataFrame({"barcode": barcodes_filt, "cluster": cell_cluster})
    cell_df.to_csv(str(OUTDIR / "cell_assignments_K{}.csv".format(best_K)), index=False)
    log.info("  Cluster sizes:\n{}".format(
        cell_df["cluster"].value_counts().to_string()))

# ── Plot ELBO ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 4))
for K in K_RANGE:
    loss = [x for x in results[K]["loss"] if x is not None]
    ax.plot(loss, label="K={}".format(K), alpha=0.85)
ax.set_xlabel("Step SVI"); ax.set_ylabel("ELBO loss")
ax.set_title("CONGASp – ELBO training curves (DV90)")
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(str(OUTDIR / "ELBO_loss_curves.png"), dpi=150)
plt.close(fig)

# ── Plot BIC ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(bic_table["K"], bic_table["BIC"], "o-", label="BIC", color="steelblue")
ax.plot(bic_table["K"], bic_table["ICL"], "s--", label="ICL", color="tomato")
ax.axvline(best_K, color="grey", linestyle=":", alpha=0.7,
           label="best K={}".format(best_K))
ax.set_xlabel("K"); ax.set_ylabel("Information criterion")
ax.set_title("Model selection – BIC/ICL (DV90)")
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(str(OUTDIR / "model_selection_BIC.png"), dpi=150)
plt.close(fig)

# ── Plot CNV heatmap ──────────────────────────────────────────────────────────
if "CNA" in best_params:
    cna  = best_params["CNA"]
    fig, axes = plt.subplots(best_K, 1, figsize=(16, 2.5 * best_K), sharex=True)
    if best_K == 1:
        axes = [axes]
    try:
        cmap = plt.cm.get_cmap("RdYlBu_r", HIDDEN_DIM)
    except AttributeError:
        cmap = plt.colormaps.get_cmap("RdYlBu_r").resampled(HIDDEN_DIM)

    for k, ax in enumerate(axes):
        im = ax.imshow(cna[k:k+1, :], aspect="auto", cmap=cmap,
                       vmin=1, vmax=HIDDEN_DIM, interpolation="nearest")
        ax.set_yticks([])
        ax.set_ylabel("Clone {}".format(k), fontsize=9)
        pos = 0
        xtick_pos, xtick_lbl = [], []
        for chrom in CHROMS:
            n_c = int((bins_filt["chr"] == chrom).sum())
            if n_c == 0: continue
            ax.axvline(pos + n_c - 0.5, color="black", linewidth=0.5, alpha=0.5)
            xtick_pos.append(pos + n_c // 2)
            xtick_lbl.append(chrom.replace("chr", ""))
            pos += n_c
    axes[-1].set_xticks(xtick_pos)
    axes[-1].set_xticklabels(xtick_lbl, fontsize=7, rotation=45)
    axes[-1].set_xlabel("Posizione genomica (bin 1 Mb)")
    plt.colorbar(im, ax=axes, label="Stato CN", shrink=0.6, pad=0.02,
                 ticks=list(range(1, HIDDEN_DIM + 1)))
    fig.suptitle("CONGASp CNV profile – DV90 (K={})".format(best_K), fontsize=12)
    plt.tight_layout()
    fig.savefig(str(OUTDIR / "CNV_heatmap_K{}.png".format(best_K)), dpi=150)
    plt.close(fig)

log.info("")
log.info("=" * 65)
log.info("  RUN COMPLETATO")
log.info("  Output: {}".format(OUTDIR))
for p in sorted(OUTDIR.iterdir()):
    if not p.is_dir():
        log.info("    {}".format(p.name))
log.info("=" * 65)