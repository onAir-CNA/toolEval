#!/usr/bin/env python3
"""
build_matrix.py – costruisce la count matrix per AtaCNV
Uso: python3 build_matrix.py FRAG_FILE BARCODES_FILE BIN_CSV OUT_RDS
"""
import sys, gzip
import numpy as np
import pandas as pd
import pyreadr

FRAG_FILE, BARCODES_FILE, BIN_FILE, OUT_RDS = sys.argv[1:5]
CHUNK_SIZE = 5_000_000
K          = 1_000_000

# Bin reference
print("Caricamento bin_info_hg38...")
bin_df = pd.read_csv(BIN_FILE)
bin_df.columns = [c.lower() for c in bin_df.columns]
bin_df["start"] = bin_df["start"].astype(int)
bin_df["end"]   = bin_df["end"].astype(int)
bin_df["bin_key"] = bin_df["chr"] + "_" + ((bin_df["start"] // K) * K).astype(str)
key_to_colidx = {row["bin_key"]: idx for idx, row in bin_df.iterrows()}
valid_chroms  = set(bin_df["chr"].unique())
print(f"Bin: {len(bin_df)} | Cromosomi validi: {len(valid_chroms)}")

# Barcodes
with gzip.open(BARCODES_FILE, "rt") as fh:
    barcodes = [l.strip() for l in fh if l.strip()]
bc_set = set(barcodes)
bc_idx = {bc: i for i, bc in enumerate(barcodes)}
n_cells, n_bins = len(barcodes), len(bin_df)
print(f"Barcodes: {n_cells}")

# Parsing
counts = np.zeros((n_cells, n_bins), dtype=np.int32)
total_read = total_used = chunk_n = 0

print("Parsing fragment file (chunked)...")
reader = pd.read_csv(
    FRAG_FILE, sep="\t", header=None, comment="#",
    names=["chr","start","end","barcode","frag_count"],
    dtype={"chr":str,"start":np.int32,"end":np.int32,
           "barcode":str,"frag_count":np.int32},
    chunksize=CHUNK_SIZE, engine="c",
)
for chunk in reader:
    chunk_n += 1; total_read += len(chunk)
    chunk = chunk[chunk["barcode"].isin(bc_set) &
                  chunk["chr"].isin(valid_chroms)].copy()
    if chunk.empty: continue
    bin_1mb  = (chunk["start"].values // K) * K
    bin_keys = np.array([f"{c}_{s}" for c, s in zip(chunk["chr"].values, bin_1mb)])
    col_idxs = np.array([key_to_colidx.get(k, -1) for k in bin_keys], dtype=np.int32)
    row_idxs = np.array([bc_idx[b] for b in chunk["barcode"].values], dtype=np.int32)
    valid = col_idxs >= 0
    np.add.at(counts, (row_idxs[valid], col_idxs[valid]), 1)
    total_used += valid.sum()
    if chunk_n % 5 == 0:
        print(f"  chunk {chunk_n}: {total_read:,} letti, {total_used:,} mappati")

print(f"\nTOTALE: {total_read:,} letti | {total_used:,} mappati ({100*total_used/max(total_read,1):.1f}%)")
print(f"Celle non-zero: {(counts.sum(axis=1)>0).sum()} / {n_cells}")
print(f"Bin  non-zero:  {(counts.sum(axis=0)>0).sum()} / {n_bins}")
print(f"Somma conteggi: {counts.sum():,}")

col_names = [f"{row['chr']}_{row['start']}_{row['end']}" for _, row in bin_df.iterrows()]
count_df  = pd.DataFrame(counts, index=barcodes, columns=col_names)
pyreadr.write_rds(OUT_RDS, count_df)
print(f"Salvato: {OUT_RDS}")
