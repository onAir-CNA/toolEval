#!/usr/bin/env Rscript
# ============================================================
# coverage_analysis.R
# AtaCNV – DV90 coverage analysis
# 4 livelli di coverage x 20 ripetizioni = 80 run
# Parallelizzato con parallel::mclapply
# ============================================================

library(AtaCNV)
library(parallelDist)
library(parallel)

# ── Percorsi ─────────────────────────────────────────────
COUNT_RDS    <- "/sharedFolder/Results/AtaCNV/DV90/count.rds"
BARCODES_F   <- "/sharedFolder/Data/filtered_feature_bc_matrix/barcodes.tsv.gz"
OUT_DIR      <- "/sharedFolder/Results/AtaCNV/coverage_analysis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Parametri ─────────────────────────────────────────────
COVERAGES  <- c(2/3, 1/2, 1/3, 1/4)
COV_LABELS <- c("cov_67", "cov_50", "cov_33", "cov_25")
N_REPS     <- 100
N_CORES    <- 8
set.seed(42)
SEEDS      <- sample(1:10000, N_REPS, replace = FALSE)

# ── Carica count matrix ───────────────────────────────────
cat("Caricamento count matrix...\n")
count_raw <- readRDS(COUNT_RDS)
barcodes  <- readLines(gzfile(BARCODES_F))
count     <- as.matrix(count_raw)
rownames(count) <- barcodes
cat(sprintf("Count matrix: %d celle x %d bin\n", nrow(count), ncol(count)))

# ── Subsampling per cella ─────────────────────────────────
subsample_count <- function(count, fraction, seed) {
  set.seed(seed)
  count_sub <- count
  for (i in 1:nrow(count)) {
    cell_total <- sum(count[i, ])
    if (cell_total == 0) next
    n_keep <- round(cell_total * fraction)
    if (n_keep == 0) next
    probs <- count[i, ] / cell_total
    count_sub[i, ] <- as.integer(rmultinom(1, size = n_keep, prob = probs))
  }
  return(count_sub)
}

# ── Singolo run AtaCNV ────────────────────────────────────
run_one <- function(args) {
  cov_frac  <- args$cov_frac
  cov_label <- args$cov_label
  rep       <- args$rep
  seed      <- args$seed
  count     <- args$count

  out_rep <- file.path(OUT_DIR, cov_label, sprintf("rep_%02d", rep))
  cr_file <- file.path(out_rep, "copy_ratio.rds")

  # Skip se già fatto
  if (file.exists(cr_file)) {
    cat(sprintf("  [SKIP] %s/rep_%02d già esistente\n", cov_label, rep))
    return(cr_file)
  }

  dir.create(file.path(out_rep, "normalize"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_rep, "cnv"),       recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("  [START] %s/rep_%02d (seed=%d, coverage=%.0f%%)\n",
              cov_label, rep, seed, cov_frac * 100))

  # Subsampling
  count_sub <- subsample_count(count, cov_frac, seed)

  # Normalize
  norm_re <- tryCatch(
    AtaCNV::normalize(
      count         = count_sub,
      genome        = "hg38",
      mode          = "none",
      output_dir    = file.path(out_rep, "normalize", ""),
      output_name   = "norm_re.rds",
      gc_correction = TRUE
    ),
    error = function(e) {
      cat(sprintf("  [ERROR normalize] %s/rep_%02d: %s\n", cov_label, rep, e$message))
      NULL
    }
  )
  if (is.null(norm_re)) return(NULL)

  # Calculate CNV
  seg_re <- tryCatch(
    AtaCNV::calculate_CNV(
      norm_count  = norm_re$norm_count,
      baseline    = norm_re$baseline,
      output_dir  = file.path(out_rep, "cnv", ""),
      output_name = "CNV_re.rds"
    ),
    error = function(e) {
      cat(sprintf("  [ERROR CNV] %s/rep_%02d: %s\n", cov_label, rep, e$message))
      NULL
    }
  )
  if (is.null(seg_re)) return(NULL)

  # Salva copy_ratio
  saveRDS(seg_re$copy_ratio, cr_file)
  cat(sprintf("  [DONE] %s/rep_%02d → copy_ratio salvato\n", cov_label, rep))
  return(cr_file)
}

# ── Costruisci lista di job ───────────────────────────────
jobs <- list()
for (ci in seq_along(COVERAGES)) {
  for (rep in 1:N_REPS) {
    jobs[[length(jobs) + 1]] <- list(
      cov_frac  = COVERAGES[ci],
      cov_label = COV_LABELS[ci],
      rep       = rep,
      seed      = SEEDS[rep],
      count     = count
    )
  }
}
cat(sprintf("\nJob totali: %d | Core: %d\n\n", length(jobs), N_CORES))

# ── Lancia in parallelo ───────────────────────────────────
results <- mclapply(jobs, run_one, mc.cores = N_CORES)

# ── Report finale ─────────────────────────────────────────
ok   <- sum(!sapply(results, is.null))
fail <- sum(sapply(results, is.null))
cat(sprintf("\n=== COMPLETATO ===\n"))
cat(sprintf("Run OK:     %d\n", ok))
cat(sprintf("Run FAILED: %d\n", fail))
cat(sprintf("Output in:  %s\n", OUT_DIR))