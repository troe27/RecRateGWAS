#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/tilman/Apis/gmmat-test"
DATA_DIR="$ROOT/real_data"
R_VERSION="${R_VERSION:-4.5.3}"
R_LIB_PROJECT="${R_LIB_PROJECT:-$ROOT/r-lib/$R_VERSION}"

export R_LIBS_USER="$R_LIB_PROJECT"

rig run -r "$R_VERSION" -f "$ROOT/BEEWAS/run_gmmat_gwas.R" \
  "$DATA_DIR/CO_per_bp_corrected_genome_size_filtered_rm_outliers.gmmat.union.tsv" \
  "$DATA_DIR/Q_D_merged_sorted_rm_bad_indv_mac1.recode.onlyqueens.union.vcf" \
  "$DATA_DIR/LDAKgmat_unimp.named.union.tsv" \
  "$DATA_DIR/gmmat_real"
