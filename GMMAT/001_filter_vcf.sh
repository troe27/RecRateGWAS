#!/usr/bin/env bash
set -euo pipefail

BCFTOOLS="/home/tilman/miniforge3/bin/bcftools"
INPUT_VCF="/home/tilman/Apis/gmmat-test/real_data/Q_D_merged_sorted_rm_bad_indv_mac1.recode.vcf.gz"
QUEEN_LIST="/home/tilman/Apis/gmmat-test/real_data/queen_list.txt"
INPUT_DIR="$(dirname "$INPUT_VCF")"
INPUT_BASENAME="$(basename "$INPUT_VCF")"
INPUT_STEM="${INPUT_BASENAME%.vcf.gz}"
OUTPUT_VCF="$INPUT_DIR/${INPUT_STEM}.onlyqueens.vcf"
SAMPLE_LIST="$INPUT_DIR/${INPUT_STEM}.onlyqueens.samples.txt"

if [[ ! -x "$BCFTOOLS" ]]; then
  echo "Error: bcftools not found at $BCFTOOLS" >&2
  exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "Error: input VCF not found at $INPUT_VCF" >&2
  exit 1
fi

if [[ ! -f "$QUEEN_LIST" ]]; then
  echo "Error: queen list not found at $QUEEN_LIST" >&2
  exit 1
fi

"$BCFTOOLS" query -l "$INPUT_VCF" | grep -Fxf "$QUEEN_LIST" > "$SAMPLE_LIST"

if [[ ! -s "$SAMPLE_LIST" ]]; then
  echo "Error: no queen-list samples were found in $INPUT_VCF" >&2
  exit 1
fi

"$BCFTOOLS" view \
  --samples-file "$SAMPLE_LIST" \
  --output-type v \
  --output-file "$OUTPUT_VCF" \
  "$INPUT_VCF"

echo "Wrote sample list: $SAMPLE_LIST"
echo "Wrote filtered VCF: $OUTPUT_VCF"
