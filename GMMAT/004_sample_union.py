#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import subprocess

import pandas as pd

ROOT = Path('/home/tilman/Apis/gmmat-test')
REAL_DATA = ROOT / 'real_data'
BCFTOOLS = Path('/home/tilman/miniforge3/bin/bcftools')
VCF_IN = REAL_DATA / 'Q_D_merged_sorted_rm_bad_indv_mac1.recode.onlyqueens.vcf'
GRM_IN = REAL_DATA / 'LDAKgmat_unimp.named.tsv'
PHENO_IN = REAL_DATA / 'CO_per_bp_corrected_genome_size_filtered_rm_outliers.gmmat.tsv'
VCF_OUT = REAL_DATA / 'Q_D_merged_sorted_rm_bad_indv_mac1.recode.onlyqueens.union.vcf'
GRM_OUT = REAL_DATA / 'LDAKgmat_unimp.named.union.tsv'
PHENO_OUT = REAL_DATA / 'CO_per_bp_corrected_genome_size_filtered_rm_outliers.gmmat.union.tsv'
SAMPLE_LIST_OUT = REAL_DATA / 'sample_union.samples.txt'


def read_vcf_samples(vcf_path: Path) -> list[str]:
    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                return line.rstrip('\n').split('\t')[9:]
    raise ValueError(f'No #CHROM header found in {vcf_path}')


def main() -> None:
    if not BCFTOOLS.exists():
        raise FileNotFoundError(f'bcftools not found at {BCFTOOLS}')

    vcf_samples = read_vcf_samples(VCF_IN)
    grm = pd.read_csv(GRM_IN, sep='\t', index_col=0)
    pheno = pd.read_csv(PHENO_IN, sep='\t')

    if 'queen_id' not in pheno.columns:
        raise ValueError('Formatted phenotype table must contain a queen_id column')

    shared = sorted(set(vcf_samples) & set(grm.index) & set(pheno['queen_id']))
    if not shared:
        raise ValueError('No shared samples found across VCF, GRM, and phenotype table')

    pd.Series(shared).to_csv(SAMPLE_LIST_OUT, index=False, header=False)

    grm_union = grm.loc[shared, shared]
    grm_union.to_csv(GRM_OUT, sep='\t', index=True, index_label=False)

    pheno_union = pheno[pheno['queen_id'].isin(shared)].copy()
    pheno_union['queen_id'] = pd.Categorical(pheno_union['queen_id'], categories=shared, ordered=True)
    sort_cols = [col for col in ['offspring_index', 'queen_id', 'indv'] if col in pheno_union.columns]
    pheno_union = pheno_union.sort_values(sort_cols).reset_index(drop=True)
    pheno_union['queen_id'] = pheno_union['queen_id'].astype(str)
    pheno_union.to_csv(PHENO_OUT, sep='\t', index=False)

    subprocess.run(
        [
            str(BCFTOOLS),
            'view',
            '--samples-file',
            str(SAMPLE_LIST_OUT),
            '--output-type',
            'v',
            '--output-file',
            str(VCF_OUT),
            str(VCF_IN),
        ],
        check=True,
    )

    print(f'Wrote sample union list: {SAMPLE_LIST_OUT}')
    print(f'Wrote union VCF: {VCF_OUT}')
    print(f'Wrote union GRM: {GRM_OUT}')
    print(f'Wrote union phenotype table: {PHENO_OUT}')
    print(f'Shared queen samples: {len(shared)}')


if __name__ == '__main__':
    main()
