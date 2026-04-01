#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import csv
import gzip

import pandas as pd

ROOT = Path('/home/tilman/Apis/gmmat-test')
REAL_DATA = ROOT / 'real_data'
PHENO_IN = REAL_DATA / 'CO_per_bp_corrected_genome_size_filtered_rm_outliers.tsv'
FAM_IN = REAL_DATA / 'all_queen_phased.fam'
VCF_IN = REAL_DATA / 'Q_D_merged_sorted_rm_bad_indv_mac1.recode.vcf.gz'
QUEEN_LIST_IN = REAL_DATA / 'queen_list.txt'
PHENO_OUT = REAL_DATA / 'CO_per_bp_corrected_genome_size_filtered_rm_outliers.gmmat.tsv'
MAP_OUT = REAL_DATA / 'queen_grm_id_map.tsv'
BP_PER_MB = 1_000_000.0


def read_vcf_samples(vcf_path: Path) -> list[str]:
    opener = gzip.open if vcf_path.suffix == '.gz' else open
    with opener(vcf_path, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                return line.rstrip('\n').split('\t')[9:]
    raise ValueError(f'No #CHROM header found in {vcf_path}')


def load_fam(fam_path: Path) -> pd.DataFrame:
    fam = pd.read_csv(
        fam_path,
        sep=r'\s+',
        header=None,
        names=['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'fam_pheno'],
        dtype=str,
    )
    fam['grm_id'] = [str(i) for i in range(1, len(fam) + 1)]
    return fam


def main() -> None:
    pheno = pd.read_csv(PHENO_IN, sep='\t', dtype={'colony': str, 'indv': str, 'indv_same': str})
    fam = load_fam(FAM_IN)
    vcf_samples = set(read_vcf_samples(VCF_IN))
    queen_list = pd.read_csv(QUEEN_LIST_IN, header=None, names=['queen_id'], dtype=str)['queen_id'].tolist()
    queen_set = set(queen_list)

    required_pheno_cols = ['colony', 'indv', 'CO_per_bp']
    missing = [col for col in required_pheno_cols if col not in pheno.columns]
    if missing:
        raise ValueError(f'Phenotype file is missing required columns: {missing}')

    fam = fam[['family_id', 'sample_id', 'grm_id']].copy()
    fam['expected_sample_id'] = fam['family_id'] + '_queen'
    fam['sample_matches_expected'] = fam['sample_id'] == fam['expected_sample_id']
    if not fam['sample_matches_expected'].all():
        bad = fam.loc[~fam['sample_matches_expected'], ['family_id', 'sample_id']]
        raise ValueError(
            'FAM sample IDs do not match the expected <family_id>_queen pattern. '
            f'Examples: {bad.head().to_dict(orient="records")}'
        )

    fam_by_colony = fam.rename(columns={'family_id': 'colony'})[['colony', 'grm_id']].copy()
    fam_by_colony['queen_id'] = fam_by_colony['colony']

    formatted = pheno.merge(fam_by_colony, on='colony', how='left', validate='many_to_one')
    formatted = formatted[formatted['colony'].isin(queen_set)].copy()

    if formatted['queen_id'].isna().any():
        missing_colonies = sorted(formatted.loc[formatted['queen_id'].isna(), 'colony'].unique())
        raise ValueError(
            'Some phenotype colonies were not found in the FAM file: '
            f'{missing_colonies[:20]}'
        )

    missing_in_vcf = sorted(set(formatted['queen_id']) - vcf_samples)
    if missing_in_vcf:
        raise ValueError(
            'Some queen sample IDs implied by the phenotype/FAM files are missing from the VCF: '
            f'{missing_in_vcf[:20]}'
        )

    missing_in_queen_list = sorted(set(formatted['queen_id']) - queen_set)
    if missing_in_queen_list:
        raise ValueError(
            'Some formatted queen IDs are missing from queen_list.txt: '
            f'{missing_in_queen_list[:20]}'
        )

    formatted['offspring_index'] = formatted.groupby('queen_id').cumcount() + 1
    formatted['phenotype'] = formatted['CO_per_bp'].astype(float) * BP_PER_MB
    formatted = formatted[
        ['queen_id', 'grm_id', 'colony', 'indv', 'offspring_index', 'phenotype']
    ].copy()
    formatted = formatted.sort_values(['offspring_index', 'queen_id', 'indv']).reset_index(drop=True)

    formatted.to_csv(PHENO_OUT, sep='\t', index=False, quoting=csv.QUOTE_NONE)
    fam_by_colony[fam_by_colony['queen_id'].isin(queen_set)][['grm_id', 'colony', 'queen_id']].to_csv(MAP_OUT, sep='\t', index=False)

    print(f'Wrote formatted phenotype table: {PHENO_OUT}')
    print(f'Wrote queen/GRM ID map: {MAP_OUT}')
    print(f'Rows in formatted phenotype table: {len(formatted)}')
    print(f'Unique queens in formatted phenotype table: {formatted["queen_id"].nunique()}')
    print(f'Phenotype units: crossovers per megabase (CO_per_bp * {BP_PER_MB:g})')
    print(f'Phenotype mean: {formatted["phenotype"].mean():.6f}')
    print(f'Phenotype SD: {formatted["phenotype"].std(ddof=1):.6f}')


if __name__ == '__main__':
    main()
