#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path('/home/tilman/Apis/gmmat-test')
REAL_DATA = ROOT / 'real_data'
FAM_IN = REAL_DATA / 'all_queen_phased.fam'
GRM_RAW_IN = REAL_DATA / 'LDAKgmat_unimp.grm.raw'
VCF_IN = REAL_DATA / 'Q_D_merged_sorted_rm_bad_indv_mac1.recode.onlyqueens.vcf'
GRM_OUT = REAL_DATA / 'LDAKgmat_unimp.named.tsv'
GRM_VCF_ORDER_OUT = REAL_DATA / 'LDAKgmat_unimp.named.vcf_order.tsv'
GRM_MAP_OUT = REAL_DATA / 'LDAKgmat_unimp.grm_id_map.tsv'


def read_vcf_samples(vcf_path: Path) -> list[str]:
    with vcf_path.open() as handle:
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
    fam['expected_sample_id'] = fam['family_id'] + '_queen'
    mismatch = fam.loc[fam['sample_id'] != fam['expected_sample_id'], ['family_id', 'sample_id']]
    if not mismatch.empty:
        raise ValueError(
            'FAM sample IDs do not match the expected <family_id>_queen pattern. '
            f'Examples: {mismatch.head().to_dict(orient="records")}'
        )
    return fam


def main() -> None:
    fam = load_fam(FAM_IN)
    grm = pd.read_csv(GRM_RAW_IN, sep=r'\s+', header=None)

    if grm.shape[0] != grm.shape[1]:
        raise ValueError(f'Raw GRM is not square: {grm.shape}')
    if grm.shape[0] != len(fam):
        raise ValueError(
            f'Raw GRM dimension {grm.shape[0]} does not match FAM rows {len(fam)}'
        )

    queen_ids = fam['family_id'].tolist()
    named_grm = pd.DataFrame(grm.to_numpy(), index=queen_ids, columns=queen_ids)
    named_grm.to_csv(GRM_OUT, sep='\t', index=True, index_label=False)

    grm_map = fam[['grm_id', 'family_id', 'sample_id']].rename(
        columns={'family_id': 'queen_id', 'sample_id': 'fam_sample_id'}
    )
    grm_map.to_csv(GRM_MAP_OUT, sep='\t', index=False)

    if VCF_IN.exists():
        vcf_samples = read_vcf_samples(VCF_IN)
        missing_from_grm = [sample for sample in vcf_samples if sample not in named_grm.index]
        if missing_from_grm:
            raise ValueError(
                'Some VCF samples are missing from the named GRM: '
                f'{missing_from_grm[:20]}'
            )
        grm_vcf_order = named_grm.loc[vcf_samples, vcf_samples]
        grm_vcf_order.to_csv(GRM_VCF_ORDER_OUT, sep='\t', index=True, index_label=False)
        print(f'Wrote VCF-ordered GRM: {GRM_VCF_ORDER_OUT}')
        print(f'VCF queen samples retained in GRM: {len(vcf_samples)}')
    else:
        print(f'Skipped VCF-ordered GRM because {VCF_IN} does not exist yet')

    print(f'Wrote named GRM: {GRM_OUT}')
    print(f'Wrote GRM/FAM ID map: {GRM_MAP_OUT}')
    print(f'Named GRM size: {named_grm.shape[0]} x {named_grm.shape[1]}')


if __name__ == '__main__':
    main()
