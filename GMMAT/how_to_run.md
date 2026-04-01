# How to run GMMAT on our data:

#### software for this specific workflow and formatting:
 - rig
 - python=3.11
 - numpy>=1.26
 - pandas>=2.2
 - pysam>=0.22
 - bcftools>=1.20 
#### required software to run GMMAT:
 - R (≥ 3.2.0)
 - [GMMAT](https://cran.r-project.org/web/packages/GMMAT/index.html)

### files in this folder:
```
├── 001_filter_vcf.sh                # dataset specific: filters drones out of the .vcf file
├── 002_format_pheno.py              # dataset specific: formats phenotype akin to GMMAT test-data, and multiplies to get CO/Mb. Also interleaves the repeat observations. 
├── 003_format_grm.py                # dataset specific: formats GRM akin to GMMAT test-data
├── 004_sample_union.py              # dataset specific: takes intersection of all datasets (GRM, VCF, Pheno) as we have more individuals than phenotypical observations
├── 005_run_gmmat.sh                 # Wrapper to run GMMAT.
├── run_gmmat_gwas.R                 # Rscript called by 005 - runs GMMAT GWAS with repeated measurements
├── gmmat_input_issues.md            # troubleshooting hints when you are having issues with your input
├── how_to_run.md                    # this document.
└── setup                            # folder containing setup scripts.
    ├── 010_install_rig_ubuntu.sh    # install rig to have a clean R independent of the system R
    ├── 011_setup_r_gmmat_env.sh     # wrapper to install packages
    ├── RIG_SETUP.md                 # LLM generated readme for the rig installation
    ├── environment.yml              # requirements for the python scripts
    └── install_r_packages.R         # r install script used by 011
```


#### Expected input:
- `concat_queen_biallelic_rm_bad_indv_mac1_RIC.vcf`
- `CO_per_bp_corrected_genome_size_filtered_rm_outliers.tsv`
- `LDAKgmat_unimp.grm.raw`
- `all_queen_phased.fam`


#### how to run: 
- run scripts in sequence: 010 -> 011 -> 001 -> 002 -> 003 -> 004 -> 005
- note that some filepaths might still be hardcoded.
