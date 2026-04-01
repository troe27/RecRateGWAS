# rig-Based R Setup For GMMAT

This project now includes a reproducible local R bootstrap based on `rig`, instead of relying on a conda/mamba R stack.

## Files

- `010_install_rig_ubuntu.sh`
  - installs `rig` from the official Ubuntu/Debian apt repository
- `011_setup_r_gmmat_env.sh`
  - installs an R version via `rig`
  - sets it as default
  - creates a project-local package library
  - installs the R packages needed for this project
- `install_r_packages.R`
  - installs the required R packages into the selected library path, including `ashr`

## Why this route

`GMMAT` repeated-measures fitting appears to be failing in the current local R stack, while simpler mixed models work on the same phenotype. A clean `rig`-managed R install is the fastest way to separate a package-stack problem from a data problem.

The upstream `rig` documentation says Ubuntu 20.04 is supported, and recommends installing `rig` from its Debian/Ubuntu repository, then using `rig add release` to install R. This project script avoids `rig resolve` during setup because that command can depend on a live network lookup.
Sources:
- https://github.com/r-lib/rig
- https://github.com/r-lib/rig/releases

## 1. Install rig

Run:

```bash
./010_install_rig_ubuntu.sh
```

This requires:
- Ubuntu/Debian
- `sudo`
- network access

## 2. Install a clean R and project packages

Run:

```bash
./011_setup_r_gmmat_env.sh
```

Defaults:
- `R_SELECTOR=release`
- project library root: `./r-lib`

You can override these, e.g.:

```bash
R_SELECTOR=4.5.1 ./011_setup_r_gmmat_env.sh
PROJECT_LIB_ROOT=/some/path/r-lib ./011_setup_r_gmmat_env.sh
```

## 3. Use the clean R install

For one-off commands:

```bash
R_LIBS_USER="$(pwd)/r-lib/$(rig list | awk '/^[*] / {print $2; exit}')" \
  rig run release -- Rscript debug_gmmat_gwas.R \
  real_data/CO_per_bp_corrected_genome_size_filtered_rm_outliers.gmmat.union.tsv \
  real_data/concat_queen_biallelic_rm_bad_indv_mac1_RIC.onlyqueens.union.vcf \
  real_data/LDAKgmat_unimp.named.union.norm_ridge.tsv
```

Or make the installed version your default and export the library path:

```bash
export R_LIBS_USER=/path/to/project/r-lib/<resolved-version>
rig default <resolved-version>
```

## Installed R packages

The bootstrap currently installs:
- `Matrix`
- `nlme`
- `data.table`
- `Rcpp`
- `RcppArmadillo`
- `GMMAT`
- `ashr`

## Notes

- `011_setup_r_gmmat_env.sh` is safe to rerun.
- The package library is version-specific: `r-lib/<R-version>`.
- Keeping the library inside the project makes it easy to reproduce on another machine.
