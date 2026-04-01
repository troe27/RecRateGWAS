#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[[1]] else Sys.getenv('R_LIBS_USER', unset = '')
if (!nzchar(lib)) {
  stop('Library path must be supplied as the first argument or via R_LIBS_USER.')
}

dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(lib, winslash = '/', mustWork = TRUE), .libPaths()))

repos <- c(CRAN = 'https://cloud.r-project.org')
pkgs <- c(
  'Matrix',
  'nlme',
  'data.table',
  'Rcpp',
  'RcppArmadillo',
  'GMMAT',
  'ashr'
)

installed <- rownames(installed.packages())
need <- setdiff(pkgs, installed)

message('Library path: ', .libPaths()[1])
message('Requested packages: ', paste(pkgs, collapse = ', '))
message('Already installed: ', paste(intersect(pkgs, installed), collapse = ', '))

if (length(need) > 0) {
  install.packages(need, repos = repos, lib = .libPaths()[1], Ncpus = max(1L, parallel::detectCores() - 1L))
} else {
  message('All requested packages are already installed.')
}

message('Final package versions:')
for (pkg in pkgs) {
  message('  ', pkg, ' ', as.character(packageVersion(pkg)))
}
