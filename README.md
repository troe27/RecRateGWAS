### Code and scripts used for analyses in the paper



#  Genetic determinants of intraspecific variation in crossover frequencies in the Honeybee, _Apis mellifera_

Turid Everitt<sup>1</sup>\*, Luong Hieu Trinh<sup>1</sup>\*, Demetris Taliadoros<sup>1</sup>,Tilman Rönneburg<sup>1</sup>, Anna Olsson<sup>1</sup>, Joachim R. de Miranda<sup>2</sup>, Bertrand Servin<sup>3</sup>, Matthew T. Webster<sup>1†</sup>

1) Dept. Medical Biochemistry and Microbiology, SciLifeLab, Uppsala University, Uppsala,
Sweden
2) Department of Ecology, Swedish University of Agricultural Sciences, Uppsala, Sweden
3) GenPhySE, Université de Toulouse, INRAE, ENVT, Castanet-Tolosan, France  
\* shared first-authorship  
† correspondence to matthew.webster@imbim.uu.se


## Abstract

Meiotic recombination facilitates natural selection and is necessary for correct chromosomal segregation in most sexually reproducing species. Crossover rates vary greatly both within and among species, but the determinants of this variation are not fully understood. The honeybee Apis mellifera has extremely high recombination rates. Honeybee males (drones) are haploid, which enables the distribution of crossovers to be directly estimated from the progeny of a single reproductive female (queen). Here we map crossover events in the honeybee using whole genome sequencing of 1509 drone progeny of 184 queens. This allows us to assay intra-specific variation in recombination rate and its genetic and non-genetic determinants. We estimate the average crossover rate as 23 cM/Mb, with between 22 and 88 crossovers events detected in individual offspring. We estimate 28% of this variation is additive heritable variation among queens. There is no effect of queen age or genetic background on crossover rate. A genome-wide association study identifies variation in the gene mlh1 as associated with mean crossover rate. We estimate that variation in the gene is associated with a 10% difference in crossover rate between the two homozygous genotypes at the most significant SNP. This gene has a well-established role in recombination and variation in the gene could affect crossover rates by affecting resolution of Holliday junctions as crossovers. This is the first gene discovered to be associated with recombination rate variation in an insect. Adaptive evolution of this gene could potentially underlie the extremely high recombination rates in honeybees.



currently a [preprint on bioRxiv](https://doi.org/10.64898/2026.02.03.702590)

## Software Requirements

### Software

- rig
- R
- Rscript
- python3
- bgzip
- tabix
- bcftools
- bwa
- samtools
- bamtools
- picard
- java
- gatk
- vcftools
- yapp
- plink
- plink2
- ldak6
- vep

### R Packages

- GMMAT
- ashr
- ggplot2
- plot3D
- MM4LMM
- dplyr
- scales
- xoi
- purrr
- readr
- tidyr
- stringr

### Python Packages

- pandas
- numpy
- pysam
- cyvcf2

