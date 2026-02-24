import sys
import argparse
from cyvcf2 import VCF, Writer


def main():
    parser = argparse.ArgumentParser(
        prog="hapgq", description="Compute GQ in haploids VCFs"
    )
    parser.add_argument("InputVCF")
    parser.add_argument("OutputVCF")
    args = parser.parse_args()

    vcfin = args.InputVCF
    vcfou = args.OutputVCF

    vcf = VCF(vcfin)
    w = Writer(vcfou, vcf)

    for v in vcf:
        newgq = abs(v.gt_phred_ll_homalt - v.gt_phred_ll_homref)
        v.set_format("GQ", newgq)
        w.write_record(v)
    w.close()
    vcf.close()


if __name__ == "__main__":
    main()
