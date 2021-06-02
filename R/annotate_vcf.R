## args = commandArgs(TRUE)
## in.vcf = args[1]
## out.vcf = args[2]

## test data
## downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/
if(!file.exists('HG002_SVs_Tier1_v0.6.vcf.gz')){
  download.file('ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz', 'HG002_SVs_Tier1_v0.6.vcf.gz')
}
in.vcf = 'HG002_SVs_Tier1_v0.6.vcf.gz'
out.vcf = 'test-annotated.vcf'

## read VCF
library(VariantAnnotation)
vcf.o = readVcf(in.vcf)

## annotate gene overlapped by SVs
source('annotate_genes.R')
vcf.o = annotate_genes(vcf.o)

## annotate frequency
source('annotate_frequency.R')
vcf.o = annotate_frequency(vcf.o)

## annotate known clinical SVs
source('annotate_known_clinical_SVs.R')
vcf.o = annotate_known_clinical_SVs(vcf.o)

## writ annotated VCF
writeVcf(vcf.o, file=out.vcf)
