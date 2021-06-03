## args = commandArgs(TRUE)
## in.vcf = args[1]
## out.vcf = args[2]
## out.csv = args[3]

## test data
## downloaded from ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/
if(!file.exists('HG002_SVs_Tier1_v0.6.vcf.gz')){
  download.file('ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz', 'HG002_SVs_Tier1_v0.6.vcf.gz')
}
in.vcf = 'HG002_SVs_Tier1_v0.6.vcf.gz'
out.vcf = 'clinical-sv-annotated.vcf'
out.csv = 'clinical-sv-table.csv'

## annotation data contains pre-formatted annotation (e.g. genes, clinical svs)
## see prepare_annotation_data.R
load('annotation_data.RData')

## read VCF
suppressWarnings(suppressMessages(library(VariantAnnotation)))
vcf.o = readVcf(in.vcf)

## make sure chromosomes are in the form 'chrX'
if(all(!grepl('chr', seqlevels(vcf.o)))){
  seqlevels(vcf.o) = paste0('chr', seqlevels(vcf.o))
}

## annotate gene overlapped by SVs
source('annotate_genes.R')
vcf.o = annotate_genes(vcf.o, genc)

## annotate frequency
source('annotate_frequency.R')
vcf.o = annotate_frequency(vcf.o, gnomad)

## annotate known clinical SVs
source('annotate_known_clinical_SVs.R')
vcf.o = annotate_known_clinical_SVs(vcf.o, clinsv)

## write annotated VCF
writeVcf(vcf.o, file=out.vcf)

## write tables
svs = tibble(gene=info(vcf.o)$GENE,
             variant_id=names(vcf.o),
             chr=as.character(seqnames(vcf.o)),
             start=start(vcf.o),
             end=end(vcf.o),
             size=abs(unlist(lapply(info(vcf.o)$SVLEN, '[', 1))),
             frequency=info(vcf.o)$AF,
             svtype=info(vcf.o)$SVTYPE,
             clinsv=info(vcf.o)$CLINSV)

write.table(svs, file=out.csv, sep=',', quote=TRUE, row.names=FALSE)
