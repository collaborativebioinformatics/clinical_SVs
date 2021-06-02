library(VariantAnnotation)

annotate_genes <- function(vcf.o){
  ## placehoder: dummy data for test
  genes = rep(NA, length(vcf.o))
  genes[sample.int(length(vcf.o), 1000)] = sample(paste0('gene', 1:200), 1:1000, replace=TRUE)

  ## IDEA: we could only keep SVs that overlap a gene to reduce the amount
  ## of variants to analyze in the next modules

  ## IDEA: maybe it will be easier later if we duplicate SVs that overlap multiple genes
  ## i.e. one record for each SV-gene pair. Then the GENE field will contain only one gene
  ## name which might be easier to filter later.
  
  ## add field to VCF object
  info(header(vcf.o)) = rbind(info(header(vcf.o)),
                              GENE=S4Vectors::DataFrame(Number='1',
                                                        Type='String',
                                                        Description='Overlapped gene(s)'))  
  info(vcf.o)$genes = genes
  return(vcf.o)
}
