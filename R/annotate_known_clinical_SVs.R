library(VariantAnnotation)

annotate_known_clinical_SVs <- function(vcf.o){
  ## placehoder: dummy data for test
  clinsvs = rep(NA, length(vcf.o))
  clinsvs[sample.int(length(vcf.o), 1000)] = paste0('clin', 1:1000)

  ## add field to VCF object
  info(header(vcf.o)) = rbind(info(header(vcf.o)),
                              CLINSV=S4Vectors::DataFrame(Number='1',
                                                          Type='String',
                                                          Description='IDs of matching known clinical SVs'))  
  info(vcf.o)$CLINSV = CLINSV
    
  return(vcf.o)
}
