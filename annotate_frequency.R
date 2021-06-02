library(VariantAnnotation)

annotate_frequency <- function(vcf.o){
  ## placehoder: dummy data for test
  afs = runif(length(vcf.o))

  ## IDEA: we could only keep SVs that are rare to reduce the amount
  ## of variants to analyze in the next modules
  
  ## add field to VCF object
  info(header(vcf.o)) = rbind(info(header(vcf.o)),
                              AF=S4Vectors::DataFrame(Number='1',
                                                      Type='Float',
                                                      Description='Allele frequency'))  
  info(vcf.o)$AF = afs
    
  return(vcf.o)
}
