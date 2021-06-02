options(stringsAsFactors = F)
suppressWarnings(suppressMessages(library(clusterProfiler, quietly = T)))
suppressWarnings(suppressMessages(library(pathview, quietly = T)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db, quietly = T)))
suppressWarnings(suppressMessages(library(enrichplot, quietly = T)))
suppressWarnings(suppressMessages(library(DOSE, quietly = T)))
suppressWarnings(suppressMessages(library(ggnewscale, quietly = T)))
suppressWarnings(suppressMessages(library(cowplot, quietly = T)))
suppressWarnings(suppressMessages(library(tidyverse, quietly = T)))

args = commandArgs(TRUE)

if(length(args) < 3)
{
 stop("#ERROR! please supply all inputs.
\nUSAGE: Rscript ./geneFunctionalAnnotation.R <listofgene> <pvalueCutoff> <GeneIDtype>\n")
}

listofgene = args[1] #list of genes
pvalueCutoff = as.numeric(args[2]) #the minimum p value for enriched GO and pathways
GeneIDtype = as.character(args[3]) #ID should be one of: SYMBOL,REFSEQ,ENSEMBL (## case sensitive)

# listofgene = "listofENSEMBLID.txt"
# pvalueCutoff = 0.05
# GeneIDtype = "ENSEMBL"

listofgene <- read.table(listofgene, header = F)
listofgene <- listofgene[,1]

## Covert list of genes to entrezID
genelist2convertID <- function(listofgene, fromType, toType){
  id <- suppressMessages(suppressWarnings(bitr(listofgene, fromType=fromType, toType=toType, OrgDb="org.Hs.eg.db", drop = TRUE)))
  ### convert as charecter vector
  id2 = as.character(id[,2])
  ## to strip the NA values
  x <- id2[!is.na(id2)] ## entrez id
  x<- unique(x)
  return(x)
}

listEntrezID<-genelist2convertID(listofgene = listofgene, fromType = GeneIDtype, toType = 'ENTREZID')

## gene Enrichment analysis (only Disease ontology as this is clinical data)
enrichedGene <- function (listEntrezID, showCategory=20){
	endgn <- DOSE::enrichDGN(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff)
	enDO <- DOSE::enrichDO(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff, qvalueCutoff = 0.1)
	bar_plot <- barplot(enDO, showCategory=showCategory) + ggtitle("Higher level of Disease Ontology")
	dot_plot <- enrichplot::dotplot(endgn, showCategory=showCategory, font.size = 8) + ggtitle("Disease associations from DisGeNET")
	edox <- DOSE::setReadable(endgn, 'org.Hs.eg.db', 'ENTREZID')
	cnet_plot <- cnetplot(edox, colorEdge = TRUE, categorySize="pvalue") + ggtitle("Top 5 Category")
	cnet_plot <- cnet_plot+theme(plot.title = element_text(hjust=0.7))
	anno_plot <- cowplot::plot_grid(bar_plot, dot_plot, cnet_plot, label_size = 8, ncol=3)
	return(anno_plot)
}

## run enrich and plot
p <- enrichedGene(listEntrezID=listEntrezID)
ggsave('geneAnnotation.png', width = 16, height = 10, plot=p)
