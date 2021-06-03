# clinical_SVs

Clinically reportable structural variant calls


![](ClinicalSVsLogo.png)


## Contributors

Jean Monlong  - `Lead, Liaison`

Rupesh Kesharwani - `Sysadmin and code developer`

Pranav Khade - `Data Guru and code Developer`

Weiyu Zhou - `App Developer`

Ahmad Al Khleifat - `Writer`



## Goals

Write a **workflow and/or app** to **annotate structural variants** (SVs) calls with **clinically relevant information**.

Eventually, the workflow could call the SVs from sequencing reads (e.g. from a BAM file) and extract some QC information from raw reads.

## Overview Diagram

![](sv-clinic-workflow.jpg)

## Notes/Documentation

For the *gene-level* metrics, the table contains the following columns.

| name       | description                                                            |
|------------|------------------------------------------------------------------------|
| gene       | names of genes overlapped, separated by `\\|`                          |
| variant_id | SV ID                                                                  |
| chr        | chromosome name                                                        |
| start      | start position                                                         |
| end        | end position                                                           |
| size       | size of the SV in bp                                                   |
| frequency  | allele frequency                                                       |
| svtype     | type of SV. E.g. DEL, DUP, INS, ...                                    |
| clinsv     | dbVar accession IDs of matching known clinical SVs (separated by `\\|` |

See [`clinical-sv-table.csv`](R/clinical-sv-table.csv) for an example on our test data (HG002 from GIAB).

We will also run a gene set enrichment and highlight SVs in enriched pathways/diseases.
This is will represent a set of *patient-level* metrics.
The output will be a set of graphs (image files) 

## Installation

## Quick Start

## Component Details

Two ways to use our tools:

1. A **command-line workflow** that can do
   - SV discovery
   - Annotation
   - Automated QC graphs of supporting reads
1. An **interactive app** to annotate a VCF with SV calls and visualize the results.
   - Building on GeneVar but now the user can upload their own VCF

### SV calling

If needed SVs can be called using parliament2. 
We will provide commands to run this variant discovery and integrate it to the workflow

### Annotation of SVs in R

The different modules of the annotation are written as functions and saved in separate files.
Then the master annotation script can read a VCF, *source* these functions and use them to annotate the SVs. 
See the current master annotation script [`annotate_vcf.R`](R/annotate_vcf.R) and the different scripts *source*d inside.

The `server.R` file of the app can also use the same approach: source the same functions and use them to annotate a VCF before displaying the info in the app.

#### Annotation modules

The resources used in the modules are downloaded and prepared by the [`prepare_annotation_data.R`](R/prepare_annotation_data.R). 

- [X] [`annotate_genes.R`](R/annotate_genes.R) 
   - Gencode v35
   - Filter: keeps only SVs overlapping CDS regions.
   - Uses a "wide" format: new field *GENE* lists all genes overlapped separated by `|`.
- [X] [`annotate_frequency.R`](R/annotate_frequency.R) 
   - Uses the position and the *SVLEN*, *SVTYPE* fields to match SVs in the gnomAD-SV catalog (i.e. make sure the SV calls contains these fields).
   - Lenient matching criteria: 10% minimum reciprocal overlap; insertions can be as distant as 100 bp.
   - New field *AF* reports the maximum frequency across all the SVs matched in gnomAD-SV.
- [x] [`annotate_known_clinical_SVs.R`](R/annotate_known_clinical_SVs.R)
   - Known clinical SVs: dbVar *nstd102* study
- [ ] [`geneFunctionalAnnotation.R`](R/geneFunctionalAnnotation.R)
   - Need a list of list of genes similar to `demo/listofENSEMBLID.txt` (Other than ENTREZ Gene ID; such as SYMBOL / REFSEQ / ENSEMBL) to output a three types of Disease Ontology plots such as `(demo/geneAnnotation.png)` including a barplot (high level catogory), a dot plot (show upto 20 diseases association) and a disease-gene network graph.
   - The list of genes can be extracted from annotated vcf based on any SV types.
- [ ] Known SVs in cancer from COSMIC?

### Visualization of aligned reads around a SV

We could show the reads around a SV as a quality control if the BAM is available.
Either invent our own graph (for example in python), or run an existing tool like [samplot](https://github.com/ryanlayer/samplot).
