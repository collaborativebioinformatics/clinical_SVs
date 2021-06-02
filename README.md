# clinical_SVs

Clinically reportable structural variant calls

LOGO

## Goals

Write a **workflow and/or app** to **annotate structural variants** (SVs) calls with **clinically relevant information**.

Eventually, the workflow could call the SVs from sequencing reads (e.g. from a BAM file).

## Overview Diagram

![](sv-clinic-workflow.jpg)

## Notes/Documentation

The new fields in the annotated VCF will be:

| name   | description                           |
|--------|---------------------------------------|
| AF     | Allele frequency                      |
| CLINSV | the ids of matching known clinical SV |
| GENE   | name of the gene(s) overlapped        |

## Installation

## Quick Start

## Component Details

Two ways to use our tools:

1. A command-line workflow that can do
   - SV discovery
   - Annotation
   - Automated QC graphs of supporting reads
1. An interactive app to annotate a VCF with SV calls and visualize the results.

### SV calling

If needed SVs can be called using parliament2. 
We will provide commands to run this variant discovery and integrate it to the workflow

### Annotation of SVs in R

The different modules of the annotation are written as functions and saved in separate files.
Then the master annotation script can read a VCF, *source* these functions and use them to annotate the SVs. 
See the current master annotation script [`annotate_vcf.R`](annotate_vcf.R) and the different scripts *source*d inside.

The `server.R` file of the app can also use the same approach: source the same functions and use them to annotate a VCF before displaying the info in the app.

#### Annotation modules

- [ ] [`annotate_genes.R`](annotate_genes.R) TODO. Currently returns dummy test values.
- [ ] [`annotate_frequency.R`](annotate_frequency.R) TODO. Currently returns dummy test values.
- [ ] [`annotate_known_clinical_SVs.R`](annotate_known_clinical_SVs.R) TODO. Currently returns dummy test values.
- [ ] ???

### Visualization of aligned reads around a SV

We could show the reads around a SV as a quality control if the BAM is available.
Either invent our own graph (for example in python), or run an existing tool like [samplot](https://github.com/ryanlayer/samplot).
