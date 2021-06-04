VCFanalyse.py
-------------
This script takes a VCF file as input. It analyzes the chromosome by dividing the chromosome into 'n' number of the equal-sized bin (parameter bins) and counts the different types of structural variants (SVs) to plot them as a figure. It allows the user to analyze what parts of the chromosomes are most affected by what type of variation.

This binning technique can also be used to create a patient SV profile to train a neural network model for patients for diagnosis and treatment (based on similarity with other successful treatments of the past) of patients

**Please note that I have added a VCF parser and data processing functions in the python file(s). In the future, it might be helpful and help other people.**

The parameters for the program are as following
-chr Chromosome ID (eg.. chr1, chr2)
-vcf VCF filename (and location if not in the same folder)
-bins As described above.
-minlen Cutoff length for the SV (Please note that SVTYPE BND does not have length defined)

Example command:\
`python VCFanalyse.py -vcf NA19461.final.manta.diploidSV.vcf -chr chr2 -bins 100 -minlen 100`

Example output for Chromosome 2 for the example command above:\
<img src="https://github.com/collaborativebioinformatics/clinical_SVs/blob/main/VCF_parser/output.png">