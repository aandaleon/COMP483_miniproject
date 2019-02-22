# Mini-Project (COMP 483)

The script `miniproject.py` is a Python wrapper for examining and comparing the uropathogenic E. coli strains *HM27*, *HM46*, *HM65*, and *HM69* to find if differential gene expression drives uropathogenicity. Genotypes and gene expression are retrieved from NCBI databases within the script, so there is no user-defined input. Output includes annotations of genomes (`Angela_Andaleon/UPEC.log`) and significant differences of transcript expression between strains (`Angela_Andaleon/TO BE NAMED`). 

This repository can simply be used with the commands `git clone https://github.com/aandaleon/COMP483_miniproject.git; python3 miniproject.py` or alternatively `wget https://raw.githubusercontent.com/aandaleon/COMP483_miniproject/master/miniproject.py`.

## Requirements
* Python dependency `os`
* [Prokka](https://www.ncbi.nlm.nih.gov/pubmed/24642063)
* [TopHat](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2672628/)
* [Cufflinks](https://www.ncbi.nlm.nih.gov/pubmed/20436464)
