# Mini-Project (COMP 483) v1.0

The script `miniproject.py` is a Python wrapper for examining and comparing the uropathogenic E. coli strains *HM27*, *HM46*, *HM65*, and *HM69* to find if differential gene expression drives uropathogenicity. Genotypes and gene expression are retrieved from NCBI databases within the script, so there is no user-defined input or sample data. Output includes annotations of genomes (`Angela_Andaleon/UPEC.log`) and gene expression of each strain sorted from highest to lowest (`Angela_Andaleon/STRAIN NAME_normalized_sorted.tsv`). 

This repository can be used with the commands `git clone https://github.com/aandaleon/COMP483_miniproject.git; cd COMP483_miniproject/; python3 miniproject.py`.

In detail, the script is a wrapper around various genome assembly, annotation, and read mapping and quantification tools for the aforementioned strains. It begins by retrieving assemblies from NCBI and reporting contigs and total read lengths, and annotates them with Prokka. It also writes discrepancies of CDS and tRNA in comparison to the *E. coli K-12* genome. It then retrieves transcriptomes from NCBI with the SRA Toolkit and builds the index, maps reads, and quantifies expression with Bowtie2, Tophat2, SAMtools, and Cufflinks, and sorts the gene expression output from highest to lowest gene expression.

## Requirements
* Python3 and libraries `os` and `pandas`
* [Prokka](https://www.ncbi.nlm.nih.gov/pubmed/24642063) (1.13.3)
* [SRA Toolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/) (2.8.2)
* [Bowtie](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334321/) (1.2.2)
* [TopHat](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2672628/) (2.1.1)
* [SAMtools](https://www.ncbi.nlm.nih.gov/pubmed/19505943) (1.7)
* [Cufflinks](https://www.ncbi.nlm.nih.gov/pubmed/20436464) (2.2.1)
