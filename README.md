# gwas_clean_script

## Script used for cleaning and quality checking GWAS Summary Statistics

## Pre-requisite Software:
* [Python 3.7](https://www.python.org/), and the following packages used:
* [numpy (*>*=1.18.1 recommded)](https://numpy.org/)
* [pandas (*>*=0.24.2 recommended)](http://pandas.pydata.org/)
* [matplotlib (*>*=3.1.3 recommended)](https://matplotlib.org/)
* [assocplots](https://github.com/khramts/assocplots)
* [manhattan_generator](https://github.com/pgxcentre/manhattan_generator)

## Procedures performed in the script:
* Check data file for P-values smaller than Python can handle. If they exist, change them to the smallest float value in Python, and save original P-values to a separate file.
* Remove P-values that are zero or non-numeric.
* Parse and create chromosome and position columns from SNP if these columns don't exist in the data file.
* Remove data that describe non-autosomal chromosomes.
* Remove empty columns or rows that are missing Chromosome, Base Pair, SNP or P values.
* Determine SNP distribution based on Chromosomes.
* Produce Allele combination table.
* Determine Biallelic SNP distribution based on Nucleotides.
* Determine Non-Biallelic SNP distribution based on Type.
* Remove Non-Biallelic SNPs.
* Create Manhattan Plot.
* Create QQ Plot.

## Command Arguments:
* **--initials** - your initials for output naming.
* **--trait** - trait name.
* **--chromosome** - if chromosome column exists, its name.
* **--position** - if chromosome column exists, its name.
* **--pval** - P-value column name.
* **--snp** - SNP column name.
* **--effect** - reference allele column name.
* **--noneffect** - alternative allele column name.
* **--sep** - the separation strategy of the file must be **enclosed in quotation marks and without a slash**.
* **--has_chr_and_pos** - has chromosome and position (Base Pair) column.
* **--slash_required** - e.g. if separation strategy is '\t', enter 't' for sep and set this flag.
* **--remove_cleaned_datafiles** - delete the final cleaned version of the data file after the script is finished.
* **--input_path** - path to the input data.
* **--output_path** - recommended to use output folder.
* **--r_script_path** - path to the R script provided for cleaning small P-values.
* **--manhattan_path** - path to the *manhattan_generator* package downloaded from github.

## Using the script
Use the shell command in the **sample_script_call.txt** file as a template.
If a dataset has features that the script does not cover, the script will have to be modified manually to accomodate for the edge case.

The **cleaning for small p-values** and **plotting manhattan** steps are very time consuming as they rely on external software. It is normal for these steps to take 10+ minutes depending on a standard personal machine.

### Note
* If chromosome and base pair position columns aren't included in the dataset, the script makes assumptions about the structure of the SNP columns.
* Script assumes that dataset has header.

## Acknowledgements:
This script uses an R package created at Vanderbilt University. 