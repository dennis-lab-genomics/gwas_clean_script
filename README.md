# GWAS Summary Statistics QC and Standardization Pipeline
  
## Pre-requisite Software:
* [Python 3.7](https://www.python.org/), and the following packages used:
* [numpy (*>*=1.18.1 recommded)](https://numpy.org/)
* [pandas (*>*=0.24.2 recommended)](http://pandas.pydata.org/)
* [matplotlib (*>*=3.1.3 recommended)](https://matplotlib.org/)
* [assocplots](https://github.com/khramts/assocplots)
* [manhattan_generator](https://github.com/pgxcentre/manhattan_generator)
* [pyliftover](https://pypi.org/project/pyliftover/)

## Running the pipeline
Sample execution definition file **args.json** provided as template.

### Local
```bash
$ python3 {path_to_args.json_file}
```

### HPC
```bash
$ python3 {path_to_args.json_file}
```
```bash
$ gwas_launch_script.sh {trait_name} {path_to_args.json_file} {approximate_size_of_summ_stat_file_in_Gb}
```
*Modify the paths in **gwas_launch_script.sh** to work on your system.

## Procedures performed by the pipeline:
* Record meta-data of the pipeline execution:
  * Genotyping array used in the source GWAS.
  * Max sample size in source GWAS.
  * Genome build used for SNP coordinates in the input summ stats.
  * Allele definition in input summ stats.
  * Version of dbSNP data used in the standardization step.
  * Analyst's name.
* Check summ stats for P-values smaller than Python can handle. If they exist, set to smallest float value in Python, and save original P-values to a separate file.
* Drop SNPs if:
  * P-value not float or not within (0,1).
  * SNP ID not valid rsid.
  * Chromosome value doesn't match human.
  * Chromosome not autosomal (unless pass argument to keep).
  * Base Pair not integer.
  * One of the alleles not SNV.
  * Not biallelic.
  * SNP is palindromic (unless effect allele frequency provided and SNP can be rescued).
* Standardization:
  * Standardize names of core columns.
  * Fill in missing rsid, chromosome, base pair information based on 1000 Genomes P3 VCF.
  * Attempt to rescue palindromic SNPs.
  * Report SNP effect with respect to the Alternative allele with the highest frequency in 1000 Genomes P3 (Primary GWAS population).
  * Perform Liftover to desired genome build, + sense strand.
* Output summ stat qc and description file:
  * How many SNPs were kept.
  * Breakdown of reasons for dropping SNPs with SNP counts.
  * Number of palindromic SNPs.
  * Mean SNP effect size before and after processing (for SNPs with P-value < 0.001).
  * Chromosome distribution.
  * Allele counts.
  * Allele combination table.
* Output summ stats: 
  * All dropped SNPs.
  * Clean SNPs with both the standardized columns and the original columns.
  * Clean SNPs with standardized columns only.
* Plots:
  * Manhattan Plot.
  * QQ Plot.

*All output file names follow format: ```{trait_code}_{genome_build}_{pmid}_{analyst_initials}_{file_type}_{date}.{extension}```

## Argument JSON file:
### General rules:
* If you don't have a value for a field, set it to ```false```.
* If parsing values (e.g. chromosome and bair pair) from a column that contains combined values, set the column name field for that value to ```false```.

### Field descriptions:
* **chr** - name of chromosome column.
* **bp** - name of base pair column.
* **rsid** - name of SNP id column. If too many rows don't have valid rsid, set to ```false```.
* **effect** - name of column containing allele with respect to which the effect estimate is reported.
* **other** - name of column containing the other allele.
* **beta** - name of beta column.
* **odds** - name of odds ratio column.
* **uci** - name of upper confidence interval column.
* **lci** - name of lower confidence interval column.
* **eaf** - name of column column containing frequency of the effect allele.
* **p** - name of pvalue column.
* **parse_cols** - this field contains instructions for parsing columns that are a combination of value types.
  * **required** -  name of such column.
  * **full** - set to ```true``` if column follows the ```{chrom}:{bp}_{A1}_{A2}``` pattern, otherwise ```false```.
  * **sep** - symbol used to separate value types.
    
    *only used if **full** is ```false```
  * **chr** - 1-based order of chromosome value as it appears in the column.
    
    *only used if **full** is ```false```
  * **bp** - 1-based order of base pair value as it appears in the column.
    
    *only used if **full** is ```false```
  * **rsid** - 1-based order of rsid value as it appears in the column.
    
    *only used if **full** is ```false```
  * **effect** - 1-based order of effect allele value as it appears in the column.
    
    *only used if **full** is ```false```
  * **other** - 1-based order of other allele value as it appears in the column.
    
    *only used if **full** is ```false```
* **bld_in** - genome build of input summ stats (e.g. HG19).
* **bld_out** - genome build of input summ stats (e.g. HG38).
* **initials** - initials of analyst.
* **trait** - trait code.
* **sep** - column separator of the summ stat file.
* **input_path** - path to the summ stat file. If running on HPC, this has to be defined in terms of the container bound directory names.
* **output_path** - path to directory where a trait name output directory will be created. If running on HPC, this has to be defined in terms of the container bound directory names.
* **harmonize** - set to ```false``` if you only want basic cleaning and plots. If ```true``` must contain **effect** and **other** columns, as well as an effect estimate.
* **palindromic** - threshold minor allele frequency for rescuing palindromic SNPs, or ```false``` to drop all palindromic SNPs.
* **keep_X** - ```true``` if you don't want X chromosome to be dropped from results.
* **keep_Y** - ```true``` if you don't want Y chromosome to be dropped from results.
* **keep_M** - ```true``` if you don't want MT chromosome to be dropped from results.
* **pop_code** - code of primary GWAS population as defined in 1000 Genomes (```[AFR, AMR, EAS, EUR, SAS]```).
* **overwrite** - ```true``` if you want to overwrite previous output for this trait.
* **array** - genotyping array used.
* **N** - max GWAS sample size.
* **pubid** - pubmed ID.
* **dbSNP** - dbSNP version of 1000 Genomes P3 VCF.

## Acknowledgements:
This script uses an R package created at Vanderbilt University. 
