# gwas_clean_script

## Script used for cleaning and quality checking GWAS Summary Statistics
## Pre-requisite Software:
### Procedures performed in the script:


### Using the script
Use the shell command in the **sample_script_call.txt** file as a template. All the arguments are required except the arguments that act as booleans.
If a dataset has features that the script does not cover, the script will have to be modified manually to accomodate for the edge case.

The **cleaning for small p-values** and **plotting manhattan** steps are very time consuming as they rely on external software. It is normal for these steps to take 10+ minutes depending on a standard personal machine.

## Acknowledgements:
