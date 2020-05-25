
# Provide sep as only the letter, don't use backslash
# back_slash_require flag set
# Script prepends an additional / to the sep argument
# If chromosome and base pair position columns aren't included in the dataset, the script makes assumptions about the structure of the SNP columns
# Script assumes that dataset has header

import os
import re
import sys
import math
import argparse
import subprocess
from datetime import datetime
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import assocplots.qqplot as qqplot


indexes = []
index = 0

def get_non_autosomal(chrom):
    global index
    if (type(chrom) == int or (type(chrom) == str and chrom.isdigit())) and 0 < int(chrom) < 23:
        index += 1
        return chrom
    else:
        indexes.append(index)
        index += 1
        return np.nan
        

def get_small_p(p):
    global index
    if p == float(0):
        indexes.append(index)
        index += 1
        return np.nan
    else:
        index += 1
        return p
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--initials", type=str)
    parser.add_argument("--trait", type=str)
    parser.add_argument("--chromosome", default="chr", type=str)
    parser.add_argument("--position", default="bp", type=str)
    parser.add_argument("--pval", type=str)
    parser.add_argument("--snp", type=str)
    parser.add_argument("--effect", type=str)
    parser.add_argument("--noneffect", type=str)
    parser.add_argument("--sep", type=str)
    parser.add_argument("--has_chr_and_pos", action="store_true")
    parser.add_argument("--back_slash_required", action="store_true")
    parser.add_argument("--remove_cleaned_datafiles", action="store_true")
    parser.add_argument("--input_path", type=str)
    parser.add_argument("--output_path", type=str)
    parser.add_argument("--r_script_path", type=str)
    parser.add_argument("--manhattan_path", type=str)

    args = parser.parse_args()

    initials = args.initials
    trait = args.trait
    chromosome = args.chromosome
    position = args.position
    p_value = args.pval
    snp = args.snp
    allele1 = args.effect
    allele2 = args.noneffect
    separation_strategy_letter = args.sep
    has_chr_and_pos = args.has_chr_and_pos
    back_slash_required = args.back_slash_required
    remove_cleaned_datafiles = args.remove_cleaned_datafiles
    data_file_path = args.input_path
    output_dir = args.output_path
    r_script_path = args.r_script_path
    manhattan_generator_path = args.manhattan_path

    date = ''.join(str(datetime.date(datetime.now())).split('-'))
    separation_strategy = f"\\{separation_strategy_letter}" if back_slash_required else separation_strategy_letter

    cleaned_data_file_path = f'{output_dir}/{trait}_small_p_clean_data_{initials}_{date}.txt'
    bad_pvalue_file_path = f'{output_dir}/{trait}_bad_pvalues_{initials}_{date}.txt'
    description_file_path = f'{output_dir}/{trait}_qc_checks_{initials}_{date}.txt'
    modified_data_path = f'{output_dir}/clean_{trait}_data.txt'
    manhattan_output_path = f'{output_dir}/{trait}_manhattan_{initials}_{date}'
    qq_plot_path = f'{output_dir}/{trait}_qqplot_{initials}_{date}.jpg'


    print("Importing dataset...")
    df = pd.read_csv(data_file_path, sep=separation_strategy)

    try:
        os.mkdir(output_dir)
    except OSError:
        print(f"Couldn't create output directory at: {output_dir}") 
    else:
        print(f"Successfully created the directory {output_dir}")


    print("Extracting small pvalues...")
    p_column_number = df.columns.get_loc(p_value) + 1
    # Allow for multiple columns to be used in script?
    os.system(f"Rscript '{r_script_path}' '{data_file_path}' '{cleaned_data_file_path}' $'{separation_strategy}' True {p_column_number}")


    print("Importing new file with small pvalues...")
    df = pd.read_csv(cleaned_data_file_path, sep=separation_strategy)

    chromosome = chromosome.replace("-", ".")
    position = position.replace("-", ".")
    p_value = p_value.replace("-", ".")
    snp = snp.replace("-", ".")
    allele1 = allele1.replace("-", ".")
    allele2 = allele2.replace("-", ".")

    required = [chromosome, snp, position, p_value]


    print("Cleaning bad pvalues...")
    file_bad_pvalues = open(bad_pvalue_file_path, 'a')

    index = 0
    df[p_value].apply(get_small_p)

    for i in indexes:
        row = df.iloc[i]
        file_bad_pvalues.write(f'Index {i} : {row[snp]} : {row[p_value]}\n')
        
    if len(indexes) == 0:
        file_bad_pvalues.write("No bad pvalues to clean.")

    del indexes[:]
    index = 0
    df[p_value] = df[p_value].apply(get_small_p)

    file_bad_pvalues.close()
            

    print("Cleaning and describing dataset...")
    file = open(description_file_path,'a')

    # Extracting chromosome and position data from SNP columns
    if not has_chr_and_pos:
        df[chromosome] = df[snp].apply(lambda x: int(x.split(':')[0][3:]) if len(x.split(':')) > 1 else np.nan)
        df[position] = df[snp].apply(lambda x: int(x.split(':')[1]) if len(x.split(':')) > 1 else np.nan)

    # Finding non-autosomal chromosomes
    index = 0
    df[chromosome].apply(get_non_autosomal)

    chromosome_number = len(indexes)
    if chromosome_number != 0:
        file.write('The following are the SNPs from non-autosomal chromosomes that were removed: \n')
        for i in indexes:
            row = df.iloc[i]
            file.write(f'{row[snp]} : {row[chromosome]}\n')
    else:
        file.write("No non-autosomal chromosomes to report.\n")

    del indexes[:]
    index = 0
    df[chromosome] = df[chromosome].apply(get_non_autosomal)

    # Cleaning bad data
    before = df.shape
    df.dropna(axis=0, subset=required, inplace=True)
    df.dropna(axis=1, how='all', inplace=True)
    after = df.shape
    file.write(f'\n{before[0] - after[0] - chromosome_number} rows removed beacuse they had nan values.')
    file.write(f'\n{before[1] - after[1]} columns removed because they were empty. \n\n')

    # Determine distribution of data in chromosomes and Alleles
    # Write results to file

    counts = df[chromosome].value_counts()
    chroms = np.sort(list(map(int, df[chromosome].unique())))

    biallelic = df[(df[allele1].isin(['a','c','g','t', 'A', 'G', 'C', 'T'])) & (df[allele2].isin(['a','c','g','t', 'A', 'G', 'C', 'T']))]
    biallelic_A1 = df[(df[allele1].isin(['a','c','g','t', 'A', 'G', 'C', 'T']))][allele1]
    biallelic_A2 = df[(df[allele2].isin(['a','c','g','t', 'A', 'G', 'C', 'T']))][allele2]

    non_biallelic1 = df[allele1].apply(lambda x: 'multi' if re.match("^[acgtACGT][acgtACGT]+", x) else (x if x not in ['a','c','g','t', 'A', 'G', 'C', 'T'] else 'biallelic'))
    non_biallelic2 = df[allele2].apply(lambda x: 'multi' if re.match("^[acgtACGT][acgtACGT]+", x) else (x if x not in ['a','c','g','t', 'A', 'G', 'C', 'T'] else 'biallelic'))

    length = len(df[chromosome])

    # Writing chromosome distribution to text file
    for i in chroms:
        file.write(f'\nCHR{i}: {counts.loc[i]}; {counts.loc[i]/length}%')

    # Writing Allele Combination Table
    file.write('\n\n\nBiallelic Combination Table:\n')
    for allele in biallelic[allele1].unique():
        df_t = biallelic[biallelic[allele1]==allele]
        file.write(f'\n  {allele}:\n')
        file.write(f'{df_t[allele2].value_counts().to_string()}')
        
    # Writing SNP distribution to text file
    file.write('\n\n\nA1 biallelic:\n')
    file.write(f'{biallelic_A1.value_counts().to_string()}\n')
    file.write('\nA1 SNP types (non-biallelic SNPs will be removed):\n')
    file.write(f'{non_biallelic1.value_counts().to_string()}\n')
        
    file.write('\n\nA2 biallelic:\n')
    file.write(f'{biallelic_A2.value_counts().to_string()}\n')
    file.write('\nA2 SNP types (non-biallelic SNPs will be removed):\n')
    file.write(f'{non_biallelic2.value_counts().to_string()}\n')

    df[allele1] = df[allele1].apply(lambda x: x if x in ['a','c','g','t', 'A', 'G', 'C', 'T'] else np.nan)
    df[allele2] = df[allele2].apply(lambda x: x if x in ['a','c','g','t', 'A', 'G', 'C', 'T'] else np.nan)
    df.dropna(axis=0, subset=[allele1, allele2], inplace=True)
        
    file.close()

 
    # Perform Bonferronni correction
    n = df[p_value].shape[0]
    bonferronni = -math.log10(0.05/n)


    print("Exporting clean data for plotting...")
    df.to_csv(modified_data_path, sep='\t')


    print("Creating manhattan plot...")
    os.system(f"python3 {manhattan_generator_path} --twopoint {modified_data_path} --col-chr {chromosome} --col-name {snp} --col-pos {position} --col-pvalue {p_value} --bp --use-pvalues --abline {bonferronni} --significant-threshold {bonferronni} --no-annotation --significant-point-size 2 --point-size 1 --graph-title '{trait} GWAS Manhattan Plot' --chr-text-size 10 --output {manhattan_output_path}")


    print("Creating QQ plot...")
    # mpl.rcParams['figure.dpi']=150
    # mpl.rcParams['savefig.dpi']=150
    # mpl.rcParams['figure.figsize']=7.375, 3.375

    mpl.rcParams['figure.dpi']=100
    mpl.rcParams['savefig.dpi']=100
    mpl.rcParams['figure.figsize']=5.375, 5.375

    qqplot.qqplot([df[p_value]], 
        [f'{trait} GWAS'], 
        color=['b'], 
        fill_dens=[0.2], 
        error_type='theoretical', 
        distribution='beta',
        title=f'{trait} GWAS QQ Plot')

    mpl.pyplot.legend(loc=0)
    plt.savefig(qq_plot_path, dpi=300)

    
    # Wrapping up
    os.remove(cleaned_data_file_path)
    if remove_cleaned_datafiles:
        os.remove(modified_data_path)
