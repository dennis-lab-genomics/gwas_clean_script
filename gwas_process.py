import os
import shutil
import re
import sys
import math
import json
import argparse
import subprocess
import multiprocessing 
import asyncio
from datetime import datetime
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import assocplots.qqplot as qqplot
from pyliftover import LiftOver

pd.options.mode.chained_assignment = None


G_BUILDS = SimpleNamespace(**{
    "hg15" : "HG15",
    "hg16" : "HG16",
    "hg17" : "HG17",
    "hg18" : "HG18",
    "hg19" : "HG19",
    "hg38" : "HG38"
})


PATHS = SimpleNamespace(**{
    "data" : "/work",
    "chain_files" : "/work/shared/data/liftover_chains",
    "vcf_files" : "/work/shared/data/vcf_variant_data",
    "vcf" : "1000GENOMES-phase_3.vcf.gz",
    "p_val_script_path" : "/app/soft/clean_file_for_R_small_values.R",
    "vcftools" : "/app/soft/vcftools",
    "manhattan_script_path" : "/app/soft/manhattan_generator.py"
})


COLS = SimpleNamespace(**{
    "chr" : "chr_JD_STANDARD",
    "bp" : "bp_JD_STANDARD",
    "rsid" : "rsid_JD_STANDARD",
    "effect" : "effect_JD_STANDARD",
    "other" : "other_JD_STANDARD",
    "beta" : "beta_JD_STANDARD",
    "odds" : "or_JD_STANDARD",
    "uci" : "uci_JD_STANDARD",
    "lci" : "lci_JD_STANDARD",
    "eaf" : "eaf_JD_STANDARD",
    "p" : "p_JD_STANDARD",
    "full_coords" : "full_coords",
    "index" : "temp_numeric_rsid",
    "drop" : "drop"
})


COLS_OUT = SimpleNamespace(**{
    "chr" : "chr",
    "bp" : "bp",
    "rsid" : "rsid",
    "effect" : "effect",
    "other" : "other",
    "p" : "p"
})


CLEAN_STEPS = SimpleNamespace(**{
    "chr" : "\t\tbad chr:",
    "rsid" : "\t\tbad rsid:",
    "bp" : "\t\tbad bp:",
    "p" : "\t\tbad p-values:",
    "pre_clean" : "\n\tInitial cleaning:",
    "lift1" : "\t\tfirst liftover:",
    "lift2" : "\n\t\tsecond liftover:",
    "update" : "\t\tharmonizing non-palindromic SNPs:",
    "drop_palindr" : "\n\t\tdropped palindromic SNPs:",
    "save_palindr" : "\n\t\tambiguous palindromic SNPs:",
    "total_palindromic" : "\n\nTotal number of palindromic SNPs:",
    "standardize" : "\n\tHarmonizing SNPs to Ensembl:",
    "chrX" : "\n\tCleaning non-autosomal (chrX):",
    "chrY" : "\n\tCleaning non-autosomal (chrY):",
    "chrM" : "\n\tCleaning non-autosomal (chrM):",
    "indel" : "\t\t\tSNPs have indels:",
    "non_biallelic" : "\t\t\tSNPs are not biallelic:",
    "bad_allele" : "\t\t\tbad allele:",
    "allele_total" : "\t\tfiltering alleles:"
})
    
    
COMPLEMENT = {
    "A" : "T",
    "T" : "A",
    "G" : "C",
    "C" : "G"
}


CHROM_CODE = {
    "1" : 1,
    "2" : 2,
    "3" : 3,
    "4" : 4,
    "5" : 5,
    "6" : 6,
    "7" : 7,
    "8" : 8,
    "9" : 9,
    "10" : 10,
    "11" : 11,
    "12" : 12,
    "13" : 13,
    "14" : 14,
    "15" : 15,
    "16" : 16,
    "17" : 17,
    "18" : 18,
    "19" : 19,
    "20" : 20,
    "21" : 21,
    "22" : 22,
    "X" : 23,
    "x" : 23,
    "Y" : 24,
    "y" : 24,
    "M" : 25,
    "MT" : 25,
    "m" : 25,
    "mt" : 25
}


CHROM_VAL = {
    1 : "1",
    2 : "2",
    3 : "3",
    4 : "4",
    5 : "5",
    6 : "6",
    7 : "7",
    8 : "8",
    9 : "9",
    10 : "10",
    11 : "11",
    12 : "12",
    13 : "13",
    14 : "14",
    15 : "15",
    16 : "16",
    17 : "17",
    18 : "18",
    19 : "19",
    20 : "20",
    21 : "21",
    22 : "22",
    23 : "X",
    24 : "Y",
    25 : "M"
}


STRANDS = SimpleNamespace(**{
    "pos" : 1,
    "neg" : 2,
    "mix" : 3
})


VCF_DTYPES = {
    "#CHROM" : str,
    "POS" : int,
    "ID" : str,
    "REF" : str,
    "ALT" : str,
    "QUAL" : str,
    "FILTER" : str,
    "INFO" : str
}
     
    
ARGS_TO_CHECK = ["effect", "other", "p", "rsid", "chr", "bp", "beta", "odds", "uci", "lci", "eaf"]


REQUIRED = [COLS_OUT.rsid, COLS_OUT.p, COLS_OUT.effect, COLS_OUT.other]


def check_required_cols(args):
    for col in REQUIRED:
        if not vars(args)[col]:
            if col == COLS_OUT.rsid and args.parse_cols.rsid:
                continue
            if col == COLS_OUT.chr and (args.parse_cols.full or args.parse_cols.chr):
                continue
            if col == COLS_OUT.bp and (args.parse_cols.full or args.parse_cols.bp):
                continue
            if col == COLS_OUT.effect and (args.parse_cols.full or args.parse_cols.effect):
                continue
            if col == COLS_OUT.other and (args.parse_cols.full or args.parse_cols.other):
                continue
            else:
                print(f"All of the following columns are required:\n{REQUIRED}")
                sys.exit()
            
    
def create_dir(path, overwrite):
    if os.path.exists(path):
        if overwrite:
            if os.path.isdir(path):
                sys_call(shutil.rmtree, path, "")
            else:
                sys_call(os.remove, path, "")
            sys_call(os.mkdir, path, f"couldn't create directory at: {path}")

    else:
        sys_call(os.mkdir, path, f"couldn't create directory at: {path}")
    
    
def check_column_names(cols, args):
    for arg in ARGS_TO_CHECK:
        if args[arg] and args[arg] not in cols:
            print(f"Input column names incorrect (check for spaces). Must be found in the following list: \n\t{list(cols)}")
            sys.exit()
  
 
def is_float(val):
    try:
        float(val)
        return True
    except ValueError:
        return False
      
         
def is_string(val):
    try:
        str(val)
        return True
    except ValueError:
        return False
    
    
def round_to_3(num):
    return round(num, -int(np.floor(np.log10(abs(num))))+2)


def write_drop_breakdown(file, tup, total):
    (snps, text) = tup        
    line = f"{text} {snps} {f'({round_to_3(snps/total*100)}%)' if snps > 0 else ''}\n"
    file.write(line)
        
        
def classify_snp(a):
    try:
        a = str(a).upper()

        if a in COMPLEMENT.keys():
            return a
        elif "," in a:
            return "non_biallelic"
        elif not a in COMPLEMENT.keys():
            return "indel"

    except ValueError:
        return "bad"
    
    
def sys_call(f, arg, msg):
    try:
        f(arg)
    except OSError:
        print(msg)
        sys.exit()

        
def get_clean_vcf(chrom, out_dir, filt, filter_path):
    cmd = f"{PATHS.vcftools} \
    --gzvcf {os.path.join(PATHS.vcf_files, PATHS.vcf)} \
    {filt} {filter_path} \
    --recode \
    --recode-INFO-all \
    --out {os.path.join(out_dir, f'ensembl_snps_{chrom}')}"
    sys_call(os.system, cmd, "vcftools failed")

    vcf_path = os.path.join(out_dir, f"ensembl_snps_{chrom}.recode.vcf")
    
    bad_lines = 0
    with open(vcf_path) as file:
        for line in file.readlines(10**6):
            if line.startswith("##"):
                bad_lines += 1
            else:
                break

    vcf = pd.read_csv(vcf_path, sep = "\t", skiprows = list(range(bad_lines)), index_col = False, dtype=VCF_DTYPES)
    vcf.rename(columns={"#CHROM":"CHROM"}, inplace = True)
    vcf["CHROM"] = vcf.CHROM.astype(str)
        
    sys_call(os.remove, vcf_path, "Couldn't delete file")
    sys_call(os.remove, filter_path, "Couldn't delete file")
    
    return vcf


def infer_strand(effect, other, vcf_alt, vcf_ref):
    alleles = vcf_alt.split(",")
    alleles += vcf_ref
    
    if effect in alleles and other in alleles:
        return STRANDS.pos
    
    return STRANDS.neg


def harmonise_snp(effect, other, eaf, beta, odds, uci, lci, i, vcf_alt, strand):
    if strand == STRANDS.neg:
        effect[i] = COMPLEMENT[effect[i]]
        other[i] = COMPLEMENT[other[i]]
        
    if effect[i] not in vcf_alt.split(","):
        temp_allele = effect[i]
        effect[i] = other[i]
        other[i] = temp_allele
        
        try:
            if eaf:
                eaf[i] = 1 - float(eaf[i])
        except ValueError:
            pass
            
        try:
            if beta:
                beta[i] = float(beta[i] * -1)
        except ValueError:
            pass
        
        try:
            if odds:
                odds[i] = float(odds[i] ** -1)
        except ValueError:
            pass
        
        try:
            if uci and lci:
                temp_uci = float(lci[i]) ** -1
                lci[i] = float(uci[i]) ** -1
                uci[i] = temp_uci
        except ValueError:
            pass


def update_consensus_strand(strand, consensus):
    if consensus < 0:
        return strand
    elif strand != consensus:
        return STRANDS.mix
    else:
        return strand


def update_snps(chrom, gwas, vcf, has_id, harmonize):
    consensus_strand = -1
    palindromic_count = 0
    
    if vcf.shape[0] == 0 and gwas.shape[0] > 0:
        gwas[COLS.drop] = True
        return consensus_strand, palindromic_count
        
    if has_id:
        gwas[COLS.bp] = 0
        gwas[COLS.chr] = 0
    
    else:
        gwas[COLS.rsid] = ""
        
    if not has_id:
        
        if chrom != "all":
            gwas[COLS.chr] = CHROM_CODE[chrom]
        else:
            gwas[COLS.chr] = gwas[COLS.chr].apply(lambda x: CHROM_CODE[x])  
            
        bp = gwas[COLS.bp].array
        gwas_chr = gwas[COLS.chr].array
        drop = gwas[COLS.drop].array
        rsid = gwas[COLS.rsid].array
        effect = gwas[COLS.effect].array if COLS.effect in gwas.columns else False
        other = gwas[COLS.other].array if COLS.other in gwas.columns else False
        eaf = gwas[COLS.eaf].array if COLS.eaf in gwas.columns else False
        beta = gwas[COLS.beta].array if COLS.beta in gwas.columns else False
        odds = gwas[COLS.odds].array if COLS.odds in gwas.columns else False
        uci = gwas[COLS.uci].array if COLS.uci in gwas.columns else False
        lci = gwas[COLS.lci].array if COLS.lci in gwas.columns else False

        vcf_chr = vcf.columns.get_loc("CHROM")
        pos = vcf.columns.get_loc("POS")
        snp_id = vcf.columns.get_loc("ID")
        alt = vcf.columns.get_loc("ALT")
        ref = vcf.columns.get_loc("REF")
        vcf = vcf.to_numpy()  
        
        j = 0
        for i in range(bp.size):
            while j < vcf.shape[0] - 1 and bp[i] > vcf[j, pos]:
                j += 1
            
            if bp[i] == vcf[j,pos]:
                rsid[i] = vcf[j, snp_id]
                
                if harmonize and effect[i] != COMPLEMENT[other[i]]:
                    strand = infer_strand(effect[i], other[i], vcf[j,alt], vcf[j,ref])
                    harmonise_snp(effect, other, eaf, beta, odds, uci, lci, i, vcf[j,alt], strand)
                    consensus_strand = update_consensus_strand(strand, consensus_strand)
                elif harmonize:
                    palindromic_count += 1
                    
            else:
                drop[i] = True
    
    else:
        gwas[COLS.index] = gwas[COLS.rsid].apply(lambda x: int(x[2:]))
        gwas.sort_values([COLS.index], inplace = True)
        gwas_index = gwas[COLS.index].array
        
        vcf[COLS.index] = vcf.ID.apply(lambda x: int(x[2:]))
        vcf.sort_values([COLS.index], inplace = True)
        
        bp = gwas[COLS.bp].array
        gwas_chr = gwas[COLS.chr].array
        drop = gwas[COLS.drop].array
        rsid = gwas[COLS.rsid].array
        effect = gwas[COLS.effect].array if COLS.effect in gwas.columns else False
        other = gwas[COLS.other].array if COLS.other in gwas.columns else False
        eaf = gwas[COLS.eaf].array if COLS.eaf in gwas.columns else False
        beta = gwas[COLS.beta].array if COLS.beta in gwas.columns else False
        odds = gwas[COLS.odds].array if COLS.odds in gwas.columns else False
        uci = gwas[COLS.uci].array if COLS.uci in gwas.columns else False
        lci = gwas[COLS.lci].array if COLS.lci in gwas.columns else False
        
        vcf_index = vcf.columns.get_loc(COLS.index)        
        vcf_chr = vcf.columns.get_loc("CHROM")
        pos = vcf.columns.get_loc("POS")
        snp_id = vcf.columns.get_loc("ID")
        alt = vcf.columns.get_loc("ALT")
        ref = vcf.columns.get_loc("REF")
        vcf = vcf.to_numpy()
        
        j = 0 
        for i in range(gwas_index.size):
            while j < vcf.shape[0] - 1 and gwas_index[i] > vcf[j,vcf_index]:
                j += 1
                
            
            if gwas_index[i] == vcf[j,vcf_index]:
                gwas_chr[i] = CHROM_CODE[vcf[j,vcf_chr]]
                bp[i] = int(vcf[j,pos])
                
                if harmonize and effect[i] != COMPLEMENT[other[i]]:
                    strand = infer_strand(effect[i], other[i], vcf[j,alt], vcf[j,ref])
                    harmonise_snp(effect, other, eaf, beta, odds, uci, lci, i, vcf[j,alt], strand)
                    consensus_strand = update_consensus_strand(strand, consensus_strand)
                elif harmonize:
                    palindromic_count += 1
                    
            else:
                drop[i] = True
                
        gwas.drop(COLS.index, axis = 1, inplace = True)
    
    return consensus_strand, palindromic_count
    
    
def extract_pop_aaf(info, pop_code):
    items = info.split(";")
    filt = f"{pop_code}="
    aaf = list(filter(lambda x: x.startswith(filt), items))
    if len(aaf) < 1:
        return -1
    else:
        freq = aaf[0].split("=")[1].split(",")[0]
        return float(freq)
    
    
def get_maf(freq):
    if freq <= 0.5:
        return freq
    else:
        return 1 - freq
    
    
def infer_strand_palindromic(eaf, info, palindromic, pop_code):
    try:
        aaf = extract_pop_aaf(info, pop_code)
        if aaf < 0:
            return STRANDS.mix
        maf_eaf = get_maf(float(eaf))
        maf_aaf = get_maf(aaf)
        if maf_eaf < palindromic and maf_aaf < palindromic:
            if (float(eaf) < 0.5 and aaf < 0.5) or (float(eaf) >= 0.5 and aaf >= 0.5):
                return STRANDS.pos
            else:
                return STRANDS.neg
        else:
            return STRANDS.mix
    except ValueError:
        return STRANDS.mix
    
    
def save_palindromic_snps(gwas, vcf, consensus, palindromic, pop_code):
    gwas.sort_values([COLS.chr, COLS.bp], inplace = True)
    vcf.sort_values(["CHROM", "POS"], inplace = True)
    
    strand = consensus
    j = 0
    
    bp = gwas[COLS.bp].array
    chrom = gwas[COLS.chr].array
    drop = gwas[COLS.drop].array
    effect = gwas[COLS.effect].array
    other = gwas[COLS.other].array
    eaf = gwas[COLS.eaf].array if COLS.eaf in gwas.columns else False
    beta = gwas[COLS.beta].array if COLS.beta in gwas.columns else False
    odds = gwas[COLS.odds].array if COLS.odds in gwas.columns else False
    uci = gwas[COLS.uci].array if COLS.uci in gwas.columns else False
    lci = gwas[COLS.lci].array if COLS.lci in gwas.columns else False

    vcf_chr = vcf.columns.get_loc("CHROM")
    pos = vcf.columns.get_loc("POS")
    info = vcf.columns.get_loc("INFO")
    alt = vcf.columns.get_loc("ALT")
    vcf = vcf.to_numpy()  
    
    for i in range(bp.size):
        while j < vcf.shape[0] - 1 and bp[i] > vcf[j, pos]:
            j += 1

        if bp[i] == vcf[j,pos] and chrom[i] == CHROM_CODE[vcf[j,vcf_chr]]:
            if effect[i] == COMPLEMENT[other[i]]:
                if consensus == STRANDS.mix:
                    strand = infer_strand_palindromic(eaf[i], vcf[j,info], palindromic, pop_code)
                    if strand == STRANDS.mix:
                        drop[i] = True
                        continue

                harmonise_snp(effect, other, eaf, beta, odds, uci, lci, i, vcf[j,alt], strand)
        

def map_to_ensembl(chrom, gwas, out_dir, from_id, palindromic, pop_code, harmonize):
    if not from_id:
        gwas.sort_values([COLS.bp], inplace = True)
        
    filter_path = os.path.join(out_dir, f"filter_list_{chrom}.temp")
    filt = "--snps" if from_id else "--positions"
    cols = [COLS.rsid] if from_id else [COLS.chr, COLS.bp]
    
    print(f"\t\tGetting vcf for chrom {chrom}...")
    
    gwas[cols].to_csv(filter_path, sep = "\t", index = False, header = False)
    
    vcf = get_clean_vcf(chrom, out_dir, filt, filter_path)
    
    print(f"\t\tUpdatings SNPs for chrom {chrom}...")
    
    drop_map = SimpleNamespace(**{})
    
    (strand, palindromic_count) = update_snps(chrom, gwas, vcf, from_id, harmonize)
    gwas_bad = gwas[gwas[COLS.drop]]
    gwas = gwas[~gwas[COLS.drop]]
    
    drop_map.update = gwas_bad.shape[0]
    drop_map.palindromic_total = palindromic_count
    
    reason = CLEAN_STEPS.drop_palindr
    
    if harmonize:
        if not palindromic or (strand == STRANDS.mix and (COLS.eaf not in gwas.columns or not pop_code)):
            print(f"\t\tDropping palindromic SNPs for chrom {chrom}...")
            gwas[COLS.drop] = gwas[[COLS.effect, COLS.other]].apply(lambda x: x[COLS.effect] == COMPLEMENT[x[COLS.other]], axis=1)

        else:
            print(f"\t\tSaving palindromic SNPs for chrom {chrom}...")
            save_palindromic_snps(gwas, vcf, strand, palindromic, pop_code)
            reason = CLEAN_STEPS.save_palindr


        gwas_bad_temp = gwas[gwas[COLS.drop]]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[~gwas[COLS.drop]]

        drop_map.palindromic_lost = (gwas_bad_temp.shape[0], reason)
        
    else:
        drop_map.palindromic_lost = (0, "Alleles Unknown.")

    del vcf
    
    return gwas, gwas_bad, drop_map
    
    
def map_parallel(shared_objects, chrom, out_dir, from_id, palindromic, pop_code, harmonize):
        (gwas_chr, gwas_bad, drop_map) = map_to_ensembl(chrom, shared_objects[chrom], out_dir, from_id, palindromic, pop_code, harmonize)
        
        result = dict()
        result["gwas"] = gwas_chr
        result["gwas_bad"] = gwas_bad
        result["drop_map"] = drop_map
        shared_objects[chrom] = result
    
    
def map_to_ensembl_procedure(gwas, drop_map, out_dir, from_id, palindromic, pop_code, harmonize):
    
    with multiprocessing.Manager() as manager:
        shared_objects = manager.dict()
        chrom_list = [f"rs{i}" for i in range(0,20)]
        by_chrom = COLS.chr in gwas.columns

        if by_chrom: 
            gwas[COLS.chr] = gwas[COLS.chr].astype(str)
            chrom_list = gwas[COLS.chr].unique()

            for chrom in chrom_list:
                shared_objects[chrom] = gwas[gwas[COLS.chr] == chrom].copy()

        else:
            interval = np.ceil(gwas.shape[0]/len(chrom_list)).astype(int)
            ranges = list(np.arange(0, gwas.shape[0], interval))
            if ranges[-1] < gwas.shape[0]:
                ranges.append(gwas.shape[0])

            for chrom in chrom_list:
                shared_objects[chrom] = gwas.iloc[ranges[int(chrom[2:])]:ranges[int(chrom[2:])+1]].copy()

        jobs = [multiprocessing.Process(target=map_parallel, args=(shared_objects, chrom, out_dir, from_id, palindromic, pop_code, harmonize)) for chrom in chrom_list]

        half = np.floor(len(jobs)/2).astype(int)
        
        for job in jobs[:half]:
            job.start()

        for job in jobs[:half]:
            job.join()

        gwas = pd.concat([shared_objects[chrom]["gwas"] for chrom in chrom_list[:half]], ignore_index = True)
        gwas_bad = pd.concat([shared_objects[chrom]["gwas_bad"] for chrom in chrom_list[:half]], ignore_index = True)
        
        for job in jobs[half:]:
            job.start()

        for job in jobs[half:]:
            job.join()

        gwas = pd.concat([gwas] + [shared_objects[chrom]["gwas"] for chrom in chrom_list[half:]], ignore_index = True)
        gwas_bad = pd.concat([gwas_bad] + [shared_objects[chrom]["gwas_bad"] for chrom in chrom_list[half:]], ignore_index = True)

        ensembl_drop_map = dict()

        palindromic_lost_reason = ""
        total_update = 0
        total_palindromic_lost = 0 
        total_palindromic = 0

        for chrom in chrom_list:
            chrom_map = shared_objects[chrom]["drop_map"]
            total_update += chrom_map.update
            total_palindromic += chrom_map.palindromic_total
            total_palindromic_lost += chrom_map.palindromic_lost[0]
            palindromic_lost_reason = chrom_map.palindromic_lost[1]
            if by_chrom:
                ensembl_drop_map[chrom] = chrom_map

        if not by_chrom:
            chrom_map = shared_objects[chrom_list[0]]["drop_map"]
            chrom_map.update = total_update
            chrom_map.palindromic_total = total_palindromic
            chrom_map.palindromic_lost = (total_palindromic_lost, chrom_map.palindromic_lost[1])
            ensembl_drop_map["all"] = chrom_map

        drop_map.map_ensembl = ensembl_drop_map

        drop_map.total_update = (total_update, CLEAN_STEPS.update)
        drop_map.total_palindromic = (total_palindromic, CLEAN_STEPS.total_palindromic)
        drop_map.total_palindromic_lost = (total_palindromic_lost, palindromic_lost_reason)

    return gwas, gwas_bad
    
    
def standardize_snps(gwas, drop_map, build_in, build_out, out_dir, from_id, palindromic, pop_code, harmonize):
    gwas_bad = pd.DataFrame(columns = gwas.columns)

    if build_in != G_BUILDS.hg38 and not from_id:
        print("\tDoing first liftover step...")
        gwas = liftover_procedure(build_in, G_BUILDS.hg38, gwas)
    
        gwas_bad = gwas[gwas[COLS.drop]]
        gwas = gwas[~gwas[COLS.drop]]
        
        drop_map.lift1 = (gwas_bad.shape[0], CLEAN_STEPS.lift1)
        
    else:
        drop_map.lift1 = (0, CLEAN_STEPS.lift1)
        
    print(f"\tMapping to Ensembl...")    
    
    (gwas, gwas_bad_temp) = map_to_ensembl_procedure(gwas, drop_map, out_dir, from_id, palindromic, pop_code, harmonize)
    
    gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
    
    if build_out != G_BUILDS.hg38:
        print("\tDoing second liftover step...")
        
        gwas[COLS.chr] = gwas[COLS.chr].astype(str)
        
        gwas = liftover_procedure(G_BUILDS.hg38, build_out, gwas, True)
        
        gwas_bad_temp = gwas[gwas[COLS.drop]]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[~gwas[COLS.drop]]

        drop_map.lift2 = (gwas_bad_temp.shape[0], CLEAN_STEPS.lift2)
        
        gwas[COLS.chr] = gwas[COLS.chr].apply(lambda x: CHROM_CODE[x])
        
    else:
        drop_map.lift2 = (0, CLEAN_STEPS.lift2)
        
    gwas.sort_values([COLS.chr, COLS.bp], inplace = True)

    return gwas, gwas_bad
        
        
def liftover(build_in, build_out, gwas, chrom, code):
    
    if code:
        chrom = f"chr{CHROM_VAL[int(chrom)]}"
    else:
        chrom = f"chr{CHROM_VAL[CHROM_CODE[chrom]]}"
        
    chain = os.path.join(PATHS.chain_files, f"{build_in}to{build_out}.chain")
    lo = LiftOver(chain)

    bp = gwas[COLS.bp].array
    drop = gwas[COLS.drop].array
    chroms = gwas[COLS.chr].array
    effect = gwas[COLS.effect].array if COLS.effect in gwas.columns else False
    other = gwas[COLS.other].array if COLS.other in gwas.columns else False
    alleles = effect and other
    
    for i in range(bp.size):
        pos_new = lo.convert_coordinate(chrom, bp[i])
        
        if len(pos_new) == 0:
            drop[i] = True
            continue
            
        elif len(pos_new) > 1:
            j = 0
            best = 0
            for index, tup in pos_new:
                if tup[3] > best:
                    best = tup[3]
                    j = index
            pos_new = pos_new[j:j+1]
            
        pos_new = pos_new[0]
        if pos_new[0][3:] in CHROM_CODE.keys():
            chroms[i] = str(pos_new[0][3:])
            bp[i] = int(pos_new[1])
        else:
            drop[i] = True
        
        if pos_new[2] == "-" and alleles:
            effect[i] = COMPLEMENT[effect[i]]
            other[i] = COMPLEMENT[other[i]]
        
    return gwas


def liftover_parallel(shared_objects, build_in, build_out, chrom, code):
    gwas = liftover(build_in, build_out, shared_objects[chrom], chrom, code)
    shared_objects[chrom] = gwas

    
def liftover_procedure(build_in, build_out, gwas, code=False):
    with multiprocessing.Manager() as manager:
        shared_objects = manager.dict()
        chrom_list = gwas[COLS.chr].unique()

        for chrom in chrom_list:
            shared_objects[chrom] = gwas[gwas[COLS.chr] == chrom].copy()

        jobs = [multiprocessing.Process(target=liftover_parallel, args=(shared_objects, build_in, build_out, chrom, code)) for chrom in chrom_list]

        for job in jobs:
            job.start()

        for job in jobs:
            job.join()

        gwas = pd.concat([shared_objects[chrom] for chrom in chrom_list], ignore_index = True)
    
    return gwas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="path to JSON file containing script arguments")

    json_input = parser.parse_args()
    
    with open(json_input.path) as f:
        args = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
    
    if not args.rsid and not args.parse_cols.rsid:
        REQUIRED.remove(COLS_OUT.rsid)
        REQUIRED += [COLS_OUT.chr, COLS_OUT.bp]
        
    if not args.harmonize:
        REQUIRED.remove(COLS_OUT.effect)
        REQUIRED.remove(COLS_OUT.other)

    check_required_cols(args)

    date = ''.join(str(datetime.date(datetime.now())).split('-'))
    
#    input_path = os.path.join(PATHS.data, args.input_path)
    input_path = args.input_path
    out_dir = os.path.join(args.output_path, args.trait)
    small_p_out_path = os.path.join(out_dir, "small_p_temp.txt")
    readme_file_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_{date}.README")
    description_file_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_qc_checks_{date}.txt")
    bad_gwas_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_dropped_snps_{date}.tsv")
    temp_plot_data_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_temp_plot_data_{date}.tsv")
    std_clean_data_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_clean_data_standardized_{date}.tsv")
    all_clean_data_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_clean_data_original_{date}.tsv")
    manhattan_output_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_manhattan_{date}")
    qq_output_path = os.path.join(out_dir, f"{args.trait}_{args.bld_out}_{args.pubid}_{args.initials}_qqplot_{date}.jpg")

    
    print("Importing dataset...")
        
    gwas = pd.read_csv(input_path, sep=args.sep, index_col=False)
    
    check_column_names(gwas.columns, vars(args))

    old_cols = gwas.columns
        
    create_dir(out_dir, args.overwrite)
    
    with open(readme_file_path, "w") as file:
        if not args.effect or not args.other:
            file.write(f"This GWAS output can't be harmonized to 1000 Genomes!\n")
        file.write(f"Array: {args.array}\n")
        file.write(f"Input genome build: {args.bld_in}\n")
        if args.effect and args.other:
            file.write(f"Effect/Noneffect column names(in): {args.effect}/{args.other}\n")
            file.write(f"Effect/Noneffect column names(out): {COLS.effect}/{COLS.other}\n")
        else:
            file.write(f"Effect/Noneffect column not provided on don't exist.\n")
        file.write(f"Analyst's name: {args.initials}\n")
        file.write(f"N max: {args.N}\n")
        file.write(f"dbSNP: {args.dbSNP}\n")
        
        
    print("Extracting small pvalues...")
    p_col_num = gwas.columns.get_loc(args.p) + 1
    
    cmd = f"Rscript '{PATHS.p_val_script_path}' '{args.input_path}' '{small_p_out_path}' '{args.sep}' True {p_col_num}"
    sys_call(os.system, cmd, "Small pvalue script failed to run")

    
    print("Importing new file with small pvalues...")
    gwas = pd.read_csv(small_p_out_path, sep=args.sep, index_col=False)
    gwas.columns = old_cols
    
    drop_map = SimpleNamespace(**{}) 
    drop_map.total = gwas.shape[0]
    
    if args.beta:
        drop_map.mean_effect = (gwas[gwas[args.p] < 0.001][args.beta].mean(skipna=True), "Beta")
    elif args.odds:
        drop_map.mean_effect = (gwas[gwas[args.p] < 0.001][args.odds].mean(skipna=True), "OR")
    else:
        drop_map.mean_effect = ("unknown")
    
    sys_call(os.remove, small_p_out_path, "Couldn't delete output of small P-value script")

    
    print("Parsing, cleaning columns and standardizing column names...")
    
    parse_args = args.parse_cols
    cols_saved = []
                
    if parse_args.required:
        
        if parse_args.full:
            
            if not args.chr:
                cols_saved.append(COLS.chr)   
                args.chr = COLS.chr
            if not args.bp:
                cols_saved.append(COLS.bp)   
                args.bp = COLS.bp
            if not args.effect:
                cols_saved.append(COLS.effect)   
                args.effect = COLS.effect
            if not args.other:
                cols_saved.append(COLS.other)   
                args.other = COLS.other
                
            df_one = gwas[parse_args.required].str.split("_", expand=True).rename(columns={0:parse_args.required, 1:COLS.effect, 2:COLS.other})
            df_two = df_one[parse_args.required].str.split(":", expand=True).rename(columns={0:COLS.chr, 1:COLS.bp})
            df_one.drop([parse_args.required], axis=1, inplace=True)
            df = df_one.join(df_two)
            
            gwas = gwas.join(df[cols_saved])
            gwas.dropna(subset=cols_saved, inplace=True)
                
        else:
            col_names = dict()   
                
            if parse_args.chr and not args.chr:
                col_names[parse_args.chr - 1] = COLS.chr
                cols_saved.append(COLS.chr)   
                args.chr = COLS.chr
            if parse_args.bp and not args.bp:
                col_names[parse_args.bp - 1] = COLS.bp
                cols_saved.append(COLS.bp)    
                args.bp = COLS.bp   
            if parse_args.rsid and not args.rsid:
                col_names[parse_args.rsid - 1] = COLS.rsid
                cols_saved.append(COLS.rsid)       
                args.rsid = COLS.rsid
            if parse_args.effect and not args.effect:
                col_names[parse_args.effect - 1] = COLS.effect
                cols_saved.append(COLS.effect)       
                args.effect = COLS.effect
            if parse_args.other and not args.other:
                col_names[parse_args.other - 1] = COLS.other
                cols_saved.append(COLS.other)   
                args.other = COLS.other
            
            df = gwas[parse_args.required].str.split(parse_args.sep, expand=True).rename(columns=col_names)
            
            gwas = gwas.join(df[cols_saved])
            gwas.dropna(subset=cols_saved, inplace=True)
        
    gwas[COLS.drop] = False

    print("\tCleaning bad SNPs based on:")

    print(f"\t\tP-values...")
            
    gwas.dropna(axis=0, subset=[args.p], inplace=True)
    gwas[COLS.drop] = gwas[args.p].apply(lambda x: not is_float(x) or (x <= 0 or x >= 1))
    gwas_bad = gwas[gwas[COLS.drop]]
    gwas = gwas[~gwas[COLS.drop]]
    
    drop_map.p = (gwas_bad.shape[0], CLEAN_STEPS.p)
    
    if args.rsid:
        print("\t\tSNP IDs...")  
        gwas.dropna(axis=0, subset=[args.rsid], inplace=True)
        gwas[COLS.drop] = gwas[args.rsid].apply(lambda x: not is_string(x) or not x.startswith("rs") or not x[2:].isnumeric())
        gwas_bad_temp = gwas[gwas[COLS.drop]]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[~gwas[COLS.drop]]
        
        drop_map.rsid = (gwas_bad_temp.shape[0], CLEAN_STEPS.rsid)
    
    else:
        drop_map.rsid = (0, CLEAN_STEPS.rsid)
        
    if args.chr:
        print("\t\tChromosomes...")
        gwas.dropna(axis=0, subset=[args.chr], inplace=True)
        gwas[args.chr] = gwas[args.chr].astype(str)
        gwas[args.chr] = gwas[args.chr].apply(lambda x: x.split(".")[0] if "." in x else x)
        gwas[args.chr] = gwas[args.chr].apply(lambda x: x[3:] if x.startswith("chr") else (x))
        gwas[COLS.drop] = gwas[args.chr].apply(lambda x: x not in CHROM_CODE.keys())
        gwas_bad_temp = gwas[gwas[COLS.drop]]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[~gwas[COLS.drop]]
        
        drop_map.chr = (gwas_bad_temp.shape[0], CLEAN_STEPS.chr)
        
    else:
        drop_map.chr = (0, CLEAN_STEPS.chr)
   
    if not args.rsid:
        print(f"\t\tBase Pair...")
        gwas.dropna(axis=0, subset=[args.bp], inplace=True)
        gwas[COLS.drop] = gwas[args.bp].apply(lambda x: not str(x).isnumeric())
        gwas_bad_temp = gwas[gwas[COLS.drop]]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[~gwas[COLS.drop]]
        
        drop_map.bp = (gwas_bad_temp.shape[0], CLEAN_STEPS.bp)
        gwas[args.bp] = gwas[args.bp].astype(int)
                                             
    else:
        drop_map.bp = (0, CLEAN_STEPS.bp)
        
    
    print(f"\t\tAlleles...")
    if args.effect and args.other:
        gwas[args.effect] = gwas[args.effect].apply(classify_snp)
        gwas[args.other] = gwas[args.other].apply(classify_snp)

        gwas_bad_temp = gwas[(gwas[args.effect] == "bad") | (gwas[args.other] == "bad")]
        gwas_bad_temp_2 = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        drop_map.bad_allele = (gwas_bad_temp.shape[0], CLEAN_STEPS.bad_allele)

        gwas_bad_temp = gwas[(gwas[args.effect] == "non_biallelic") | (gwas[args.other] == "non_biallelic")]
        gwas_bad_temp_2 = pd.concat([gwas_bad_temp_2, gwas_bad_temp], ignore_index=True)
        drop_map.non_biallelic = (gwas_bad_temp.shape[0], CLEAN_STEPS.non_biallelic)

        gwas_bad_temp = gwas[(gwas[args.effect] == "indel") | (gwas[args.other] == "indel")]
        gwas_bad_temp_2 = pd.concat([gwas_bad_temp_2, gwas_bad_temp], ignore_index=True)
        drop_map.indel = (gwas_bad_temp.shape[0], CLEAN_STEPS.indel)

        drop_map.allele_total = (drop_map.indel[0] + drop_map.non_biallelic[0] + drop_map.bad_allele[0], CLEAN_STEPS.allele_total)
        
        if gwas_bad_temp_2.shape[0] > gwas.shape[0] * 0.7:
            args.harmonize = False
            args.effect = False
            args.other = False
            drop_map.bad_allele = (0, "Alleles are bad")
            drop_map.non_biallelic = (0, "Alleles are bad")
            drop_map.indel = (0, "Alleles are bad")
        
        else:
            gwas_bad = gwas_bad_temp_2
            gwas = gwas[(gwas[args.effect].str.len() == 1) & (gwas[args.other].str.len() == 1)]
        
    else:
        drop_map.bad_allele = (0, "Alleles Unknown")
        drop_map.non_biallelic = (0, "Alleles Unknown")
        drop_map.indel = (0, "Alleles Unknown")
        drop_map.allele_total = (0, "Alleles Unknown")
        
    drop_map.pre_clean = (gwas_bad.shape[0], CLEAN_STEPS.pre_clean)
    
    old_cols = gwas_bad.columns
    
                                            
    print("\tStandardizing column names")
    
    gwas[COLS.p] = gwas[args.p]
                        
    if args.effect and args.effect not in cols_saved:
        gwas[COLS.effect] = gwas[args.effect]       
    if args.other and args.other not in cols_saved:
        gwas[COLS.other] = gwas[args.other]        
    if args.rsid and args.rsid not in cols_saved:
        gwas[COLS.rsid] = gwas[args.rsid]
    if args.chr and args.chr not in cols_saved:
        gwas[COLS.chr] = gwas[args.chr]
    if args.bp and args.bp not in cols_saved:
        gwas[COLS.bp] = gwas[args.bp]
    if args.beta:
        gwas[COLS.beta] = gwas[args.beta]
    if args.odds:
        gwas[COLS.odds] = gwas[args.odds]
    if args.uci:
        gwas[COLS.uci] = gwas[args.uci]
    if args.lci:
        gwas[COLS.lci] = gwas[args.lci]
    if args.eaf:
        gwas[COLS.eaf] = gwas[args.eaf]
    
    print("Standardizing SNPs...")
    
    (gwas, gwas_bad_temp) = standardize_snps(gwas, drop_map, args.bld_in, args.bld_out, out_dir, args.rsid, args.palindromic, args.pop_code, args.harmonize)
    
    gwas_bad = pd.concat([gwas_bad, gwas_bad_temp[old_cols]], ignore_index=True)
    
    drop_map.standardize = (gwas_bad_temp.shape[0], CLEAN_STEPS.standardize)
    
    
    print("Checking for non-autosomal chromosomes...")
    if not args.keep_X:
        gwas_bad_temp = gwas[gwas[COLS.chr] == CHROM_CODE["X"]][old_cols]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[gwas[COLS.chr] != CHROM_CODE["X"]]
        
        drop_map.chrX = (gwas_bad_temp.shape[0], CLEAN_STEPS.chrX)
        
    if not args.keep_Y:
        gwas_bad_temp = gwas[gwas[COLS.chr] == CHROM_CODE["Y"]][old_cols]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[gwas[COLS.chr] != CHROM_CODE["Y"]]
        
        drop_map.chrY = (gwas_bad_temp.shape[0], CLEAN_STEPS.chrY)
        
    if not args.keep_M:
        gwas_bad_temp = gwas[gwas[COLS.chr] == CHROM_CODE["M"]][old_cols]
        gwas_bad = pd.concat([gwas_bad, gwas_bad_temp], ignore_index=True)
        gwas = gwas[gwas[COLS.chr] != CHROM_CODE["M"]]
        
        drop_map.chrM = (gwas_bad_temp.shape[0], CLEAN_STEPS.chrM)

        
    print("Describing GWAS...")
    with open(description_file_path,'w') as file:
        
        file.write(f"Input GWAS had {drop_map.total} SNPs.") 
        file.write(f" The average SNP had {drop_map.mean_effect[0]} effect size {f'({drop_map.mean_effect[1]})' if args.beta or args.odds else ''}. \n")  
        
        if args.beta:
            drop_map.mean_effect = (gwas[gwas[COLS.p] < 0.001][COLS.beta].mean(skipna=True), "Beta")
        elif args.odds:
            drop_map.mean_effect = (gwas[gwas[COLS.p] < 0.001][COLS.odds].mean(skipna=True), "OR")
        else:
            drop_map.mean_effect = ("unknown")
            
        file.write(f"Output GWAS has {gwas.shape[0]} ({round_to_3(gwas.shape[0]/drop_map.total*100)}%) SNPs.")
        file.write(f" The average SNP has {drop_map.mean_effect[0]} effect size {f'({drop_map.mean_effect[1]})' if args.beta or args.odds else ''}. \n")
        
        file.write(f"\n\nBreakdown of dropped SNPs\n")
        
        ORDER = ["pre_clean", "p", "rsid", "chr", "bp", "allele_total", "bad_allele", "non_biallelic", "indel", "standardize", "lift1", "total_update"]
        
        for key in ORDER:
            write_drop_breakdown(file, vars(drop_map)[key], drop_map.total)
        
        for chrom in drop_map.map_ensembl.keys():
            snp_count = drop_map.map_ensembl[chrom].update
            file.write(f"\t\t\tchr{chrom}: {snp_count} {f'({round_to_3(snp_count/drop_map.total*100)}%)' if snp_count > 0 else ''}\n")
            
        write_drop_breakdown(file, drop_map.total_palindromic_lost, drop_map.total)
        for chrom in drop_map.map_ensembl.keys():
            snp_count = drop_map.map_ensembl[chrom].palindromic_lost[0]
            file.write(f"\t\t\tchr{chrom}: {snp_count} {f'({round_to_3(snp_count/drop_map.total*100)}%)' if snp_count > 0 else ''}\n")
            
        ORDER = ["lift2"]
        
        if not args.keep_X:
            ORDER.append("chrX")
        if not args.keep_Y:
            ORDER.append("chrY")
        if not args.keep_M:
            ORDER.append("chrM")
        
        for key in ORDER:
            write_drop_breakdown(file, vars(drop_map)[key], drop_map.total)
        
        file.write(f"\n\nChromosome distribution\n")
        
        chrom_dist = gwas[COLS.chr].value_counts()
        for chrom in np.sort(gwas[COLS.chr].unique()):
            file.write(f'\tchr{CHROM_VAL[chrom]} {chrom_dist.loc[chrom]} ({round_to_3(chrom_dist.loc[chrom]/gwas.shape[0]*100)}%)\n')
            
        write_drop_breakdown(file, drop_map.total_palindromic, drop_map.total)
        
        if args.effect:
            file.write(f"\n\nEffect allele\n")
            file.write(f"{gwas[COLS.effect].value_counts().to_string()}\n")
        
        if args.other:
            file.write(f"\nOther allele\n")
            file.write(f"{gwas[COLS.other].value_counts().to_string()}\n")
        
        if args.effect and args.other:
            file.write('\n\nAllele Combination Table\n')
            for allele in gwas[COLS.effect].unique():
                file.write(f"\n{allele}\n")
                file.write(f"{gwas[gwas[COLS.effect]==allele][COLS.other].value_counts().to_string()}\n")
    
    gwas[COLS.chr] = gwas[COLS.chr].apply(lambda x: CHROM_VAL[x])
        
        
    print("Saving GWAS files...")
    
    gwas_bad.drop(COLS.drop, axis = 1, inplace = True)
    gwas_bad.to_csv(bad_gwas_path, sep="\t", index=False)
    
    del gwas_bad, gwas_bad_temp
    
    gwas.drop(COLS.drop, axis = 1, inplace = True)
    gwas[COLS.p] = gwas[COLS.p].astype(float)
    gwas.to_csv(all_clean_data_path, sep="\t", index=False)
    
    orig_cols = [args.p]
    
    if args.effect and args.effect not in cols_saved:
        orig_cols.append(args.effect)
    if args.other and args.other not in cols_saved:
        orig_cols.append(args.other)
    if args.rsid and args.rsid not in cols_saved:
        orig_cols.append(args.rsid)
    if args.chr and args.chr not in cols_saved:
        orig_cols.append(args.chr)
    if args.bp and args.bp not in cols_saved:
        orig_cols.append(args.bp)
        
    if args.beta:
        orig_cols.append(args.beta)
    if args.odds:
        orig_cols.append(args.odds)
    if args.uci:
        orig_cols.append(args.uci)
    if args.lci:
        orig_cols.append(args.lci)
    if args.eaf:
        orig_cols.append(args.eaf)
        
    gwas.drop(orig_cols, axis = 1, inplace = True)
    final_cols = [col.split("_")[0] if col.endswith("_JD_STANDARD") else col for col in gwas.columns]
    gwas.columns = final_cols
    gwas[COLS_OUT.p] = gwas[COLS_OUT.p].astype(float)
    gwas.to_csv(std_clean_data_path, sep="\t", index=False)
    
    gwas = gwas[(gwas[COLS_OUT.p] > 0) & (gwas[COLS_OUT.p] < 1)]
    gwas[[COLS_OUT.p, COLS_OUT.rsid, COLS_OUT.chr, COLS_OUT.bp]].to_csv(temp_plot_data_path, sep="\t", index=False)
    
    bonferronni = -math.log10(0.05/gwas.shape[0])

    
    print("Creating manhattan plot...")
    cmd = f"python3 {PATHS.manhattan_script_path} \
    --twopoint {temp_plot_data_path} \
    --col-chr {COLS_OUT.chr} \
    --col-name {COLS_OUT.rsid} \
    --col-pos {COLS_OUT.bp} \
    --col-pvalue {COLS_OUT.p} \
    --bp \
    --use-pvalues \
    --abline {bonferronni} \
    --significant-threshold {bonferronni} \
    --no-annotation \
    --significant-point-size 2 \
    --point-size 1 \
    --graph-title '{args.trait} GWAS Manhattan Plot' \
    --chr-text-size 10 \
    --output {manhattan_output_path}"

    sys_call(os.system, cmd, "Couldn't generate manhattan plot.")
    
    print("Creating QQ plot...")

    mpl.rcParams['figure.dpi']=100
    mpl.rcParams['savefig.dpi']=100
    mpl.rcParams['figure.figsize']=5.375, 5.375

    qqplot.qqplot([gwas[COLS_OUT.p]],  
        [f'{args.trait} GWAS'], 
        color=['b'], 
        fill_dens=[0.2], 
        error_type='theoretical', 
        distribution='beta',
        title=f'{args.trait} GWAS QQ Plot')

    mpl.pyplot.legend(loc=0)
    plt.savefig(qq_output_path, dpi=300)
    
    sys_call(os.remove, temp_plot_data_path, "Couldn't remove temp data file for plotting.")

 
