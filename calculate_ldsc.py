#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu


Usage: ./calculate_ldsc.py \
--vcf /storage/s1saini/hipstr_allfilters/phased_feb18/hipstr.chr22.phased.vcf.gz \
--bim /storage/s1saini/pgc_analysis/ldsc/calculate_ldsc/str_bim/str.chr22.bial.bim \
--out . \
--ld-wind-cm 1 \
--region 22:16788134-17788134


To do:
1. Add support for STRs
    - Done
2. Output M and M_5_50 file
    - Done
3. Add support for annotations needed for partioned heritability
    - Done

"""
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import allel

import argparse
import sys

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def normalize_gt(X):
    ii = ~np.isnan(X)
    avg = np.mean(X[ii])
    X[np.logical_not(ii)] = avg
    std = np.std(X)
    if std==0:
        std=1
    normalized_gt = ((X-avg)/std)
    return normalized_gt

def calcMAF(counts):
    return np.sum(np.sort(counts)[:-1])


def __l2_unbiased__(x, n):
    denom = n-2 if n > 2 else n  # allow n<2 for testing purposes
    sq = np.square(x)
    return sq - (1-sq) / denom


def snp_getter(b, current_snp, bim_data, VCF_FILE, samples=None):
    CHROM = bim_data['CHR'][0]
    POS = list(bim_data['BP'].values)[current_snp:current_snp+b]
    start_bp = POS[0]
    end_bp = POS[-1]
    vcf_window = '%d:%d-%d'%(CHROM,start_bp,end_bp)
    genotypes = []

    vcf = VCF(VCF_FILE, samples=samples)
    for v in vcf(vcf_window):
        alleles, counts = np.unique(np.array(v.genotypes)[:,0:2].flatten(), return_counts=True)
        if 1 in counts: ### exclude singletons
            continue

        if len(vcf.samples) in counts: ### Exclude variants with MAF=0
            continue

        genotypes.append(normalize_gt(list(np.sum(np.array(v.genotypes)[:,0:2], axis=1))))


    return np.array(genotypes).T, current_snp+b


def snp_getter_from_array(b, current_snp, geno_array):
    return geno_array[:, current_snp:current_snp+b], current_snp+b


def getBlockLefts(coords, max_dist):
    M = len(coords)
    j = 0
    block_left = np.zeros(M)
    for i in range(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1
        block_left[i] = j
    return (block_left)


def isSTR(record_alleles_to_bases):
    if record_alleles_to_bases[0] == 1 and record_alleles_to_bases[1]==1 and len(record_alleles_to_bases.keys())==2:
        return False
    else:
        return True


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="Input VCF File", required=True, type=str)
    parser.add_argument("--bim", help="Input BIM file with CM information, if using with --ld-wind-cm", required=False, type=str)
    parser.add_argument("--out", help="Output directory name", required=True, type=str)
    parser.add_argument("--ld-wind-cm", help="Specify the window size to be used for estimating LD Scores in units of centiMorgans (cM)", required=False, type=float)
    parser.add_argument("--ld-wind-kb", help="Specify the window size to be used for estimating LD Scores in units of kilobase-pairs (kb)", required=False, type=int)
    parser.add_argument("--ld-wind-snps", help="Specify the window size to be used for estimating LD Scores in units of # of SNPs.", required=False, type=int)
    parser.add_argument("--keep", help="File with individuals to include in LD Score estimation", required=False, type=str)
    parser.add_argument("--chunk-size", help="Chunk size for LD Score calculation. Use the default.", required=False, type=int, default=50)
    parser.add_argument("--annot", help="Filename prefix for annotation file for partitioned LD Score estimation.", required=False, type=str)
    parser.add_argument("--region", help="Fetch only chr:start-end from the VCF", required=False, type=str)
    parser.add_argument("--min-maf", help="Minimum MAF to use for filtering alleles. Default: 0.05", required=False, type=float, default=0.05)
    args = parser.parse_args()



    VCF_FILE = args.vcf

    if args.keep:
        SAMPLES_FILE = args.keep
        samples = []
        with open(SAMPLES_FILE) as f:
            for line in f:
                samples.append(line.strip().split()[0])
        print("Restricting to %d samples"%len(samples))
    else:
        samples=None

    chunk_size = args.chunk_size
    directory = args.out

    if args.region:
        region = args.region
    else:
        region=None

    x = np.array((args.ld_wind_kb, args.ld_wind_cm, args.ld_wind_snps), dtype=bool)
    if np.sum(x) != 1:
        sys.exit('Must specify exactly one --ld-wind option')


    if args.ld_wind_cm:
        if not args.bim:
            print("--ld-wind-cm requires a BIM file with CM information")
            sys.exit("Terminated")
        else: ### Load the BIM file ###
            BIM_FILE = args.bim
            bim_data = pd.read_csv(BIM_FILE, names=["CHR","SNP","cM","BP","REF","ALT"], delim_whitespace=True)
            bim_snps = list(bim_data['SNP'].values)


### Load meta data Information per record
### Used for finding STR allele lengths
### Used to create bim_data DF when BIM file not used
    print("Loading meta data information")
    meta_data = []
    vcf = VCF(VCF_FILE)
    if args.region:
        vcf_iter = vcf(args.region)
    else:
        vcf_iter = vcf()
    alleles_to_bases = {}
    for v in vcf_iter:
        alleles = [len(v.REF)]+[len(i) for i in v.ALT]
        alleles_to_bases[str(v.ID)] = dict(zip(range(len(alleles)), alleles))
        meta_data.append([v.CHROM, v.ID, v.POS])

    if not args.bim:
        bim_data = pd.DataFrame(meta_data, columns=["CHR","SNP","BP"])
        bim_snps = list(bim_data['SNP'].values)


### Load Genotypes from the VCF ###
    print("Loading Genotypes")
    keep_snps = []
    M_out = 0
    M_5_50_out = 0
    callset = allel.read_vcf(VCF_FILE, region=region, samples=samples, fields=['samples', 'variants/ID', 'calldata/GT'])
    genotype_array = callset['calldata/GT']
    genotypes = np.zeros((genotype_array.shape[0],genotype_array.shape[1]))
    count = 0

    for i in range(genotype_array.shape[0]):
        if i%1000 == 0:
            print("Processed %d variants"%(i))
        ID = callset['variants/ID'][i]
        if ID not in bim_snps:
            continue
        if ID in keep_snps:
            continue

        record_alleles_to_bases = alleles_to_bases[ID]
        is_str = isSTR(record_alleles_to_bases)

        if not is_str:
            numSamples = callset['samples'].shape[0]
            alleles, counts = np.unique(genotype_array[i,:,:], return_counts=True)
            if 1 in counts: ### exclude singletons
                continue
            if callset['samples'].shape[0] in counts: ### Exclude variants with MAF=0
                continue
            maf = calcMAF(counts/numSamples)
            if maf >= 0.05:
                M_out += 1
                M_5_50_out += 1
            else:
                M_out += 1
            genotypes[count] = (normalize_gt(np.sum(genotype_array[i,:,:], axis=1)))
            keep_snps.append(ID)
            count += 1

        else:
            numSamples = callset['samples'].shape[0]
            alleles, counts = np.unique(genotype_array[i,:,:], return_counts=True)
            counts = counts/numSamples
            maf = calcMAF(counts)
            if maf >= 0.05:
                M_out += 1
                M_5_50_out += 1
            else:
                M_out += 1
            rm_alleles = alleles[np.argwhere(counts<args.min_maf)]
            rm_alleles = np.insert(rm_alleles, 0, -1) ## remove missing alleles
            filtered_gt = np.where(np.isin(genotype_array[i,:,:], rm_alleles), np.nan, genotype_array[i,:,:]) # remove alleles with MAF < args.min-maf
            u,inv = np.unique(filtered_gt,return_inverse = True) # map allele index to allele length
            filtered_gt = np.array([record_alleles_to_bases.get(x,x) for x in u])[inv].reshape(filtered_gt.shape) # map allele index to allele length
            genotypes[count] = (normalize_gt(np.sum(filtered_gt, axis=1)))
            keep_snps.append(ID)
            count += 1

    keep_snps = list(set(keep_snps))
    genotypes = genotypes[0:count]
    genotypes = np.array(genotypes).T
    bim_data = bim_data[bim_data.SNP.isin(keep_snps)]
    bim_data.drop_duplicates(subset='SNP', inplace=True)

    if bim_data.shape[0] != genotypes.shape[1]:
        sys.exit("The meta data object must contain the same SNPs")


### Find the left most snp to be included in LD score calculation of snp i ###
    if args.ld_wind_cm:
        coords = list(bim_data['cM'].values)
        block_left = getBlockLefts(coords, args.ld_wind_cm)
    elif args.ld_wind_kb:
        coords = list(bim_data['BP'].values)
        block_left = getBlockLefts(coords, args.ld_wind_kb*1000)
    elif args.ld_wind_snps:
        coords = np.array(range(bim_data.shape[0]))
        block_left = getBlockLefts(coords, args.ld_wind_snps)




### Loading Annotation File ###
    if not args.annot:
        annot=None
    else:
        print("Loading Annotation File")
        annot_df = pd.read_csv(args.annot, compression="gzip", delim_whitespace=True)
        annot_df = annot_df[annot_df.SNP.isin(keep_snps)]
        annot_df.drop_duplicates(subset='SNP', inplace=True)

        annot_df.SNP = annot_df.SNP.astype("category")
        annot_df.SNP.cat.set_categories(keep_snps, inplace=True)
        annot_df.sort_values(["SNP"], inplace=True)

        n_annot, ma = len(annot_df.columns) - 4, len(annot_df)
        annot = np.array(annot_df.iloc[:,4:])
        annot_colnames = annot_df.columns[4:]

        print("Loaded %d annotations for %d variants"%(n_annot, ma))

        if annot.shape[0] != genotypes.shape[1]:
            sys.exit("The .annot file must contain the same SNPs")



##### Following code is borrowed from the original LDSC codebase
##### The code finds LD in windows
#####

    print("Calculating LD Scores")

    m, n = len(coords), callset['samples'].shape[0]
    current_snp = 0
    block_sizes = np.array(np.arange(m) - block_left)
    block_sizes = np.ceil(block_sizes / chunk_size)*chunk_size
    if annot is None:
        annot = np.ones((m, 1))


    n_a = annot.shape[1]
    cor_sum = np.zeros((m, n_a))
    b = np.nonzero(block_left > 0)
    if np.any(b):
        b = b[0][0]
    else:
        b = m
    b = int(np.ceil(b/chunk_size)*chunk_size)
    if b > m:
        c = 1
        b = m
    l_A = 0
    A, current_snp = snp_getter_from_array(b, current_snp, genotypes)


    rfuncAB = np.zeros((b, chunk_size))
    rfuncBB = np.zeros((chunk_size, chunk_size))

    for l_B in range(0, b, chunk_size):
        B = A[:, l_B:l_B+chunk_size]
        np.dot(A.T, B / n, out=rfuncAB)
        rfuncAB = __l2_unbiased__(rfuncAB, n)
        cor_sum[l_A:l_A+b, :] += np.dot(rfuncAB, annot[l_B:l_B+chunk_size, :])

    b0 = b
    md = int(chunk_size*np.floor(m/chunk_size))
    end = md + 1 if md != m else md
    for l_B in range(b0, end, chunk_size):
        old_b = b
        b = int(block_sizes[l_B])
        if l_B > b0 and b > 0:
            A = np.hstack((A[:, old_b-b+chunk_size:old_b], B))
            l_A += old_b-b+chunk_size
        elif l_B == b0 and b > 0:
            A = A[:, b0-b:b0]
            l_A = b0-b
        elif b == 0:  # no SNPs to left in window, e.g., after a sequence gap
            A = np.array(()).reshape((n, 0))
            l_A = l_B
        if l_B == md:
            chunk_size = m - md
            rfuncAB = np.zeros((b, chunk_size))
            rfuncBB = np.zeros((chunk_size, chunk_size))
        if b != old_b:
            rfuncAB = np.zeros((b, chunk_size))

        B, current_snp = snp_getter_from_array(chunk_size, current_snp, genotypes)
        p1 = np.all(annot[l_A:l_A+b, :] == 0)
        p2 = np.all(annot[l_B:l_B+chunk_size, :] == 0)
        if p1 and p2:
            continue

        np.dot(A.T, B / n, out=rfuncAB)
        rfuncAB = __l2_unbiased__(rfuncAB, n)
        cor_sum[l_A:l_A+b, :] += np.dot(rfuncAB, annot[l_B:l_B+chunk_size, :])
        cor_sum[l_B:l_B+chunk_size, :] += np.dot(annot[l_A:l_A+b, :].T, rfuncAB).T
        np.dot(B.T, B / n, out=rfuncBB)
        rfuncBB = __l2_unbiased__(rfuncBB, n)
        cor_sum[l_B:l_B+chunk_size, :] += np.dot(rfuncBB, annot[l_B:l_B+chunk_size, :])
######
######
######

    print("Writing output file")
    output_cols = []
    if args.annot:
        for i in range(len(annot_colnames)):
            colname=annot_colnames[i]+"L2"
            output_cols.append(colname)
            bim_data[colname] = cor_sum[:,i]
    else:
        bim_data['L2'] = cor_sum.flatten()
        output_cols.append('L2')

    bim_data_output = bim_data[['CHR','SNP','BP']+output_cols]


    filename = "%s/%d.l2.ldscore.gz"%(directory,int(np.unique(bim_data_output['CHR'])))
    bim_data_output.to_csv(filename, sep="\t", index=False, compression='gzip')


    filename = "%s/%d.l2.M"%(directory,int(np.unique(bim_data_output['CHR'])))
    f = open(filename, "w")
    f.write("%d"%M_out)
    f.close()

    filename = "%s/%d.l2.M_5_50"%(directory,int(np.unique(bim_data_output['CHR'])))
    f = open(filename, "w")
    f.write("%d"%M_5_50_out)
    f.close()

if __name__ == "__main__":
    main()
