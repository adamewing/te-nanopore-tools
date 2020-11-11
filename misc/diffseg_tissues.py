#!/usr/bin/env python3

from collections import defaultdict as dd
from itertools import product

import os
import pysam
import argparse
import statsmodels.stats.multitest as ssm

import pandas as pd
import numpy as np
import scipy.stats as ss

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):
    data = pd.read_csv(args.segmeth, sep='\t', header=0, index_col=0)

    #samples = {}

    #for c in data.columns:
    #    if c.endswith('_meth_calls'):
    #        samples[c.replace('_meth_calls', '')] = True

    #assert len(samples) == 3

    samples = ['hc5413_merged_tagreads','li5413_merged_tagreads','he5413_merged_tagreads']

    for s in samples:
        data[s+'_methfrac'] = (data[s+'_meth_calls']/(data[s+'_meth_calls']+data[s+'_unmeth_calls']))


    for col in data.columns:
        if col.startswith('hcc'):
            data.drop(col, axis=1)

    useable = []

    for seg in data.index:
        use_seg = True

        for s in samples:
            if data[s+'_meth_calls'].loc[seg] < 5 or data[s+'_unmeth_calls'].loc[seg] < 5:
                use_seg = False
                continue

        if use_seg:
            useable.append(seg)

    data = data.loc[useable]

    data['var'] = np.var(data[[s+'_methfrac' for s in samples]], axis=1)


    # need e.g. k-means for > 3 samples

    dist = {}

    for i, j in product(samples, repeat=2):
        if i != j:
            cmp_name = ','.join((i,j))
            cmp_name2 = ','.join((j,i))

            if cmp_name not in dist and cmp_name2 not in dist:
                dist[cmp_name] = abs(data[i+'_methfrac'] - data[j+'_methfrac'])

    #print(dist.keys())

    diff_sample = []
    fisher_p = []
    genes = []
    gdists = []
    gtypes = []

    hc_exprs = []
    he_exprs = []
    li_exprs = []

    highest_exprs = []

    gtex = pysam.Tabixfile('gtexGene.hg38.txt.gz')

    for seg in data.index:
        min_cmp = None
        min_cmp_pair = None

        for mf_cmp in dist:
            if min_cmp is None or min_cmp > dist[mf_cmp].loc[seg]:
                min_cmp = dist[mf_cmp].loc[seg]
                min_cmp_pair = mf_cmp.split(',')

        max_cmp_name = [s for s in samples if s not in min_cmp_pair][0]
        diff_sample.append(max_cmp_name)

        s0_meth   = data[max_cmp_name+'_meth_calls'].loc[seg]
        s0_unmeth = data[max_cmp_name+'_unmeth_calls'].loc[seg]

        s1_meth   = sum([data[m+'_meth_calls'].loc[seg] for m in min_cmp_pair])
        s1_unmeth = sum([data[m+'_unmeth_calls'].loc[seg] for m in min_cmp_pair])

        c_table = [[s0_meth,s0_unmeth],[s1_meth,s1_unmeth]]

        o, p = ss.fisher_exact(c_table)

        fisher_p.append(p)

        # gene info

        chrom = data.loc[seg]['seg_chrom'] 
        start = data.loc[seg]['seg_start']
        end = data.loc[seg]['seg_end']

        gene = []
        gdist = []
        gtype = []

        hc_exp = []
        he_exp = []
        li_exp = []

        high_exp = []

        if chrom in gtex.contigs:
            start_window = start-100000
            if start_window < 1:
                start_window = 1

            for rec in gtex.fetch(chrom, start_window, end+100000):
                c = rec.split()

                intronic = False


                if start > int(c[1]) and end < int(c[2]):
                    intronic = True

                d = 0
                if not intronic:
                    d = min(abs(int(c[1]) - start), abs(int(c[1])-end), abs(int(c[2])-start), abs(int(c[2])-end))

                # pick the closest protein-coding gene (non-coding if no protein-coding gene within 50kbp)
                if len(gdist) > 0:
                    if d < int(gdist[-1]) or (gtype[-1] != 'protein_coding' and c[7] == 'protein_coding'):
                        gene = []
                        gdist = []
                        gtype = []
                        hc_exp = []
                        he_exp = []
                        li_exp = []
                        high_exp = []

                    # if multiple genes at gdist == 0 (intronic), pick the one with the highest max expression
                    elif d == int(gdist[-1]) == 0:
                        hc = float(exprs[14])
                        he = float(exprs[32])
                        li = float(exprs[34])

                        max_expr = max(hc, he, li)
                        if max_expr > max(float(hc_exp[-1]), float(he_exp[-1]), float(li_exp[-1])):
                            gene = []
                            gdist = []
                            gtype = []
                            hc_exp = []
                            he_exp = []
                            li_exp = []
                            high_exp = []

                        else:
                            continue

                    else:
                        continue

                gdist.append(str(d))

                if c[3] not in gene:
                    gene.append(c[3])
                    gtype.append(c[7])

                    exprs = c[-1].split(',')

                    # GTEx tissues are ordered alphabetically in UCSC table
                    hc_exp.append(exprs[14]) # hippocampus
                    he_exp.append(exprs[32]) # left ventricle
                    li_exp.append(exprs[34]) # liver

                    hc = float(exprs[14])
                    he = float(exprs[32])
                    li = float(exprs[34])

                    if max(hc, he, li) < 1.0:
                        high_exp.append('NA')

                    elif max(hc, he, li) == hc:
                        high_exp.append('hc')

                    elif max(hc, he, li) == he:
                        high_exp.append('he')

                    elif max(hc, he, li) == li:
                        high_exp.append('li')

        if len(gene) == 0:
            gene = ['NA']
            gdist = ['NA']
            gtype = ['NA']
            hc_exp = ['0']
            he_exp = ['0']
            li_exp = ['0']
            high_exp = ['NA']

        genes.append(','.join(gene))
        gdists.append(','.join(gdist))
        gtypes.append(','.join(gtype))

        hc_exprs.append(','.join(hc_exp))
        he_exprs.append(','.join(he_exp))
        li_exprs.append(','.join(li_exp))

        highest_exprs.append(','.join(high_exp))


    data['diff_sample']  = diff_sample
    data['fisher_p']     = fisher_p
    data['Bonferroni_p'] = ssm.multipletests(data['fisher_p'], alpha=0.05, method='bonferroni')[1]
    data['nearest_gene'] = genes
    data['gene_dist']    = gdists
    data['gene_type']    = gtypes
    data['hc_exp']       = hc_exprs
    data['he_exp']       = he_exprs
    data['li_exp']       = li_exprs
    data['highest_exp']  = highest_exprs

    data.to_csv(args.segmeth + '.diffseg_tissues.tsv', sep='\t', na_rep='NA')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-s', '--segmeth', required=True, help='output from segmeth.py')


    args = parser.parse_args()
    main(args)
