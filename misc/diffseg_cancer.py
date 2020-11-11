#!/usr/bin/env python3

from collections import defaultdict as dd

import os
import pysam
import argparse
import statsmodels.stats.multitest as ssm

import pandas as pd
import numpy as np
import scipy.stats as ss


def main(args):
    data = pd.read_csv(args.segmeth, sep='\t', header=0, index_col=0)

    #samples = {}

    #for c in data.columns:
    #    if c.endswith('_meth_calls'):
    #        samples[c.replace('_meth_calls', '')] = True

    #assert len(samples) == 2

    samples = ['hcc33t_merged_tagreads', 'hcc33nt_merged_tagreads']

    for s in samples:
        data[s+'_methfrac'] = data[s+'_meth_calls']/(data[s+'_meth_calls']+data[s+'_unmeth_calls'])


    for col in data.columns:
        if '5413' in col:
            data.drop(col, axis=1)

    useable = []

    for seg in data.index:
        use_seg = True

        for s in samples:
            if data[s+'_meth_calls'].loc[seg] < 10 or data[s+'_unmeth_calls'].loc[seg] < 10:
                use_seg = False
                continue

        if use_seg:
            useable.append(seg)

    data = data.loc[useable]


    
    f_res = []
    highmeth = []
    genes = []
    gdists = []
    gtypes = []

    gtex = pysam.Tabixfile('gtexGene.hg38.txt.gz')

    s = list(samples)[0]

    for seg in data.index:
        not_s = [n for n in samples if n != s][0]

        s_meth   = data[s+'_meth_calls'].loc[seg]
        s_unmeth = data[s+'_unmeth_calls'].loc[seg]

        not_s_meth   = data[not_s+'_meth_calls'].loc[seg]
        not_s_unmeth = data[not_s+'_unmeth_calls'].loc[seg]

        c_table = [[s_meth,s_unmeth],[not_s_meth,not_s_unmeth]]

        o, p = ss.fisher_exact(c_table)

        f_res.append(p)

        highest = s

        if data[s+'_methfrac'].loc[seg] < data[not_s+'_methfrac'].loc[seg]:
            highest = not_s

        highmeth.append(highest)

        # gene info
        #chrom, interval = seg.split(':')
        #start, end = map(int, interval.split('-'))

        chrom = data.loc[seg]['seg_chrom']
        start = data.loc[seg]['seg_start']
        end = data.loc[seg]['seg_end']

        gene = []
        gdist = []
        gtype = []

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

                # pick the closest protein-coding gene

                if len(gdist) > 0:
                    if d < int(gdist[-1]) or (gtype[-1] != 'protein_coding' and c[7] == 'protein_coding'):
                        gene = []
                        gdist = []
                        gtype = []

                    else:
                        continue

                gdist.append(str(d))

                if c[3] not in gene:
                    gene.append(c[3])
                    gtype.append(c[7])


        genes.append(','.join(gene))
        gdists.append(','.join(gdist))
        gtypes.append(','.join(gtype))

    data['highest_meth'] = highmeth
    data['fisher_p'] = f_res
    data['Bonferroni_p'] = ssm.multipletests(data['fisher_p'], alpha=0.05, method='bonferroni')[1]
    data['nearest_gene'] = genes
    data['gene_dist']    = gdists
    data['gene_type']    = gtypes
    data['methfrac_diff'] = data['hcc33t_merged_tagreads_methfrac'] - data['hcc33nt_merged_tagreads_methfrac']


    data.to_csv(args.segmeth + '.diffseg_cancer.tsv', sep='\t')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-s', '--segmeth', required=True, help='output from segmeth.py')


    args = parser.parse_args()
    main(args)
