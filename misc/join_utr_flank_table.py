#!/usr/bin/env python

import argparse
import csv

import numpy as np

def main(args):

    utr_meth = {}

    with open(args.utrs) as utrs:
        csv_reader = csv.DictReader(utrs, delimiter='\t')

        for rec in csv_reader:
            if int(rec['hcc33nt_merged_tagreads_meth_calls']) + int(rec['hcc33nt_merged_tagreads_unmeth_calls']) < int(args.mincalls):
                continue

            if int(rec['hcc33t_merged_tagreads_meth_calls']) + int(rec['hcc33t_merged_tagreads_unmeth_calls']) < int(args.mincalls):
                continue

            if 0 in (int(rec['hcc33nt_merged_tagreads_meth_calls']), int(rec['hcc33nt_merged_tagreads_unmeth_calls']), int(rec['hcc33t_merged_tagreads_meth_calls']), int(rec['hcc33t_merged_tagreads_unmeth_calls'])):
                continue


            key = ''

            if rec['seg_strand'] == '+':
                key = rec['seg_chrom'] + ':' + rec['seg_start']

            if rec['seg_strand'] == '-':
                key = rec['seg_chrom'] + ':' + rec['seg_end']

            hcc33nt_frac = float(rec['hcc33nt_merged_tagreads_meth_calls'])/(float(rec['hcc33nt_merged_tagreads_meth_calls']) + float(rec['hcc33nt_merged_tagreads_unmeth_calls']))
            hcc33t_frac  = float(rec['hcc33t_merged_tagreads_meth_calls'])/(float(rec['hcc33t_merged_tagreads_meth_calls']) + float(rec['hcc33t_merged_tagreads_unmeth_calls']))

            utr_meth[key] = (rec['seg_id'],
                                int(rec['hcc33nt_merged_tagreads_meth_calls']),
                                int(rec['hcc33nt_merged_tagreads_unmeth_calls']), 
                                int(rec['hcc33t_merged_tagreads_meth_calls']), 
                                int(rec['hcc33t_merged_tagreads_unmeth_calls']), 
                                hcc33nt_frac, hcc33t_frac)


    columns = ('UTR', 'Flank', 'Subfamily', 'L1_Strand', 'hcc33nt_UTR_mCpG', 'hcc33nt_Flank_mCpG', 'hcc33t_UTR_mCpG', 'hcc33t_Flank_mCpG', 'UTR_T_NT_Ratio', 'Flank_T_NT_Ratio', 'D_Identity', 'z_Score')
    print('\t'.join(columns))
    output = []

    d_idents = []

    with open(args.flanks) as flanks:
        csv_reader = csv.DictReader(flanks, delimiter='\t')

        for rec in csv_reader:

            if int(rec['hcc33nt_merged_tagreads_meth_calls']) + int(rec['hcc33nt_merged_tagreads_unmeth_calls']) < int(args.mincalls):
                continue

            if int(rec['hcc33t_merged_tagreads_meth_calls']) + int(rec['hcc33t_merged_tagreads_unmeth_calls']) < int(args.mincalls):
                continue

            if 0 in (int(rec['hcc33nt_merged_tagreads_meth_calls']), int(rec['hcc33nt_merged_tagreads_unmeth_calls']), int(rec['hcc33t_merged_tagreads_meth_calls']), int(rec['hcc33t_merged_tagreads_unmeth_calls'])):
                continue

            key = ''

            if rec['seg_strand'] == '+':
                key = rec['seg_chrom'] + ':' + rec['seg_end']

            if rec['seg_strand'] == '-':
                key = rec['seg_chrom'] + ':' + rec['seg_start']

            flank_hcc33nt_frac = float(rec['hcc33nt_merged_tagreads_meth_calls'])/(float(rec['hcc33nt_merged_tagreads_meth_calls']) + float(rec['hcc33nt_merged_tagreads_unmeth_calls']))
            flank_hcc33t_frac  = float(rec['hcc33t_merged_tagreads_meth_calls'])/(float(rec['hcc33t_merged_tagreads_meth_calls']) + float(rec['hcc33t_merged_tagreads_unmeth_calls']))

            flank_hcc33nt_meth_count   = int(rec['hcc33nt_merged_tagreads_meth_calls'])
            flank_hcc33nt_unmeth_count = int(rec['hcc33nt_merged_tagreads_unmeth_calls'])
            flank_hcc33t_meth_count    = int(rec['hcc33t_merged_tagreads_meth_calls'])
            flank_hcc33t_unmeth_count  = int(rec['hcc33t_merged_tagreads_unmeth_calls'])


            if key in utr_meth:
                utr_id, utr_hcc33nt_meth_count, utr_hcc33nt_unmeth_count, utr_hcc33t_meth_count, utr_hcc33t_unmeth_count, utr_hcc33nt_frac, utr_hcc33t_frac = utr_meth[key]

                utr_tnt_ratio   = utr_hcc33t_frac/utr_hcc33nt_frac
                flank_tnt_ratio = flank_hcc33t_frac/flank_hcc33nt_frac

                d_ident = abs(utr_tnt_ratio-flank_tnt_ratio) / np.sqrt(2)

                if args.elts is not None:
                    if rec['seg_name'].split('_')[0] not in args.elts.split(','):
                        continue

                out_rec = [utr_id,
                    rec['seg_id'],
                    rec['seg_name'].split('_')[0],
                    rec['seg_strand'],
                    utr_hcc33nt_frac,
                    flank_hcc33nt_frac,
                    utr_hcc33t_frac,
                    flank_hcc33t_frac,
                    utr_tnt_ratio,
                    flank_tnt_ratio,
                    d_ident]

                d_idents.append(d_ident)
                output.append(out_rec)

    for out_rec in output:
        d = out_rec[-1]
        z = (d-np.mean(d_idents))/np.std(d_idents)
        out_rec.append(z)

        out_rec = map(str, out_rec)
        print('\t'.join(out_rec))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-u', '--utrs', required=True)
    parser.add_argument('-f', '--flanks', required=True)
    parser.add_argument('--mincalls', default=20)
    parser.add_argument('--elts', default=None, help='comma-delimited')
    args = parser.parse_args()
    main(args)