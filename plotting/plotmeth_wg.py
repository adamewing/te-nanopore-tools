#!/usr/bin/env python3

import argparse

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

from collections import defaultdict as dd



def main(args):

    data = pd.read_csv(args.meth, sep='\t', header=0, index_col=0)

    sample_basenames = []

    if args.samplebasenames is None:

        for c in data.columns:
            if c.endswith('_methfrac'):
                sample_basenames.append(c.replace('_methfrac', ''))

    else:
        sample_basenames = args.samplebasenames.split(',')

    print(sample_basenames)
    
    chr_order = ['chr'+chr for chr in map(str, list(range(1,23)) + ['X'])]

    binstarts = dd(list)

    for chrom in chr_order:
        for seg in data[data['seg_name'] == chrom].index:
            start = int(seg.split(':')[1].split('-')[0])
            end   = int(seg.split(':')[1].split('-')[1])

            binsize = end-start

            if len(binstarts[chrom]) == 0:
                binstarts[chrom].append(start)

            else:
                binstarts[chrom].append(binstarts[chrom][-1]+binsize)


    color = {}

    cmap = plt.get_cmap('bwr', len(sample_basenames))

    for i, s in enumerate(sample_basenames):
        color[s] = cmap(i)

    fig = plt.figure()
    fig, ax = plt.subplots(len(chr_order), figsize=(12,3+0.65*len(chr_order)), sharex=True)

    ymin = 0.0
    ymax = 1.0

    if len(chr_order) == 1:
        ax = [ax]

    last_bs = 0
    longest_chrom = chrom

    for chrom in chr_order:
        for bs in binstarts[chrom]:
            if bs > last_bs:
                last_bs = bs
                longest_chrom = chrom

    print('longest chrom: %s' % longest_chrom)

    binsize = binstarts[longest_chrom][1] - binstarts[longest_chrom][0]

    print('inferred binsize: %d' % binsize)

    xmax = last_bs + binsize

    print('xmax: %d' % xmax)

    for i, chrom in enumerate(chr_order):

        ax[i].set_ylim(ymin, ymax)
        ax[i].set_xlim(0,xmax)

        ax[i].set_yticks([])
        ax[i].set_xticks([])
        ax[i].set_ylabel(chrom)

        for s in sample_basenames:
            ax[i].scatter(binstarts[chrom], data[data['seg_name'] == chrom][s+'_methfrac'], color=color[s], alpha=0.05, marker='.')

    ticklocs = list(range(0,xmax,int(20e6)))

    plt.xticks(ticklocs, map(lambda x: '%.1E' % x, ticklocs), rotation=45)

    out_fn = '.'.join(args.meth.split('.')[:-1]) + '.wgplot.png'

    if args.svg:
        out_fn = '.'.join(args.meth.split('.')[:-1]) + '.wgplot.svg'

    plt.savefig(out_fn, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='plotter')
    parser.add_argument('-m', '--meth', required=True, help='methylation segments')
    parser.add_argument('-s', '--samplebasenames', default=None, help='sample base names')
    parser.add_argument('--svg', action='store_true', default=False)

    args = parser.parse_args()
    main(args)

