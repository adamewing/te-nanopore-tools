#!/usr/bin/env python3

from collections import defaultdict as dd
from itertools import product

import os
import pysam
import argparse

import pandas as pd
import numpy as np
import scipy.stats as ss

import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

# Illustrator compatibility
new_rc_params = {'text.usetex': False, "svg.fonttype": 'none'}
matplotlib.rcParams.update(new_rc_params)

import matplotlib.pyplot as plt
import seaborn as sns

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):
    data = pd.read_csv(args.segmeth, sep='\t', header=0, index_col=0)

    samples = args.samples.split(',')

    for s in samples:
        data[s + '_methfrac'] = data[s+'_meth_calls']/(data[s+'_meth_calls']+data[s+'_unmeth_calls'])

    useable = []

    for seg in data.index:
        use_seg = True

        for s in samples:
            if (data[s+'_meth_calls'].loc[seg] + data[s+'_unmeth_calls'].loc[seg]) < int(args.mincalls):
                use_seg = False
                continue

        if use_seg:
            useable.append(seg)

    data = data.loc[useable]

    logger.info('useable sites: %d' % len(useable))

    #kw = ss.kruskal(*[data[s+'_methfrac'] for s in samples])
    #print(kw)

    plot_data = dd(dict)

    order = []

    for seg in data.index:
        for s in samples:
            uid = seg + ':' + s
            plot_data[uid]['sample'] = s
            plot_data[uid]['mCpG']   = data[s+'_methfrac'].loc[seg]
            plot_data[uid]['group']  = data['seg_name'].loc[seg]

            if plot_data[uid]['group'] not in order:
                order.append(plot_data[uid]['group'])

    plot_data = pd.DataFrame.from_dict(plot_data).T
    plot_data = pd.DataFrame(plot_data.to_dict())

    basename = '.'.join(args.segmeth.split('.')[:-1])

    plot_data.to_csv(basename+'.segplot_data.csv')
    logger.info('plot data written to %s.segplot_data.csv' % basename)

    #order = ['WG', 'L1HS', 'AluYa5,b8', 'SVA_E,F', 'LTR5_Hs']

    pt_sz = int(args.pointsize)

    if args.categories is not None:
        order = args.categories.split(',')

    if args.violin:
        sns_plot = sns.violinplot(x='group', y='mCpG', data=plot_data, hue='sample', dodge=True, jitter=True, order=order, hue_order=samples)

    else:
        sns_plot = sns.stripplot(x='group', y='mCpG', data=plot_data, hue='sample', dodge=True, jitter=True, size=pt_sz, order=order, hue_order=samples)

    if args.tiltlabel:
        sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=45)

    sns_plot.set_ylim(float(args.ymin),float(args.ymax))

    fig = sns_plot.figure
    fig.set_size_inches(int(args.width), int(args.height)) # TE

    if args.svg:
        fig.savefig(basename+'.segplot.svg', bbox_inches='tight')
    else:
        fig.savefig(basename+'.segplot.png', bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-s', '--segmeth', required=True, help='output from segmeth.py')
    parser.add_argument('-m', '--samples', required=True, help='samples, comma delimited')
    parser.add_argument('-c', '--categories', default=None, help='categories, comma delimited, need to match seg_name column from input')
    parser.add_argument('-v', '--violin', default=False, action='store_true')
    parser.add_argument('-n', '--mincalls', default=10, help='minimum number of calls to include site (methylated + unmethylated) (default=10)')
    parser.add_argument('--width', default=12, help='figure width (default = 12)')
    parser.add_argument('--height', default=6, help='figure height (default = 6)')
    parser.add_argument('--pointsize', default=1, help='point size for scatterplot (default = 1)')
    parser.add_argument('--ymin', default=-0.05, help='ymin (default = -0.05)')
    parser.add_argument('--ymax', default=1.05, help='ymax (default = 1.05)')
    parser.add_argument('--tiltlabel', default=False, action='store_true')
    parser.add_argument('--svg', default=False, action='store_true')

    args = parser.parse_args()
    main(args)
