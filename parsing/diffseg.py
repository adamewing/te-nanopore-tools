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

    samples = args.samples.split(',')

    for s in samples:
        data[s+'_methfrac'] = (data[s+'_meth_calls']/(data[s+'_meth_calls']+data[s+'_unmeth_calls']))

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

    data['var'] = np.var(data[[s+'_methfrac' for s in samples]], axis=1)

    data.to_csv(args.segmeth + '.diffseg.tsv', sep='\t', na_rep='NA')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-s', '--segmeth', required=True, help='output from segmeth.py')
    parser.add_argument('-m', '--samples', required=True, help='samples, comma-delimited')
    parser.add_argument('-n', '--mincalls', default=10, help='minimum number of calls to include site (methylated + unmethylated) (default=10)')
    args = parser.parse_args()
    main(args)
