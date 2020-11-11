#!/usr/bin/env python

import sys
import csv
import argparse
import pandas as pd
import numpy as np
import scipy.stats as ss

from collections import defaultdict as dd
from itertools import product

def main(args):
    data = {}

    with open(args.data) as d:
        csv_reader = csv.DictReader(d)

        for rec in csv_reader:
            if rec['group'] not in data:
                data[rec['group']] = dd(list)

            data[rec['group']][rec['sample']].append(float(rec['mCpG']))

    stats = dd(dict)

    for group in list(data.keys()):
        print("Group %s has size %d" % (group, len(data[group])))

        if len(data[group]) > 2:
            s_array = [data[group][sample] for sample in data[group]]
            H, p = ss.kruskal(*s_array)
            print("Kruskal-Wallis H: %.3f p: %.3E" % (H, p))

        else:
            print("Skipping Kruskal-Wallis test on group size 2")

        print()
        table = dd(dict)
        for sample in list(data[group].keys()):
            table[sample]['group'] = group
            table[sample]['N'] = len(data[group][sample])
            table[sample]['median'] = np.median(data[group][sample])

            for other in list(data[group]):
                p = 1.0
                if other != sample:
                    u, p = ss.mannwhitneyu(data[group][sample], data[group][other])

                table[sample]['vs_'+other+'_pval'] = p

        table = pd.DataFrame.from_dict(table).T

        print('uncorrected pairwise Mann-Whitney U tests:')
        table.to_csv(sys.stdout, index_label='sample')
        print('\n')





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-d', '--data', required=True, help='csv file from segplot.py')
    args = parser.parse_args()
    main(args)