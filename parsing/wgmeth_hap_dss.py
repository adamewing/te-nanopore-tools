#!/usr/bin/env python3

from __future__ import print_function

from collections import defaultdict as dd
from collections import Counter

import os
import pysam
import argparse

import multiprocessing as mp

import pandas as pd
import numpy as np
import scipy.stats as ss

from uuid import uuid4
import gzip

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Read:
    def __init__(self, read_name, cpg_loc, llr, phase=None):
        self.read_name = read_name
        self.llrs = {}
        self.meth_calls = {}
        self.phase = phase

        self.add_cpg(cpg_loc, llr)

    def add_cpg(self, cpg_loc, llr, cutoff = 2.5):
        #assert abs(llr) > cutoff

        self.llrs[cpg_loc] = llr

        if llr > cutoff:
            self.meth_calls[cpg_loc] = 1
        elif llr < -1*cutoff:
            self.meth_calls[cpg_loc] = -1
        else:
            self.meth_calls[cpg_loc] = 0


def exclude_ambiguous_reads(fn, chrom, start, end, min_mapq=10, tag_untagged=False):
    reads = {}

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch(chrom, start, end):
        p = read.get_reference_positions()
        if p[0] < start or p[-1] > end:
            if read.mapq >= min_mapq:
                phase = None

                if tag_untagged:
                    phase = 'unphased'

                HP = None
                PS = None
                for tag in read.get_tags():
                    if tag[0] == 'HP':
                        HP = tag[1]
                    if tag[0] == 'PS':
                        PS = tag[1]

                if None not in (HP, PS):
                    phase = str(PS) + ':' + str(HP)

                reads[read.query_name] = phase

    return reads


def get_reads(fn, chrom, start, end, min_mapq=10, tag_untagged=False):
    reads = {}

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch(chrom, start, end):
        if read.mapq >= min_mapq:    
            phase = None

            if tag_untagged:
                phase = 'unphased'

            HP = None
            PS = None
            for tag in read.get_tags():
                if tag[0] == 'HP':
                    HP = tag[1]
                if tag[0] == 'PS':
                    PS = tag[1]

            if None not in (HP, PS):
                phase = str(PS) + ':' + str(HP)

            reads[read.query_name] = phase

    return reads


def get_meth_calls(bam_fn, meth_fn, chrom, seg_start, seg_end):
    meth_tbx = pysam.Tabixfile(meth_fn)

    if chrom not in meth_tbx.contigs:
        return None

    tmp_methdata = str(uuid4()) +'.tmp.methdata.tsv'

    with open(tmp_methdata, 'w') as meth_out:
        # header
        with gzip.open(meth_fn, 'rt') as _:
            for line in _:
                assert line.startswith('chromosome')
                meth_out.write(line)
                break

        for rec in meth_tbx.fetch(chrom, seg_start, seg_end):
            meth_out.write(str(rec)+'\n')

    # index by read_name
    methdata = pd.read_csv(tmp_methdata, sep='\t', header=0, index_col=4)

    os.remove(tmp_methdata)

    reads = {}
    if args.excl_ambig:
        reads = exclude_ambiguous_reads(bam_fn, chrom, seg_start, seg_end)
    else:
        reads = get_reads(bam_fn, chrom, seg_start, seg_end)

    readnames = []
    for r in reads.keys():
        if r in methdata.index:
            readnames.append(r)

    methdata = methdata.loc[readnames]


    seg_methreads = {}

    for index, row in methdata.iterrows():
        r_start = row['start']
        r_end   = row['end']
        llr     = row['log_lik_ratio']
        seq     = row['sequence']

        # get per-CG position (nanopolish/calculate_methylation_frequency.py)
        cg_pos = seq.find("CG")
        first_cg_pos = cg_pos

        while cg_pos != -1:
            cg_start = r_start + cg_pos - first_cg_pos
            cg_pos = seq.find("CG", cg_pos + 1)

            cg_seg_start = cg_start - seg_start

            if cg_start >= seg_start and cg_start <= seg_end:
                if index not in seg_methreads:
                    seg_methreads[index] = Read(index, cg_seg_start, llr, phase=reads[index])
                else:
                    seg_methreads[index].add_cpg(cg_seg_start, llr)


    meth_data = {}
    meth_data[1] = dd(list)
    meth_data[2] = dd(list)

    for name, read in seg_methreads.items():
        if read.phase is None:
            continue

        if read.phase.split(':')[1] not in ['1','2']:
            continue

        for loc in read.llrs.keys():
            phase = int(read.phase.split(':')[1])
            assert phase in (1,2)

            meth_data[phase][loc].append(read.meth_calls[loc])

    meth_table = {}
    meth_table[1] = dd(dict)
    meth_table[2] = dd(dict)

    for phase in (1,2):
        for loc in meth_data[phase]:
            pos = loc+seg_start
            meth_table[phase][pos]['chr'] = chrom
            meth_table[phase][pos]['N'] = len([call for call in meth_data[phase][loc] if call != 0])
            meth_table[phase][pos]['X'] = len([call for call in meth_data[phase][loc] if call == 1])


    return [meth_table[1], meth_table[2]]



def main(args):

    meth_table = [None, None]
    seg_meth_table_store = dd(list)

    pool = mp.Pool(processes=int(args.procs))
    bin_size = int(args.binsize)

    results = []

    with open(args.fai) as fai:
        for line in fai:
            chrom, chrlen = line.strip().split()[:2]
            chrlen = int(chrlen)

            for seg_start in range(0, chrlen, bin_size):
                seg_end = seg_start + bin_size

                if seg_end > chrlen:
                    seg_end = chrlen

                seg_start = int(seg_start)
                seg_end = int(seg_end)

                res = pool.apply_async(get_meth_calls, [args.bam, args.methdata, chrom, seg_start, seg_end])

                results.append(res)


    for res in results:
        seg_meth_table = res.get()

        if seg_meth_table is None:
            continue

        for phase in (0,1):

            seg_meth_table_store[phase].append(pd.DataFrame.from_dict(seg_meth_table[phase]).T)


    for phase in (0,1):
        meth_table[phase] = pd.concat(seg_meth_table_store[phase])


    for phase in (0,1):
        meth_table[phase]['pos'] = meth_table[phase].index
        meth_table[phase] = meth_table[phase].sort_values(['chr', 'pos'])

        outfn = '.'.join(os.path.basename(args.bam).split('.')[:-1]) + '.phase_%d.DSS.txt' % phase

        logger.info('writing %s' % outfn)
        meth_table[phase].to_csv(outfn, columns=['chr','pos','N','X'], index=False, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-b', '--bam', required=True, help='bam used for nanopolish methylation calling')
    parser.add_argument('-m', '--methdata', required=True, help='whole genome nanopolish methylation tabix')
    parser.add_argument('-s', '--binsize', required=True, help='bin size')
    parser.add_argument('-f', '--fai', required=True, help='fasta index (.fai)')
    parser.add_argument('-p', '--procs', default=1, help='multiprocessing')
    parser.add_argument('--excl_ambig', action='store_true', default=False)

    args = parser.parse_args()
    main(args)
