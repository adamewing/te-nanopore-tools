#!/usr/bin/env python3

from collections import defaultdict as dd
from collections import Counter

import os
import pysam
import argparse

import multiprocessing as mp

import pandas as pd
import numpy as np

from uuid import uuid4
import gzip


import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class Read:
    def __init__(self, read_name, cpg_loc, llr):
        self.read_name = read_name
        self.llrs = {}
        self.meth_calls = {}

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


def fetch_reads(fn, chrom, start, end, min_mapq=10):
    reads = []

    bam = pysam.AlignmentFile(fn)

    for read in bam.fetch(chrom, start, end):
        if read.mapq >= min_mapq:
            reads.append(read.query_name)

    return reads


def get_meth_calls(bam_fn, meth_fn, chrom, seg_start, seg_end, seg_name):
    meth_tbx = pysam.Tabixfile(meth_fn)

    if chrom not in meth_tbx.contigs:
        return None, (chrom, seg_start, seg_end, seg_name)

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

    reads = fetch_reads(bam_fn, chrom, seg_start, seg_end)

    reads = list(set(reads).intersection(set(methdata.index)))

    methdata = methdata.loc[reads]

    reads = {}

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
                if index not in reads:
                    reads[index] = Read(index, cg_seg_start, llr)
                else:
                    reads[index].add_cpg(cg_seg_start, llr)

    meth_calls = dd(int)

    for name, read in reads.items():
        for loc, call in read.meth_calls.items():
            meth_calls[call] += 1

    logger.info('finished segment: %s %s:%d-%d' % (bam_fn, chrom, seg_start, seg_end))

    return meth_calls, (chrom, seg_start, seg_end, seg_name)


stats = [
'meth_calls',
'unmeth_calls',
'no_calls',
'methfrac'
]

def main(args):

    data = {}

    with open(args.data) as _:
        for line in _:
            bam, meth = line.strip().split()
            data[bam] = meth

    base_names = ['.'.join(bam.split('.')[:-1]) for bam in data]

    header = ['bin', 'seg_name']

    for bn in base_names:
        for s in stats:
            header.append('%s_%s' % (bn, s))

    print('\t'.join(header))

    pool = mp.Pool(processes=int(args.procs))
    bin_size = int(args.binsize)

    results = []

    for bam_fn, meth_fn in data.items():
        base_name = '.'.join(bam_fn.split('.')[:-1])

        with open(args.fai) as fai:
            for line in fai:
                chrom, chrlen = line.strip().split()[:2]
                chrlen = int(chrlen)

                for bin_start in range(0, chrlen, bin_size):
                    bin_end = bin_start + bin_size

                    if bin_end > chrlen:
                        bin_end = chrlen

                    bin_name = chrom

                    res = pool.apply_async(get_meth_calls, [bam_fn, meth_fn, chrom, bin_start, bin_end, bin_name])

                    results.append((res, base_name))


    meth_segs = dd(dict)

    for res, base_name in results:
        meth_data, seg = res.get()

        if meth_data is None:
            continue

        seg_id = '%s:%d-%d' % seg[:3]

        seg_name = seg[3]

        no_calls = 0
        meth_calls = 0
        unmeth_calls = 0

        if -1 in meth_data:
            unmeth_calls = meth_data[-1]

        if 0 in meth_data:
            no_calls = meth_data[0]

        if 1 in meth_data:
            meth_calls = meth_data[1]

        meth_segs[seg_id]['seg_name'] = seg_name
        meth_segs[seg_id][base_name + '_meth_calls'] = meth_calls
        meth_segs[seg_id][base_name + '_unmeth_calls'] = unmeth_calls
        meth_segs[seg_id][base_name + '_no_calls'] = no_calls
        meth_segs[seg_id][base_name + '_methfrac'] = 0.0
        if meth_calls + unmeth_calls > 0:
            meth_segs[seg_id][base_name + '_methfrac'] = float(meth_calls)/float(meth_calls + unmeth_calls)

    for seg in meth_segs:
        output = [seg]
        for h in header[1:]:
            if h not in meth_segs[seg]:
                logger.warning('Warning: no output for %s in segment %s, skipped segment' % (h, seg))
                break

            output.append(str(meth_segs[seg][h]))
        print('\t'.join(output))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-d', '--data', required=True, help='text file with .bam filename and corresponding nanopolish call-methylation tabix (whitespace-delimited)')
    parser.add_argument('-p', '--procs', default=1, help='multiprocessing')
    parser.add_argument('-b', '--binsize', required=True, help='bin size')
    parser.add_argument('-f', '--fai', required=True, help='fasta index (.fai)')

    args = parser.parse_args()
    main(args)
