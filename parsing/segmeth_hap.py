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


def get_meth_calls(bam_fn, meth_fn, chrom, seg_start, seg_end, seg_name, seg_strand):
    meth_tbx = pysam.Tabixfile(meth_fn)

    tmp_methdata = str(uuid4()) +'.tmp.methdata.tsv'

    with open(tmp_methdata, 'w') as meth_out:
        # header
        with gzip.open(meth_fn, 'rt') as _:
            for line in _:
                assert line.startswith('chromosome')
                meth_out.write(line)
                break

        assert chrom in meth_tbx.contigs

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

    meth_calls = dd(list)

    for name, read in seg_methreads.items():
        for loc, call in read.meth_calls.items():
            meth_calls[call].append(read.phase)

    return meth_calls, (chrom, seg_start, seg_end, seg_name, seg_strand)


stats = [
'meth_calls',
'meth_calls_phase1',
'meth_calls_phase2',
'unmeth_calls',
'unmeth_calls_phase1',
'unmeth_calls_phase2',
'phase_ratio',
'phase_p',
'no_calls'
]


def main(args):

    data = {}

    with open(args.data) as _:
        for line in _:
            bam, meth = line.strip().split()
            data[bam] = meth

    base_names = ['.'.join(bam.split('.')[:-1]) for bam in data]

    header = ['seg_id', 'seg_chrom', 'seg_start', 'seg_end', 'seg_name', 'seg_strand']

    for bn in base_names:
        for s in stats:
            header.append('%s_%s' % (bn, s))

    print('\t'.join(header))

    pool = mp.Pool(processes=int(args.procs))

    results = []

    for bam_fn, meth_fn in data.items():
        base_name = '.'.join(bam_fn.split('.')[:-1])

        with open(args.intervals) as _:
            for line in _:
                c = line.strip().split()
                chrom, seg_start, seg_end = c[:3]

                seg_name = 'NA'
                seg_strand = 'NA'

                if len(c) > 3:
                    seg_name = c[3]

                if len(c) > 4:
                    seg_strand = c[4]

                seg_start = int(seg_start)
                seg_end = int(seg_end)

                res = pool.apply_async(get_meth_calls, [bam_fn, meth_fn, chrom, seg_start, seg_end, seg_name, seg_strand])

                results.append((res, base_name))


    meth_segs = dd(dict)
    bad = {}

    for res, base_name in results:
        meth_data, seg = res.get()
        #print(meth_data)

        seg_id = '%s:%d-%d' % seg[:3]

        if seg_id in bad:
            logger.warning('segment %s marked as bad in another sample, skipping.' % seg_id)
            continue

        seg_chrom, seg_start, seg_end, seg_name, seg_strand = map(str, seg)

        no_calls = 0
        meth_calls = 0
        unmeth_calls = 0

        meth_calls_phase1 = 0
        meth_calls_phase2 = 0

        unmeth_calls_phase1 = 0
        unmeth_calls_phase2 = 0

        phases = []

        if -1 in meth_data:
            unmeth_calls = len(meth_data[-1])
            for p in meth_data[-1]:
                if p is not None:
                    if p not in phases:
                        phases.append(p)

                    if p.split(':')[-1] == '1':
                        unmeth_calls_phase1 += 1

                    if p.split(':')[-1] == '2':
                        unmeth_calls_phase2 += 1

        if 0 in meth_data:
            no_calls = len(meth_data[0])

        if 1 in meth_data:
            meth_calls = len(meth_data[1])
            for p in meth_data[1]:
                if p is not None:
                    if p not in phases:
                        phases.append(p)

                    if p.split(':')[-1] == '1':
                        meth_calls_phase1 += 1

                    if p.split(':')[-1] == '2':
                        meth_calls_phase2 += 1



        if len(phases) > 2:
            logger.warning('segment %s has > 2 phase calls, skipping.' % seg_id)
            bad[seg_id] = True
            continue
        meth_segs[seg_id]['seg_id'] = seg_id
        meth_segs[seg_id]['seg_chrom'] = seg_chrom
        meth_segs[seg_id]['seg_start'] = seg_start
        meth_segs[seg_id]['seg_end'] = seg_end
        meth_segs[seg_id]['seg_name'] = seg_name
        meth_segs[seg_id]['seg_strand'] = seg_strand

        meth_segs[seg_id][base_name + '_meth_calls'] = meth_calls
        meth_segs[seg_id][base_name + '_unmeth_calls'] = unmeth_calls

        meth_segs[seg_id][base_name + '_meth_calls_phase1'] = meth_calls_phase1
        meth_segs[seg_id][base_name + '_meth_calls_phase2'] = meth_calls_phase2

        meth_segs[seg_id][base_name + '_unmeth_calls_phase1'] = unmeth_calls_phase1
        meth_segs[seg_id][base_name + '_unmeth_calls_phase2'] = unmeth_calls_phase2

        fisher_ratio = 'NA'
        fisher_p = 'NA'


        if meth_calls_phase1 >= 10 and meth_calls_phase2 >=10 and unmeth_calls_phase1 >= 10 and unmeth_calls_phase2 >= 10:
            f_table = [[meth_calls_phase1, unmeth_calls_phase1], [meth_calls_phase2, unmeth_calls_phase2]]
            fisher_ratio, fisher_p = ss.fisher_exact(f_table)


        meth_segs[seg_id][base_name + '_phase_ratio'] = fisher_ratio
        meth_segs[seg_id][base_name + '_phase_p'] = fisher_p


        meth_segs[seg_id][base_name + '_no_calls'] = no_calls


    out_meth_segs = dd(dict)

    for seg_id in meth_segs:
        if seg_id not in bad:
            out_meth_segs[seg_id] = meth_segs[seg_id]

    for seg in out_meth_segs:
        output = []
        for h in header:
            output.append(str(out_meth_segs[seg][h]))
        print('\t'.join(output))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-d', '--data', required=True, help='text file with .bam filename and corresponding nanoport call-methylation tabix (whitespace-delimited)')
    parser.add_argument('-i', '--intervals', required=True, help='.bed file')
    parser.add_argument('-p', '--procs', default=1, help='multiprocessing')
    parser.add_argument('--excl_ambig', action='store_true', default=False)

    args = parser.parse_args()
    main(args)
