#!/usr/bin/env python

import os
import argparse
import subprocess
import pysam

import multiprocessing as mp
from uuid import uuid4

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class ExonerateAlignment:
    def __init__(self, alignment):
        c = alignment.strip().split('\t')
        self.original  = c
        self.subfamily = c[1]
        self.score     = int(c[2])
        self.qry_start = int(c[3])
        self.qry_end   = int(c[4])
        self.tgt_start = int(c[5])
        self.tgt_end   = int(c[6])
        self.match     = float(c[7])
        self.orient    = c[9]

        if self.orient == '-':
            self.tgt_start, self.tgt_end = self.tgt_end, self.tgt_start

    def __lt__(self, other):
        return self.qry_start < other.qry_start


class AlignmentGroup:
    def __init__(self, alignment):
        self.alignments = []
        self.add_alignment(alignment)

    def add_alignment(self, alignment):
        if len(self.alignments) == 0:
            self.alignments.append(alignment)

        else:
            non_redundant = True

            for a in self.alignments:
                if alignment.qry_start > a.qry_start and alignment.qry_end < a.qry_end:
                    non_redundant = False

            if non_redundant:    
                self.alignments.append(alignment)

    def subfamily(self):
        return self.alignments[0].subfamily

    def best_score(self):
        return max([a.score for a in self.alignments])

    def best_match(self):
        return max([a.match for a in self.alignments])

    def min_qry_coord(self):
        return min([a.qry_start for a in self.alignments])

    def max_qry_coord(self):
        return max([a.qry_end for a in self.alignments])

    def min_tgt_coord(self):
        return min([a.tgt_start for a in self.alignments])

    def max_tgt_coord(self):
        return max([a.tgt_end for a in self.alignments])

    def min_qry_orient(self):
        ''' return orientation of leftmost alignment '''
        return sorted(self.alignments)[0].orient

    def max_qry_orient(self):
        ''' return orientation of rightmost alignment '''
        return sorted(self.alignments)[-1].orient


def te_align(te_fa, refseq, minmatch=80.0, te_fa_is_seq=False, max_segs=3):
    rnd = str(uuid4())
    tgt_fa = 'tmp.' + rnd + '.tgt.fa'

    with open(tgt_fa, 'w') as tgt:
        tgt.write('>ref\n%s\n' % refseq)

    if te_fa_is_seq:
        te_seq = te_fa
        te_fa = 'tmp.' + rnd + '.te.fa'
        with open(te_fa, 'w') as te:
            te.write('>te\n%s\n' % te_seq)

    cmd = ['exonerate', '--bestn', str(max_segs), '--model', 'affine:local', '--showalignment','0', '--ryo', 'TE\t%qi\t%s\t%qab\t%qae\t%tab\t%tae\t%pi\t%qS\t%tS\n', te_fa, tgt_fa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    align_groups = {}

    for pline in p.stdout.readlines():
        pline = pline.decode()
        if pline.startswith('TE'):
            a = ExonerateAlignment(pline.strip())
            if a.match >= minmatch:
                if a.subfamily in align_groups:
                    align_groups[a.subfamily].add_alignment(a)

                else:
                    align_groups[a.subfamily] = AlignmentGroup(a)

    os.remove(tgt_fa)

    if te_fa_is_seq:
        os.remove(te_fa)

    best_group = None

    for subfam in align_groups:
        if best_group is None:
            best_group = best_group = align_groups[subfam]

        elif align_groups[subfam].best_score() > best_group.best_score():
            best_group = align_groups[subfam]

    return best_group


def main(args):
    ref = pysam.Fastafile(args.ref)

    pool = mp.Pool(processes=int(args.procs))
    results = []

    with open(args.bed) as bed:
        for line in bed:
            c = line.strip().split()
            assert len(c) > 2

            chrom = c[0]
            start = int(c[1])
            end   = int(c[2])

            refseq = ref.fetch(chrom, start, end)

            res = pool.apply_async(te_align, [args.elts, refseq])
            results.append((res, c))


        for res, c in results:
            best_align = res.get()

            if best_align is None:
                logger.info('no good alignment for record: %s' % line.strip())
                continue

            logger.info('record: %s, best align: %s' % (str(c), best_align.subfamily()))
            c[3] = best_align.subfamily()
            c.append(str(best_align.best_match()))

            print('\t'.join(c))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='yee')
    parser.add_argument('-b', '--bed', required=True, help='bed file')
    parser.add_argument('-e', '--elts', help='reference elements .fa', required=True)
    parser.add_argument('-r', '--ref', help='samtools faidx-indexed ref', required=True)
    parser.add_argument('-p', '--procs', default=1, help='multiprocessing')
    args = parser.parse_args()
    main(args)
