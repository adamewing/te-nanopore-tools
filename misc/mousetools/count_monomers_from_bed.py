#!/usr/bin/env python

import os
import subprocess
import argparse
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
        self.seqid     = c[1]
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


def rc(dna):
    ''' reverse complement '''
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def align(te_seq, monomer_fa):
    rnd = str(uuid4())

    te_fa = 'tmp.' + rnd + '.te.fa'
    with open(te_fa, 'w') as te:
        te.write('>te\n%s\n' % te_seq)

    cmd = ['exonerate', '--bestn', '1', '--model', 'affine:local', '--showalignment','0', '--ryo', 'TE\t%qi\t%s\t%qab\t%qae\t%tab\t%tae\t%pi\t%qS\t%tS\n', monomer_fa, te_fa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = None

    for pline in p.stdout.readlines():
        pline = pline.decode()
        if pline.startswith('TE'):
            a = ExonerateAlignment(pline.strip())

            if best is None:
                best = a

            else:
                if a.score > best.score:
                    best = a
    
    os.remove(te_fa)

    return best


def load_falib(infa):
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip().split()[0]
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def count_monomers(args, te_seq):
    monomer_intervals = []

    while True:
        m = align(te_seq, args.monomers)

        if m is None:
            break

        # mask monomer
        te_seq = list(te_seq)
        te_seq[m.tgt_start:m.tgt_end] = 'N'*(m.tgt_end-m.tgt_start)
        te_seq = ''.join(te_seq)

        monomer_intervals.append((m.seqid, m.tgt_start, m.tgt_end))

    monomer_intervals.sort(key=lambda x: x[1])

    return monomer_intervals


def main(args):
    ref = pysam.Fastafile(args.ref)
    results = []

    pool = mp.Pool(processes=int(args.procs))

    with open(args.bed) as bed:
        for line in bed:
            c = line.strip().split()
            chrom, start, end = c[:3]
            start = int(start)
            end = int(end)

            te_seq = ref.fetch(chrom, start, end)

            if c[-1] in ('+', '-'):
                if c[-1] == '-':
                    te_seq = rc(te_seq)

            res = pool.apply_async(count_monomers, [args, te_seq])
            results.append((res, line))

    for res, line in results:
        monomer_intervals = res.get()
        logger.info('processed %s' % line.strip())
        m_out = ['%s:%d-%d' % m for m in monomer_intervals]
        print('%s\t%d\t%s' % (line.strip(), len(monomer_intervals), ','.join(m_out)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-m', '--monomers', required=True, help='monomer fasta')
    parser.add_argument('-r', '--ref', required=True, help='samtools faidx-indexed reference genome fasta')
    parser.add_argument('-b', '--bed', required=True)
    parser.add_argument('-p', '--procs', default=1)
    args = parser.parse_args()
    main(args)
