#!/usr/bin/env python3

from collections import defaultdict as dd
from itertools import product
from uuid import uuid4

import os
import gzip
import pysam
import argparse
import random

import multiprocessing as mp

import skbio.alignment as skalign
import skbio.sequence as skseq

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

from matplotlib import gridspec

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


def rc(dna):
    ''' reverse complement '''
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def exclude_ambiguous_reads(fn, chrom, start, end, min_mapq=10):
    reads = []

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch(chrom, start, end):
        p = read.get_reference_positions()
        if p[0] < start or p[-1] > end:
            if read.mapq >= min_mapq:
                reads.append(read.query_name)

    return reads


def single_seq_fa(fn):
    with open(fn, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                assert seq == '', 'input fa must have only one entry'
            else:
                seq = seq + line.strip()

    return seq


def get_reads(fn, chrom, start, end, min_mapq=10):
    reads = []

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch(chrom, start, end):
        if read.mapq >= min_mapq:
            reads.append(read.query_name)

    return reads


def slide_window(meth_table, sample, width=20, slide=2):
    midpt_min = min(meth_table['loc'])
    midpt_max = max(meth_table['loc'])

    sample_table = meth_table.loc[meth_table['sample'] == sample]

    win_start = int(midpt_min - width/2)
    win_end = win_start + width

    meth_frac = {}
    meth_n = {}

    while int((win_start+win_end)/2) < midpt_max:
        win_start += slide
        win_end += slide

        meth_count = len(meth_table.loc[(meth_table['sample'] == sample) & (meth_table['loc'] > win_start) & (meth_table['loc'] < win_end) & (meth_table['call'] == 1)])
        unmeth_count = len(meth_table.loc[(meth_table['sample'] == sample) & (meth_table['loc'] > win_start) & (meth_table['loc'] < win_end) & (meth_table['call'] == -1)])

        midpt = int((win_start+win_end)/2)

        if meth_count + unmeth_count > 0:
            meth_frac[midpt] = meth_count/(meth_count+unmeth_count)
            meth_n[midpt] = meth_count+unmeth_count

    return meth_frac, meth_n


def smooth(x, window_len=8, window='hanning'):
    ''' modified from scipy cookbook: https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html '''

    assert window_len % 2 == 0, '--smoothwindowsize must be an even number'
    assert x.ndim == 1
    assert x.size > window_len

    if window_len<3:
        return x

    assert window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    return y[(int(window_len/2)-1):-(int(window_len/2))]


def get_meth_profile(args, seg_chrom, seg_start, seg_end, seg_name, seg_strand):
    logger.info('profiling %s %s:%d-%d:%s' % (seg_name, seg_chrom, seg_start, seg_end, seg_strand))

    te_ref_seq = single_seq_fa(args.teref)
    ref = pysam.Fastafile(args.ref)

    meth_tbx = pysam.Tabixfile(args.meth)

    tmp_methdata = str(uuid4()) +'.tmp.methdata.tsv'

    with open(tmp_methdata, 'w') as meth_out:
        # header
        with gzip.open(args.meth, 'rt') as _:
            for line in _:
                assert line.startswith('chromosome')
                meth_out.write(line)
                break

        assert seg_chrom in meth_tbx.contigs

        for rec in meth_tbx.fetch(seg_chrom, seg_start, seg_end):
            meth_out.write(str(rec)+'\n')

    # index by read_name
    methdata = pd.read_csv(tmp_methdata, sep='\t', header=0, index_col=4)

    os.remove(tmp_methdata)

    reads = []
    if args.excl_ambig:
        reads = exclude_ambiguous_reads(args.bam, seg_chrom, seg_start, seg_end)
    else:
        reads = get_reads(args.bam, seg_chrom, seg_start, seg_end)

    reads = list(set(reads).intersection(set(methdata.index)))

    methdata = methdata.loc[reads]

    seg_reads = {}

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
                if index not in seg_reads:
                    seg_reads[index] = Read(index, cg_seg_start, llr)
                else:
                    seg_reads[index].add_cpg(cg_seg_start, llr)

    meth_table = dd(dict)
    sample = '.'.join(args.bam.split('.')[:-1])

    for name, read in seg_reads.items():
        for loc in read.llrs.keys():
            uuid = str(uuid4())
            meth_table[uuid]['loc'] = loc
            meth_table[uuid]['llr'] = read.llrs[loc]
            meth_table[uuid]['read'] = name
            meth_table[uuid]['sample'] = sample
            meth_table[uuid]['call'] = read.meth_calls[loc]

    meth_table = pd.DataFrame.from_dict(meth_table).T
    meth_table['loc'] = pd.to_numeric(meth_table['loc'])
    meth_table['llr'] = pd.to_numeric(meth_table['llr'])

    meth_table['orig_loc'] = meth_table['loc']
    meth_table['loc'] = ss.rankdata(meth_table['loc'], method='dense')

    coord_to_cpg = {}
    cpg_to_coord = {}
    for orig_loc, new_loc in zip(meth_table['orig_loc'], meth_table['loc']):
        coord_to_cpg[orig_loc] = new_loc
        cpg_to_coord[new_loc]  = orig_loc

    windowed_methfrac, meth_n = slide_window(meth_table, sample, width=int(args.slidingwindowsize), slide=int(args.slidingwindowstep))

    if len(windowed_methfrac) <= int(args.smoothwindowsize):
        logger.warning('too few sites after windowing: %s:%d-%d' % (seg_chrom, seg_start, seg_end))
        return [], []

    smoothed_methfrac = smooth(np.asarray(list(windowed_methfrac.values())), window_len=int(args.smoothwindowsize))

    coord_meth_pos = []

    cpg_meth_pos = list(windowed_methfrac.keys())

    for cpg in cpg_meth_pos:
        if seg_strand == '+':
            coord_meth_pos.append(cpg_to_coord[cpg])
        if seg_strand == '-':
            coord_meth_pos.append((seg_end-seg_start)-cpg_to_coord[cpg])

    # alignment to ref elt

    elt_seq = ref.fetch(seg_chrom, seg_start, seg_end)
    if seg_strand == '-':
        elt_seq = rc(elt_seq)

    te_ref_seq = te_ref_seq.upper()
    elt_seq = elt_seq.upper()

    s_ref = skseq.DNA(te_ref_seq)
    s_elt = skseq.DNA(elt_seq)

    aln_res = []

    try:
        if args.globalign:
            aln_res = skalign.global_pairwise_align_nucleotide(s_ref, s_elt)
        else:
            aln_res = skalign.local_pairwise_align_ssw(s_ref, s_elt)
    except IndexError: # scikit-bio throws this if no bases align  >:|
        logger.warning('no align on seg: %s:%d-%d' % (seg_chrom, seg_start, seg_end))
        return [], []
    
    coord_ref, coord_elt = aln_res[2]
    
    len_ref = coord_ref[1] - coord_ref[0]
    len_elt = coord_elt[1] - coord_elt[0]

    if len_ref / len(te_ref_seq) < float(args.lenfrac):
        logger.warning('ref align too short on seg: %s:%d-%d (%f)' % (seg_chrom, seg_start, seg_end, len_ref / len(te_ref_seq)))
        return [], []

    if len_elt / len(elt_seq) < float(args.lenfrac):
        logger.warning('elt align too short on seg: %s:%d-%d (%f)' % (seg_chrom, seg_start, seg_end, len_elt / len(elt_seq)))
        return [], []

    tab_msa = aln_res[0]

    elt_to_ref_coords = {}

    pos_ref = coord_ref[0]
    pos_elt = coord_elt[0]

    for pos in tab_msa.iter_positions():
        pos = list(pos)
        b_ref = pos[0]
        b_elt = pos[1]

        if '-' not in pos:
            elt_to_ref_coords[pos_elt] = pos_ref
            pos_ref += 1
            pos_elt += 1

        if b_elt == '-':
            pos_ref += 1

        if b_ref == '-':
            elt_to_ref_coords[pos_elt] = 'na'
            pos_elt += 1

    revised_coord_meth_pos = []
    meth_profile = []

    for pos, meth in zip(coord_meth_pos, smoothed_methfrac):
        if pos not in elt_to_ref_coords:
            continue

        revised_pos = elt_to_ref_coords[pos]

        if revised_pos != 'na':
            revised_coord_meth_pos.append(revised_pos)
            meth_profile.append(meth)


    return revised_coord_meth_pos, meth_profile


def main(args):
    fams = ['L1HS']
    te_ref_seq = single_seq_fa(args.teref)

    assert os.path.exists(args.ref + '.fai'), 'ref fasta must be indexed'

    if args.fams:
        fams = args.fams.split(',')

    logger.info('fams: %s' % ','.join(fams))

    segdata = pd.read_csv(args.segdata, sep='\t', header=0, index_col=0)

    sample = '.'.join(os.path.basename(args.bam).split('.')[:-1])

    if args.sample is not None:
        sample = args.sample

    assert sample+'_meth_calls' in segdata
    assert sample+'_unmeth_calls' in segdata

    useable = []

    for seg in segdata.index:
        use_seg = True

        if segdata['seg_name'].loc[seg] not in fams:
            use_seg = False
            continue

        if (segdata[sample+'_meth_calls'].loc[seg] + segdata[sample+'_unmeth_calls'].loc[seg]) < int(args.mincalls):
            use_seg = False
            continue

        if use_seg:
            useable.append(seg)

    if len(useable) > int(args.maxelts):
        useable = random.sample(useable, int(args.maxelts))

    segdata = segdata.loc[useable]

    logger.info('useable TEs: %d' % len(useable))


    pool = mp.Pool(processes=int(args.procs))

    results = []

    for seg in segdata.index:
        seg_chrom  = str(segdata['seg_chrom'].loc[seg])
        seg_start  = int(segdata['seg_start'].loc[seg])
        seg_end    = int(segdata['seg_end'].loc[seg])
        seg_strand = str(segdata['seg_strand'].loc[seg])
        seg_name   = str(segdata['seg_name'].loc[seg])

        res = pool.apply_async(get_meth_profile, [args, seg_chrom, seg_start, seg_end, seg_name, seg_strand])

        results.append(res)

    fig = plt.figure()
    gs = gridspec.GridSpec(2,1,height_ratios=[1,8])

    # cpg
    ax0 = plt.subplot(gs[0])

    cpg_start = 0
    cpg_end = len(te_ref_seq)

    if args.start:
        cpg_start = int(args.start)

    if args.end:
        cpg_end = int(args.end)

    assert cpg_start < cpg_end

    ax0.set_xlim((cpg_start, cpg_end))

    box = matplotlib.patches.Rectangle([0, 0], cpg_end-cpg_start, 1.0, edgecolor='#555555', facecolor='#cfcfcf', zorder=1)
    ax0.add_patch(box)

    if args.blocks:
        with open(args.blocks) as blocks:
            for line in blocks:
                b_start, b_end, b_name, b_col = line.strip().split()
                b_start = int(b_start)
                b_end = int(b_end)

                box = matplotlib.patches.Rectangle([b_start, 0], b_end-b_start, 1.0, edgecolor='#555555', facecolor=b_col, zorder=2)
                ax0.add_patch(box)


    cpg_locs = []

    for i in range(len(te_ref_seq)-1):
        if i >= cpg_start and i <= cpg_end:
            if te_ref_seq[i] == 'C' and te_ref_seq[i+1] == 'G':
                cpg_locs.append(i)

    ax0.vlines(cpg_locs, 0, 1, lw=1, colors=('#FF4500'), zorder=3)

    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.xaxis.set_ticks_position('top')

    # wiggles
    ax1 = plt.subplot(gs[1])

    out_res = [] # cache for --outelts

    for res in results:
        coord_meth_pos, meth_profile = res.get()

        if len(coord_meth_pos) == 0:
            continue

        out_res.append((coord_meth_pos, meth_profile))


    if len(out_res) < int(args.outelts):
        logger.info('available profiles (%d) is less than --outelts (%d)' % (len(out_res), int(args.outelts)))
        args.outelts = len(out_res)

    for coord_meth_pos, meth_profile in random.sample(out_res, int(args.outelts)):
        ax1.plot(coord_meth_pos, meth_profile, lw=float(args.linewidth), alpha=float(args.alpha), color=args.colour)

    ax1.set_xlim((0,len(te_ref_seq)))
    ax1.set_ylim((-0.05,1.05))

    ax1.set_xlim((cpg_start, cpg_end))

    fig.set_size_inches(16, 6)

    fn_base = sample + '.' + '_'.join(fams)

    if args.svg:
        plt.savefig('%s.composite.svg' % fn_base, bbox_inches='tight')
        logger.info('plotted to %s.composite.svg' % fn_base)
    else:
        plt.savefig('%s.composite.png' % fn_base, bbox_inches='tight')
        logger.info('plotted to %s.composite.png' % fn_base)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-s', '--segdata', required=True, help='segmeth output')
    parser.add_argument('-b', '--bam', required=True)
    parser.add_argument('-m', '--meth', required=True)
    parser.add_argument('-f', '--fams', default=None, help='families, comma delimited, default=L1HS')
    parser.add_argument('-r', '--ref', required=True, help='ref genome fasta')
    parser.add_argument('-t', '--teref', required=True, help='TE ref fasta')
    parser.add_argument('-p', '--procs', default=1, help='multiprocessing')
    parser.add_argument('-c', '--colour', default='#ff4500', help='colour (default: #ff4500')
    parser.add_argument('-a', '--alpha', default=0.3, help='alpha (default: 0.3)')
    parser.add_argument('-w', '--linewidth', default=1, help='line width (default: 1)')
    parser.add_argument('-l', '--lenfrac', default=0.95, help='fraction of TE length that must align (default 0.95)')
    parser.add_argument('--blocks', default=None, help='blocks to highlight (txt file with start, end, name, hex colour)')
    parser.add_argument('--start', default=None, help='start plotting at this base (default None)')
    parser.add_argument('--end', default=None, help='end plotting at this base (default None)')
    parser.add_argument('--sample', default=None, help='specify sample name (default = infer from .bam)')
    parser.add_argument('--mincalls', default=100, help='minimum call count to include elt (default = 100)')
    parser.add_argument('--maxelts', default=300, help='maximum elements, if > max random.sample() (default = 300)')
    parser.add_argument('--outelts', default=100, help='maximum output elements, if > max random.sample() (default = 300)')
    parser.add_argument('--slidingwindowsize', default=10, help='size of sliding window for meth frac (default 10)')
    parser.add_argument('--slidingwindowstep', default=1, help='step size for meth frac (default 1)')
    parser.add_argument('--smoothwindowsize', default=8, help='size of window for smoothing (default 8)')
    parser.add_argument('--globalign', action='store_true', default=False, help='experimental')
    parser.add_argument('--excl_ambig', action='store_true', default=False)
    parser.add_argument('--svg', action='store_true', default=False)


    args = parser.parse_args()
    main(args)
