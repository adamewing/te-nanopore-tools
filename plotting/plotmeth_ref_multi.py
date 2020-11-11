#!/usr/bin/env python3

from __future__ import print_function

from collections import defaultdict as dd
from collections import Counter

import os
import pysam
import argparse

from operator import itemgetter

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
from matplotlib.patches import ConnectionPatch

from uuid import uuid4
import gzip


class Gene:
    def __init__(self, ensg, name):
        self.ensg = ensg
        self.name = name
        self.tx_start = None
        self.tx_end = None
        self.cds_start = None
        self.cds_end = None
        self.exons = []

    def add_exon(self, block):
        assert len(block) == 2
        assert block[0] < block[1]
        self.exons.append(block)
        self.exons = sorted(self.exons, key=itemgetter(0))

    def add_tx(self, block):
        assert len(block) == 2
        assert block[0] < block[1]
        if self.tx_start is None or self.tx_start > block[0]:
            self.tx_start = block[0]

        if self.tx_end is None or self.tx_end < block[1]:
            self.tx_end = block[1]

    def add_cds(self, block):
        assert len(block) == 2
        assert block[0] < block[1]

        if self.cds_start is None or self.cds_start > block[0]:
            self.cds_start = block[0]

        if self.cds_end is None or self.cds_end < block[1]:
            self.cds_end = block[1]

    def has_tx(self):
        return None not in (self.tx_start, self.tx_end)

    def has_cds(self):
        return None not in (self.cds_start, self.cds_end)

    def merge_exons(self):
        new_exons = []
        if len(self.exons) == 0:
            return

        last_block = self.exons[0]

        for block in self.exons[1:]:
            if min(block[1], last_block[1]) - max(block[0], last_block[0]) > 0: # overlap
                last_block = [min(block[0], last_block[0]), max(block[1], last_block[1])]

            else:
                new_exons.append(last_block)

            last_block = block

        new_exons.append(last_block)

        self.exons = new_exons


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


def exclude_ambiguous_reads(fn, chrom, start, end, min_mapq=10):
    reads = []

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch(chrom, start, end):
        p = read.get_reference_positions()
        if p[0] < start or p[-1] > end:
            if read.mapq >= min_mapq:
                reads.append(read.query_name)

    return reads

def get_ambiguous_reads(fn, chrom, start, end, min_mapq=10, w=50):
    reads = []

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch(chrom, start, end):
        p = read.get_reference_positions()
        if read.mapq < min_mapq or (p[0] > start-w and p[-1] < end+w):
            reads.append(read.query_name)

    return reads

def get_reads(fn, chrom, start, end, min_mapq=10):
    reads = []

    bam = pysam.AlignmentFile(fn)
    for read in bam.fetch(chrom, start, end):
        if read.mapq >= min_mapq:
            reads.append(read.query_name)

    return reads

def getmeth(args, bam, meth):

    # set up
    assert ':' in args.interval
    assert '-' in args.interval

    chrom, pos = args.interval.split(':')
    elt_start, elt_end = map(int, pos.split('-'))
    #window = int(args.window)

    bamname = '.'.join(os.path.basename(bam).split('.')[:-1])
    fn_prefix = bamname + '.' + '_'.join(args.interval.split(':')[:2])

    meth_tbx = pysam.Tabixfile(meth)

    tmp_methdata = fn_prefix+'.tmp.methdata.tsv'

    with open(tmp_methdata, 'w') as meth_out:
        # header
        with gzip.open(meth, 'rt') as _:
            for line in _:
                assert line.startswith('chromosome')
                meth_out.write(line)
                break

        assert chrom in meth_tbx.contigs

        for rec in meth_tbx.fetch(chrom, elt_start, elt_end):
            meth_out.write(str(rec)+'\n')

    # index by read_name
    methdata = pd.read_csv(tmp_methdata, sep='\t', header=0, index_col=4)

    if not args.keep_tmp_table:
        os.remove(tmp_methdata)

    # get list of relevant reads (exludes reads not anchored outside interval)
    reads = []
    if args.excl_ambig:
        reads = exclude_ambiguous_reads(bam, chrom, elt_start, elt_end)
    else:
        reads = get_reads(bam, chrom, elt_start, elt_end)
        
    reads = list(set(reads).intersection(set(methdata.index)))

    if args.unambig_highlight and args.highlight:
        h_coords = []
        for h in args.highlight.split(','):
            if ':' in h:
                h = h.split(':')[-1]
                
            h_coords += map(int, h.split('-'))

        h_coords.sort()

        h_start, h_end = h_coords[0], h_coords[-1]

        excl_reads = get_ambiguous_reads(bam, chrom, h_start, h_end)
        #print(excl_reads)

        new_reads = []
        for read in reads:
            if read not in excl_reads:
                new_reads.append(read)

        reads = new_reads


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

            cg_elt_start = cg_start - (elt_start)

            if cg_start >= elt_start and cg_start <= elt_end:
                #print (cg_start, cg_elt_start, llr, index)
                if index not in reads:
                    reads[index] = Read(index, cg_elt_start, llr)
                else:
                    reads[index].add_cpg(cg_elt_start, llr)

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


def mask_methfrac(data, cutoff=20):
    data = np.asarray(data)
    data = data > int(cutoff)

    segs = []

    in_seg = False
    seg_start = 0

    for i in range(len(data)):
        if data[i]:
            if in_seg:
                segs.append(list(range(seg_start, i)))

            in_seg = False

        else:
            if not in_seg:
                seg_start = i

            in_seg = True

    if in_seg:
        segs.append(list(range(seg_start, len(data))))

    return segs


def build_genes(gtf, chrom, start, end):
    genes = {}

    for line in gtf.fetch(chrom, start, end):

        chrom, source, feature, start, end, score, strand, frame, attribs = line.split('\t')

        block = [int(start), int(end)]

        attribs = attribs.strip()

        attr_dict = {}

        for attrib in attribs.split(';'):
            if attrib:
                key, val = attrib.strip().split()[:2]
                key = key.strip()
                val = val.strip().strip('"')
                attr_dict[key] = val

        if 'gene_id' not in attr_dict:
            continue

        if 'gene_name' not in attr_dict:
            continue

        ensg = attr_dict['gene_id']
        name = attr_dict['gene_name']

        if ensg not in genes:
            genes[ensg] = Gene(ensg, name)

        if feature == 'exon':
            genes[ensg].add_exon(block)

        if feature == 'CDS':
            genes[ensg].add_cds(block)

        if feature == 'transcript':
            genes[ensg].add_tx(block)

    return genes


def main(args):
    # set up
    assert ':' in args.interval
    assert '-' in args.interval

    chrom, pos= args.interval.split(':')
    elt_start, elt_end = map(int, pos.split('-'))
    
    # highlights

    h_start = []
    h_end = []
    h_cpg_start = []
    h_cpg_end = []

    h_colors = ['#3080ff','#ff4500']

    data = {}

    with open(args.data) as _:
        for line in _:
            bam, meth = line.strip().split()
            data[bam] = meth

    # table for plotting
    meth_table = dd(dict)

    sample_order = []

    reads = {}

    for bam, meth in data.items():
        bamname = '.'.join(os.path.basename(bam).split('.')[:-1])
        reads[bamname] = getmeth(args, bam, meth)

        for name, read in reads[bamname].items():
            for loc in read.llrs.keys():
                uuid = str(uuid4())
                meth_table[uuid]['loc'] = loc
                meth_table[uuid]['llr'] = read.llrs[loc]
                meth_table[uuid]['read'] = name
                meth_table[uuid]['sample'] = bamname
                meth_table[uuid]['call'] = read.meth_calls[loc]

                if bamname not in sample_order:
                    sample_order.append(bamname)

    meth_table = pd.DataFrame.from_dict(meth_table).T
    meth_table['loc'] = pd.to_numeric(meth_table['loc'])
    meth_table['llr'] = pd.to_numeric(meth_table['llr'])

    meth_table['orig_loc'] = meth_table['loc']
    meth_table['loc'] = ss.rankdata(meth_table['loc'], method='dense')

    coord_to_cpg = {}
    for orig_loc, new_loc in zip(meth_table['orig_loc'], meth_table['loc']):
        coord_to_cpg[orig_loc] = new_loc

    if args.highlight:
        for h in args.highlight.split(','):
            if ':' in h:
                h = h.split(':')[-1]
                
            h_s, h_e = map(int, h.split('-'))
            h_start.append(h_s)
            h_end.append(h_e)

            h_start[-1] -= elt_start
            h_end[-1] -= elt_start

            h_cpg_start.append(coord_to_cpg[min(meth_table['orig_loc'], key=lambda x:abs(x-h_start[-1]))])
            h_cpg_end.append(coord_to_cpg[min(meth_table['orig_loc'], key=lambda x:abs(x-h_end[-1]))])


    fig = plt.figure()
    gs = gridspec.GridSpec(4,1,height_ratios=[1,1,3,3])

    # plot genes
    ax0 = plt.subplot(gs[0])

    #ax0.set_xlim(elt_start, elt_end)
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.xaxis.set_ticks_position('top')

    gtf = None
    genes = [] 

    if args.gtf is not None:
        gtf = pysam.Tabixfile(args.gtf)
        genes = build_genes(gtf, chrom, elt_start, elt_end)

    exon_patches = []
    tx_lines = []

    genes_of_interest = []

    if args.genes is not None:
        genes_of_interest = args.genes.strip().split(',')

    i = 0
    for ensg in genes:
        if genes_of_interest:
            if genes[ensg].name not in genes_of_interest:
                continue

        if genes[ensg].has_tx():

            tx_lines.append(matplotlib.lines.Line2D([genes[ensg].tx_start-elt_start, genes[ensg].tx_end-elt_start], [0.4+i, 0.4+i], zorder=1))

            print('transcript: %d-%d %s' % (genes[ensg].tx_start, genes[ensg].tx_end, genes[ensg].name))

        #if genes[ensg].has_cds():
            #print('CDS: %d-%d %s' % (genes[ensg].cds_start, genes[ensg].cds_end, genes[ensg].name))

        genes[ensg].merge_exons()
        for exon_start, exon_end in genes[ensg].exons:
            exon_len = exon_end - exon_start

            exon_patches.append(matplotlib.patches.Rectangle([exon_start-elt_start, i], exon_len, 0.8, edgecolor='#777777', facecolor='#ff4500', zorder=2))


        blocks_str = ','.join(['%d-%d' % (s,e) for s, e in genes[ensg].exons])
        print('%s exons: %s' % (genes[ensg].name, blocks_str))

        i += 1

    if i < 3:
        i = 3

    ax0.set_ylim(0,i)
    ax0.set_yticks([])

    for p in exon_patches:
        ax0.add_patch(p)

    for tx in tx_lines:
        ax0.add_line(tx)


    # plot correspondence between genome space and cpg space
    ax1 = plt.subplot(gs[1])
    ax2 = ax1.twiny()

    ax1.set_ylim(0,10)
    ax1.set_yticklabels([])

    #ax1.set_xlim(elt_start, elt_end)

    x1 = []
    x2 = []

    step = int(args.topspacing)

    for i, x in enumerate(meth_table['orig_loc']):
        if i in (0, len(meth_table['orig_loc'])-1):
            x2.append(x)
            x1.append(coord_to_cpg[x])

        elif i % step == 0:
            x2.append(x)
            x1.append(coord_to_cpg[x])

    
    ax1.vlines(x1, 0, 1, color='#777777', zorder=1)
    ax2.vlines(x2, 9, 10, color='#777777', zorder=1)

    if args.highlight:
        for i in range(len(h_start)):
            orig_highlight_box = matplotlib.patches.Rectangle((h_start[i],9), h_end[i]-h_start[i], 1.0, lw=1, edgecolor='#777777', facecolor=h_colors[i], zorder=2)
            cpg_highlight_box = matplotlib.patches.Rectangle((h_cpg_start[i],0), h_cpg_end[i]-h_cpg_start[i], 1.0, lw=1, edgecolor='#777777', facecolor=h_colors[i], zorder=3)

            ax2.add_patch(orig_highlight_box)
            ax1.add_patch(cpg_highlight_box)

    for x1_x, x2_x in zip(x1, x2):
        link_end1 = (x1_x, 1)
        link_end2 = (x2_x, 9)

        l_col = '#777777'

        for i in range(len(h_start)):
            if x2_x >= h_start[i] and x2_x <= h_end[i]:
                l_col = h_colors[i]

        con = ConnectionPatch(xyA=link_end1, xyB=link_end2, coordsA="data", coordsB="data", axesA=ax1, axesB=ax2, color=l_col)
        ax2.add_artist(con)

    ax0.set_xlim(ax2.get_xlim()) # sync axes between orig coords and gtf plot
    ax2.set_xticks([])

    xt_labels = [str(int(t+elt_start)) for t in list(ax0.get_xticks())]
    ax0.set_xticklabels(xt_labels)


    # llr plot

    ax3 = plt.subplot(gs[2])

    ax3.axhline(y=2.5, c='k', linestyle='--',lw=1)
    ax3.axhline(y=0, c='#bbbbbb', linestyle='--',lw=1)
    ax3.axhline(y=-2.5, c='k', linestyle='--',lw=1)


    ax3 = sns.lineplot(x='loc', y='llr', hue='sample', data=meth_table)
    ax3.set_xlim(ax1.get_xlim())

    #ax.set_title('%s: CpG Methylation log-likelihood ratio' % args.interval)
    #ax.set_xlabel('')


    sample_color = {}
    for i, sample in enumerate(sample_order):
        sample_color[sample] = sns.color_palette(n_colors=len(sample_order))[i]


    # meth frac plot

    ax4 = plt.subplot(gs[3])

    order_stack = 1

    for sample in sample_order:
        windowed_methfrac, meth_n = slide_window(meth_table, sample, width=int(args.slidingwindowsize), slide=int(args.slidingwindowstep))

        smoothed_methfrac = smooth(np.asarray(list(windowed_methfrac.values())), window_len=int(args.smoothwindowsize))

        masked_segs = mask_methfrac(list(meth_n.values()), cutoff=args.maskcutoff)

        ax4.plot(list(windowed_methfrac.keys()), smoothed_methfrac, marker='', lw=4, color=sample_color[sample])

        order_stack += 1

        for seg in masked_segs:
            if len(seg) > 2:
                mf_seg = np.asarray(smoothed_methfrac)[seg]
                pos_seg = np.asarray(list(windowed_methfrac.keys()))[seg]
            
                ax4.plot(pos_seg, mf_seg, marker='', lw=4, color='#ffffff', alpha=0.8, zorder=order_stack)

                order_stack += 1

    ax4.set_xlim(ax1.get_xlim())
    ax4.set_ylim((-0.05,1.05))

    fn_prefix = '.'.join(args.data.split('.')[:-1]) + '.' + '_'.join(args.interval.split(':')[:2])

    if args.genes is not None:
        fn_prefix = '_'.join(args.genes.split(',')) + fn_prefix


    fig.set_size_inches(16, 8)
    if args.svg:
        plt.savefig('%s.multi.meth.svg' % fn_prefix, bbox_inches='tight')

    else:
        plt.savefig('%s.multi.meth.png' % fn_prefix, bbox_inches='tight')
    #meth_table.to_csv('%s.multi.meth.table.csv' % fn_prefix)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-d', '--data', required=True, help='text file with .bam filename and corresponding nanoport call-methylation tabix (whitespace-delimited)')
    parser.add_argument('-i', '--interval', required=True, help='chr:start-end')
    parser.add_argument('-g', '--gtf', default=None, help='genes or intervals to display in gtf format')
    parser.add_argument('-l', '--highlight', default=None, help='format: start-end')
    parser.add_argument('-s', '--slidingwindowsize', default=20, help='size of sliding window for meth frac (default=20)')
    parser.add_argument('-t', '--slidingwindowstep', default=2, help='step size for meth frac (default=2)')
    parser.add_argument('--smoothwindowsize', default=8, help='size of window for smoothing (default=8)')
    parser.add_argument('--maskcutoff', default=20, help='windowed read count masking cutoff (default=20)')
    parser.add_argument('--topspacing', default=10, help='spacing between links in top panel (default=10)')
    parser.add_argument('--genes', default=None, help='genes of interest (comma delimited)')
    parser.add_argument('--methcall_ymax', default=None)
    parser.add_argument('--methcall_xmin', default=None)
    parser.add_argument('--methcall_xmax', default=None)
    parser.add_argument('--keep_tmp_table', action='store_true', default=False)
    parser.add_argument('--excl_ambig', action='store_true', default=False)
    parser.add_argument('--unambig_highlight', action='store_true', default=False)
    parser.add_argument('--svg', action='store_true')

    args = parser.parse_args()
    main(args)
