#!/usr/bin/env python

import csv
import argparse
import pysam
import numpy as np
import pandas as pd

from collections import defaultdict as dd
from bx.intervals.intersection import Intersecter, Interval

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def avgmap(maptabix, chrom, start, end):
    ''' return average mappability across chrom:start-end region; maptabix = pysam.Tabixfile'''
    scores = []

    if None in (start, end): return None

    if chrom in maptabix.contigs:
        for rec in maptabix.fetch(chrom, int(start), int(end)):
            mchrom, mstart, mend, mscore = rec.strip().split()
            mstart, mend = int(mstart), int(mend)
            mscore = float(mscore)

            while mstart < mend and mstart:
                mstart += 1
                if mstart >= int(start) and mstart <= int(end):
                    scores.append(mscore)

        if len(scores) > 0:
            return sum(scores) / float(len(scores))
        else:
            return 0.0
    else:
        return 0.0


def main(args):

    tldr_forest    = dd(Intersecter)
    tebreak_forest = dd(Intersecter)
    map_tbx = pysam.Tabixfile(args.mappability)

    with open(args.tldr) as fh:
        csv_reader = csv.DictReader(fh, delimiter='\t')
        for rec in csv_reader:
            tldr_forest[rec['Chrom']].add_interval(Interval(int(rec['Start']), int(rec['End']), value=rec))

    with open(args.tebreak) as fh:
        csv_reader = csv.DictReader(fh, delimiter='\t')
        for rec in csv_reader:
            tebreak_forest[rec['Chromosome']].add_interval(Interval(int(rec['Left_Extreme']), int(rec['Right_Extreme']), value=rec))

    matches = []
    tldr_nomatch = []
    tebreak_nomatch = []

    tldr_count = 0
    tldr_nr_count = 0
    tldr_nr_2_count = 0
    tldr_pass_count = 0
    tldr_nr_pass_count = 0

    tebreak_count = 0
    tebreak_nr_count = 0
    tebreak_nr_2_count = 0
    tebreak_pass_count = 0
    tebreak_nr_pass_count = 0


    with open(args.tldr) as fh:
        csv_reader = csv.DictReader(fh, delimiter='\t')
        for rec in csv_reader:
            matched = False

            for tb in tebreak_forest[rec['Chrom']].find(int(rec['Start']), int(rec['End'])):
                matches.append([rec, tb.value])
                matched = True

            if not matched and rec['Filter'] == 'PASS':
                tldr_nomatch.append(rec)

            tldr_count += 1

            if rec['NonRef'] != 'NA':
                tldr_nr_count += 1

            if len(rec['NonRef'].split(',')) > 1:
                tldr_nr_2_count += 1

            if rec['Filter'] == 'PASS':
                if rec['NonRef'] != 'NA':
                    tldr_nr_pass_count += 1
                tldr_pass_count += 1


    with open(args.tebreak) as fh:
        csv_reader = csv.DictReader(fh, delimiter='\t')
        for rec in csv_reader:
            matched = False

            for tldr in tldr_forest[rec['Chromosome']].find(int(rec['Left_Extreme']), int(rec['Right_Extreme'])):
                matched = True

            if not matched and rec['Filter'] == 'PASS':
                tebreak_nomatch.append(rec)

            tebreak_count += 1

            if rec['NonRef'] != 'NA':
                tebreak_nr_count += 1

            if len(rec['NonRef'].split(',')) > 1:
                tebreak_nr_2_count += 1

            if rec['Filter'] == 'PASS':
                if rec['NonRef'] != 'NA':
                    tebreak_nr_pass_count += 1
                tebreak_pass_count += 1

    nr_matches = [m for m in matches if (m[0]['NonRef'] != 'NA' or m[1]['NonRef'] != 'NA')]
    nr_2_matches = [m for m in matches if (len(m[0]['NonRef'].split(',')) > 1 or len(m[1]['NonRef'].split(',')) > 1)]

    pass_AND_matches = [m for m in matches if (m[0]['Filter'] == 'PASS' and m[1]['Filter'] == 'PASS')]
    pass_OR_matches = [m for m in matches if (m[0]['Filter'] == 'PASS' or m[1]['Filter'] == 'PASS')]
    
    tebreak_nomatch_mapscores = []
    tldr_nomatch_mapscores = []

    match_mapscores = []

    mapscore_table = dd(dict)

    for rec in tldr_nomatch:
        chrom = rec['Chrom']
        start = int(rec['Start'])
        end   = int(rec['End'])

        flank = (200-(end-start))/2

        score = avgmap(map_tbx, chrom, start-flank, end+flank)

        if score == 0.5:  # usually due to alt haplotype in assembly rather than actual mappability
            continue

        tldr_nomatch_mapscores.append(score)

        mapscore_table[rec['UUID']]['category'] = 'tldr'
        mapscore_table[rec['UUID']]['mapscore'] = score


    for rec in tebreak_nomatch:
        chrom = rec['Chromosome']
        start = int(rec['5_Prime_End'])
        end   = int(rec['3_Prime_End'])

        if start > end:
            start, end = end, start

        flank = (200-(end-start))/2

        score = avgmap(map_tbx, chrom, start-flank, end+flank)

        if score == 0.5:  # usually due to alt haplotype in assembly rather than actual mappability
            continue

        tebreak_nomatch_mapscores.append(score)

        mapscore_table[rec['UUID']]['category'] = 'tebreak'
        mapscore_table[rec['UUID']]['mapscore'] = score

    for m in matches:
        if m[0]['Filter'] != 'PASS' or m[1]['Filter'] != 'PASS':
            continue

        rec = m[0]

        chrom = rec['Chrom']
        start = int(rec['Start'])
        end   = int(rec['End'])

        flank = (200-(end-start))/2

        score = avgmap(map_tbx, chrom, start-flank, end+flank)

        if score == 0.5:  # usually due to alt haplotype in assembly rather than actual mappability
            continue
            
        match_mapscores.append(score)

        mapscore_table[rec['UUID']]['category'] = 'both'
        mapscore_table[rec['UUID']]['mapscore'] = score

    mapscore_table = pd.DataFrame.from_dict(mapscore_table).T
    mapscore_table = pd.DataFrame(mapscore_table.to_dict())

    sns_plot = sns.stripplot(x='category', y='mapscore', data=mapscore_table, size=2, dodge=True, jitter=True)

    fig = sns_plot.figure
    fig.set_size_inches(4, 8)
    fig.savefig('mapscore_cmp.svg', bbox_inches='tight')


    print('tldr count:', tldr_count)
    print('tldr nr count:', tldr_nr_count)
    print('tldr nr(2) count:', tldr_nr_2_count)
    print('tldr pass count:', tldr_pass_count)
    print('tldr nr pass count:', tldr_nr_pass_count)
    print()

    print('tebreak count:', tebreak_count)
    print('tebreak nr count:', tebreak_nr_count)
    print('tebreak nr(2) count:', tebreak_nr_2_count)
    print('tebreak pass count:', tebreak_pass_count)
    print('tebreak nr pass count:', tebreak_nr_pass_count)
    print()

    print('overall matches:', len(matches))
    print('NonRef matches:', len(nr_matches))
    print('NonRef 2 or more studies matches:', len(nr_2_matches))
    print('PASS AND matches:', len(pass_AND_matches))
    print('PASS OR matches:', len(pass_OR_matches))
    print()

    print('tldr PASS nomatch:', len(tldr_nomatch))
    print('tebreak PASS nomatch:', len(tebreak_nomatch))
    print('tldr PASS nomatch avgmap:', np.mean(tldr_nomatch_mapscores))
    print('tebreak PASS nomatch avgmap:', np.mean(tebreak_nomatch_mapscores))
    print('matched PASS match avgmap:', np.mean(match_mapscores))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-t', '--tldr', required=True, help='tldr table')
    parser.add_argument('-b', '--tebreak', required=True, help='tebreak table')
    parser.add_argument('-m', '--mappability', required=True)
    args = parser.parse_args()
    main(args)
