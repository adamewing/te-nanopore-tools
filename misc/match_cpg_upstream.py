#!/usr/bin/env python

import argparse
import pysam

def count_cpg(ref, chrom, start, end):
    seq = ref.fetch(chrom, start, end)
    return seq.count('CG')


def main(args):
    ref = pysam.Fastafile(args.ref)

    with open(args.bed) as bed:
        for line in bed:
            chrom, start, end, name, strand = line.strip().split()
            start = int(start)
            end = int(end)

            upstream_start = start-int(args.minlen)
            upstream_end = start

            if strand == '-':
                upstream_start = end
                upstream_end = end+int(args.minlen)


            target_count = count_cpg(ref, chrom, start, end)
            upstream_count = count_cpg(ref, chrom, upstream_start, upstream_end)

            if target_count < int(args.mincpgcount):
                continue

            tries = 0

            while upstream_count < target_count:
                if tries > int(args.maxtries):
                    break

                tries += 1

                if strand == '+':
                    upstream_start -= int(args.stepsize)

                if strand == '-':
                    upstream_end += int(args.stepsize)

                upstream_count = count_cpg(ref, chrom, upstream_start, upstream_end)


            if upstream_count >= target_count:
                print('\t'.join((chrom, str(upstream_start), str(upstream_end), name+'_UPSTREAM', strand, str(target_count), str(upstream_count))))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='giant bucket')
    parser.add_argument('-b', '--bed', required=True)
    parser.add_argument('-r', '--ref', required=True)
    parser.add_argument('--minlen', default=1000)
    parser.add_argument('--stepsize', default=10)
    parser.add_argument('--maxtries', default=10000)
    parser.add_argument('--mincpgcount', default=10)
    args = parser.parse_args()
    main(args)