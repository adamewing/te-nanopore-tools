# te-nanopore-tools
A cache of scripts for manipulating nanopore/nanopolish output aimed at studying TEs

## Plotting scripts

| script                | description |
|-----------------------|-------------|
| composite_meth.py     | Plots multiple CpG methylation profiles against a reference TE, generally a consensus TE. |
| plotmeth_ref_hap.py   | Plots allele-specific CpG methylation profiles given a phased .bam, tabix-indexed nanopolish call-methylation output, and an interval of interest |
| plotmeth_ref_multi.py | Plots one or more CpG metylation profiles given one or more phased .bams, each paired to a tabix-indexed nanopolish output, and an interval of interest |
| plotmeth_wg.py        | Plots the whole-genome methylation profile given output from wgmeth.py (see parsing tools) |
| segplot.py            | Makes strip plots or violin plots for specified samples given output from segmeth.py (see parsing tools) |