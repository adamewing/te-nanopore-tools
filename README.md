# te-nanopore-tools
A cache of scripts and data for manipulating nanopore/nanopolish output aimed at studying TEs. See individual folders for further information.

### If possible I would suggest using TMNT instead of these scripts: https://github.com/adamewing/tmnt
TMNT includes updates and improvements to these tools, supports additional formats and base modfications, and is a more unified framework for exploring nanopore-based methylation data.


## Parsing scripts

| script            | description |
|-------------------|-------------|
| segmeth.py        | Given a list of segments and one or more nanopore .bams, each paired to a tabix-indexed nanopolish output, compiles information on methylation calls for each segment. Output is used as input to other plotting tools |
| segmeth_hap.py    | Similar to segmeth.py but takes phased .bams as input and outputs allele-specific information |
| wgmeth_hap_dss.py | Generates output compatible with DSS for calling DMRs |
| wgmeth.py         | Given one or more nanopore .bams, each paired to a tabix-indexed nanopolish output, generate binned methylation data across the whole genome. Used as input to plotting scripts. |
| diffseg.py        | Adds columns for methylation fraction given output from segmeth.py |


## Plotting scripts

| script                | description |
|-----------------------|-------------|
| composite_meth.py     | Plots multiple CpG methylation profiles against a reference TE, generally a consensus TE. |
| plotmeth_ref_hap.py   | Plots allele-specific CpG methylation profiles given a phased .bam, tabix-indexed nanopolish call-methylation output, and an interval of interest |
| plotmeth_ref_multi.py | Plots one or more CpG metylation profiles given one or more phased .bams, each paired to a tabix-indexed nanopolish output, and an interval of interest |
| plotmeth_wg.py        | Plots the whole-genome methylation profile given output from wgmeth.py (see parsing tools) |
| segplot.py            | Makes strip plots or violin plots for specified samples given output from segmeth.py (see parsing tools) |


Citation: Adam D. Ewing, Nathan Smits, Francisco J. Sanchez-Luque, Sandra R. Richardson, Seth W. Cheetham, Geoffrey J. Faulkner. Nanopore Sequencing Enables Comprehensive Transposable Element Epigenomic Profiling. 2020. Molecular Cell, Online ahead of print: https://doi.org/10.1016/j.molcel.2020.10.024
