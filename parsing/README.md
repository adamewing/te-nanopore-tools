# te-nanopore-tools
A cache of scripts for manipulating nanopore/nanopolish output aimed at studying TEs

## Parsing scripts

| script            | description |
|-------------------|-------------|
| segmeth.py        | Given a list of segments and one or more nanopore .bams, each paired to a tabix-indexed nanopolish output, compiles information on methylation calls for each segment. Output is used as input to other plotting tools |
| segmeth_hap.py    | Similar to segmeth.py but takes phased .bams as input and outputs allele-specific information |
| wgmeth_hap_dss.py | Generates output compatible with DSS for calling DMRs |
| wgmeth.py         | Given one or more nanopore .bams, each paired to a tabix-indexed nanopolish output, generate binned methylation data across the whole genome. Used as input to plotting scripts. |
| diffseg.py        | Adds columns for methylation fraction given output from segmeth.py |