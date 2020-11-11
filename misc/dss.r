#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) == 1) {

    library(DSS)
    require(bsseq)

    li1_fn <- paste(args[1], "phase_0.DSS.sorted.txt", sep='.')
    li2_fn <- paste(args[1], "phase_1.DSS.sorted.txt", sep='.')

    dmr_fn <- paste(args[1], "DMRs.txt", sep='.')
    dml_fn <- paste(args[1], "DMLs.txt", sep='.')

    li1 = read.table(li1_fn, header=TRUE)
    li2 = read.table(li2_fn, header=TRUE)
    BSobj = makeBSseqData(list(li1, li2), c('li1', 'li2'))

    dmlTest = DMLtest(BSobj, group1=c('li1'), group2=c('li2'), smoothing=TRUE)
    dmrs = callDMR(dmlTest, p.threshold=0.05)

    write.table(dmrs, file=dmr_fn, quote=F, sep='\t')

    dmls = callDML(dmlTest, p.threshold=1)
    write.table(dmls, file=dml_fn, quote=F, sep='\t')

}

if(length(args) != 1) {
    cat("usage: dss.r <sample name>\n")
}
