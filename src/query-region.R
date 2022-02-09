#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(doMC)
library(ggplot2)
registerDoMC(cores=4)

options(scipen=999)


args <- commandArgs(trailingOnly = TRUE)

headerText <- colnames(fread('data/header.txt'))

chromosome <- args[1]
start <- args[2]
stop <- args[3]


importGenotypes <- function(chromosome, start, stop) {
    vcfFilename <- paste('data/genomes/chromosome', chromosome, '.vcf.gz', sep='')
    tabixCommand <- paste('tabix ', vcfFilename, ' chromosome',chromosome, ':', start, '-', stop, sep='')
    genotypes <- fread(cmd=tabixCommand)
    if(nrow(genotypes) == 0) {
        return(NULL)
    }
    setnames(genotypes, headerText)
    strainNames <- colnames(genotypes)[10:ncol(genotypes)]
    genotypes <- genotypes[, (strainNames) := lapply(.SD, function(x) factor(x, levels=c('0/0', './.', '1/1'))), .SDcols=strainNames][]
    genotypes <- genotypes[, (strainNames) := lapply(.SD, function(x) as.numeric(x)), .SDcols=strainNames]
    replaceValues(genotypes, 2, NA)
    replaceValues(genotypes, 1, 0)
    replaceValues(genotypes, 3, 1)
    return(genotypes[])
}


matingTypes <- fread('data/mating-types.tsv')
matingTypes <- matingTypes[Mating %in% c('a','b')]

#if (any("1/1" == DT[,10:length(colnames(DT)), with=FALSE])) {
    # convert to factor, then to numeric, such that "0/0" is now 1, "1/1" is now 3

fillNA <- function(DT, x) {
    # replaces all NA values with x in data.table DT
    for (j in seq_len(ncol(DT)))
        set(DT,which(is.na(DT[[j]])),j,x)
}

replaceValues <- function(DT, x1, x2) {
    # replaces all NA values with x in data.table DT
    for (j in seq_len(ncol(DT)))
        set(DT, which(DT[[j]] == x1), j, x2)
}


testRegion <- function(matingTypes, chromosome, start, stop) {
    DT <- importGenotypes(chromosome, start, stop)
    if(nrow(DT) == 0) {
        print('no variants within region')
    }
    setkey(DT, POS)
    print(paste('testing region ', chromosome, ':', start, '-', stop, sep=''))
    o <- foreach(matingType = c('a','b'), .combine='rbind') %do% {
        print(paste('mating type: ', matingType, sep=''))
        starting.strain <- 'S288C'
        strains <- matingTypes[Mating == matingType, Standardized_name]
        remaining <- strains[strains != starting.strain]
        chosen <- DT[, c('#CHROM', 'POS'), with=FALSE]
        chosen[, 'S288C' := 0]
        chosen[, c('#CHROM', 'POS') := NULL]
        #chosen[, nRef := 1]
        #chosen[, nAlt := 0]
        # while there are still candidates remaining to be checked
        while(length(remaining) > 0) {
            # randomly take one of the remaining strains and extract its genotypes
            candidate <- sample(remaining, size=1)
            candidate.genotypes <- DT[, get(candidate)]
            # if candidate isn't all NAs...
            if(! all(is.na(candidate.genotypes))) {
                # if not any chosen strain for which all genotypes are equal to candidate
                if(!any(chosen[, lapply(.SD, function(x) all(candidate.genotypes==x, na.rm=T))])) {
                    # add candidate
                    chosen[, eval(quote(candidate)) := candidate.genotypes]
                } 
            }
            # eliminate candidate from remaining regardless of whether it was added to the pot
            remaining <- remaining[remaining != candidate]
        }
        n.strains <- ncol(chosen)
    chosenStrains <- paste(colnames(chosen)[! colnames(chosen) %in% c('#CHROM','POS', 'nRef','nAlt')], collapse=',')
    return(data.table(chromosome, start, stop, matingType, n.strains, chosenStrains))
    }
    return(o)
}

output <- testRegion(matingTypes, chromosome, start, stop)

fwrite(output, file=paste('reports/', 'chr', chromosome, '-', start, '-', stop, '.tsv', sep=''), quote=F, col.names=T, row.names=F, sep="\t")