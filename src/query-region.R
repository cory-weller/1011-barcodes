#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(doMC)
library(ggplot2)
registerDoMC(cores=4)

args <- commandArgs(trailingOnly = TRUE)

headerText <- colnames(fread('header.txt'))

chromosome <- args[1]
start <- args[2]
stop <- args[3]




importGenotypes <- function(chromosome, start, stop) {
    vcfFilename <- paste('genomes/chromosome', chromosome, '.vcf.gz', sep='')
    tabixCommand <- paste('tabix ', vcfFilename, ' chromosome',chromosome, ':', start, '-', stop, sep='')
    output <- fread(cmd=tabixCommand)
    setnames(output, headerText)
    return(output[])
}


matingTypes <- fread('mating-types.tsv')
matingTypes <- matingTypes[Mating %in% c('a','b')]

#if (any("1/1" == DT[,10:length(colnames(DT)), with=FALSE])) {
    # convert to factor, then to numeric, such that "0/0" is now 1, "1/1" is now 3

fillNA <- function(DT, x) {
    # replaces all NA values with x in data.table DT
    for (j in seq_len(ncol(DT)))
        set(DT,which(is.na(DT[[j]])),j,x)
}


# ref = 0
# no data = 1
# alt = 2
# anything else = 9


# a.genotypes[, nRef := apply(.SD, 1, function(x) sum(x=="0/0", na.rm=TRUE)), .SDcols=cols]
# a.genotypes[, nAlt := apply(.SD, 1, function(x) sum(x=="1/1", na.rm=TRUE)), .SDcols=cols]
# a.genotypes <- a.genotypes[nRef != 0 & nAlt != 0]

testRegion <- function(DT, matingTypes, chromosome, start, stop) {
    print(paste('testing region ', chromosome, ':', start, '-', stop, sep=''))
    o <- foreach(matingType = c('a','b'), .combine='rbind') %do% {
        print(paste('mating type: ', matingType, sep=''))
        starting.strain <- ifelse(matingType=='a', 'AHF', ifelse(matingType=='b', 'ABE', NA))
        strains <- matingTypes[Mating == matingType, Standardized_name]
        foreach(divergence.threshold = 1:3, .combine='rbind') %do% {
            print(paste('divergence threshold: ', divergence.threshold, sep=''))
            foreach(iteration = 1:20, .combine='rbind') %dopar% {
                remaining <- strains[strains != starting.strain]
                chosen <- DT[, c('#CHROM', 'POS', starting.strain), with=FALSE]
                chosen[, nRef := fifelse(get(starting.strain) == 0, 1, 0)]
                chosen[, nAlt := fifelse(get(starting.strain) == 2, 1, 0)]
                # while there are still candidates remaining to be checked
                while(length(remaining) > 0) {
                    # randomly take one of the remaining strains and extract its genotypes
                    candidate <- sample(remaining, size=1)
                    candidate.genotypes <- DT[, get(candidate)]
                    # if candidate differs at any sites that are currently fixed in the population
                    #if(any(chosen$nRef == 0 & candidate.genotypes == 0) | any(chosen$nAlt == 0 & candidate.genotypes == 2)) {
                    if(sum(chosen$nRef == 0 & candidate.genotypes == 0) + sum(chosen$nAlt == 0 & candidate.genotypes == 2) >= divergence.threshold  ) {
                        # add candidate
                        chosen[, eval(quote(candidate)) := candidate.genotypes]
                        # pdate ref and alt counts
                        chosen[get(candidate) == 0, nRef := nRef + 1]
                        chosen[get(candidate) == 2, nAlt := nAlt + 1]
                    } 
                    # eliminate candidate from remaining regardless of whether it was added to the pot
                    remaining <- remaining[remaining != candidate]
                }
                n.strains <- ncol(chosen)-4
            chosenStrains <- paste(colnames(chosen)[! colnames(chosen) %in% c('#CHROM','POS', 'nRef','nAlt')], collapse=',')
            return(data.table(chromosome, start, stop, matingType, iteration, divergence.threshold, n.strains, chosenStrains))
            }
        }
    }
    g <- ggplot(o, aes(x=factor(divergence.threshold), y=n.strains, color=matingType)) +
        geom_boxplot() +
        labs(x='# of unique SNPs required', y='# of strains uniquely identifiable at (x) SNPs')
        ggsave(g, file=paste('chr', chromosome, '-', start, '-', stop, '-unique-strains.png', sep=''))
    return(o)
}

genotypes <- importGenotypes(chromosome, start, stop)
strainNames <- colnames(genotypes)[10:ncol(genotypes)]
genotypes <- genotypes[, (strainNames) := lapply(.SD, function(x) factor(x, levels=c('0/0', './.', '1/1'))), .SDcols=strainNames][]
genotypes <- genotypes[, (strainNames) := lapply(.SD, function(x) as.numeric(x)), .SDcols=strainNames]
fillNA(genotypes, 9)

output <- testRegion(genotypes, matingTypes, chromosome, start, stop)

fwrite(output, file=paste('chr', chromosome, '-', start, '-', stop, '.tsv', sep=''), quote=F, col.names=T, row.names=F, sep="\t")
