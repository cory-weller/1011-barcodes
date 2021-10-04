#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(parallel)
library(doMC)
registerDoMC(cores=4)


args <- commandArgs(trailingOnly=TRUE)

chromosome <- args[1]
windowSize <- as.numeric(args[2])
windowStep <- as.numeric(args[3])


loadVCF <- function(chromosome) {
    fread(cmd=paste("zcat genomes/chromosome", chromosome, ".vcf.gz", sep=""))
}

# Load full VCF
vcf <- loadVCF(chromosome)

# Load header names
headerText <- colnames(fread("header.txt"))

# Reformat VCF header
setnames(vcf, headerText)
setnames(vcf, "#CHROM", "CHROM")

# Get list of haploid strains
haploids <- colnames(fread("haploids.txt"))

# Subset info columns (1-9) plus genotypes of haploid strains
infoCols <- colnames(vcf)[1:9]
subsetCols <- c(infoCols, haploids)
vcf <- vcf[, subsetCols, with=FALSE]
setkey(vcf, POS)

# Get table strain pairs (no repeats)
strain_pairs <- CJ("strainA" = haploids, "strainB" = haploids)[strainA != strainB]

# Get list of windows start positions
windowStarts <- seq(1, max(vcf$POS), windowStep)

# Filter to only include biallelic SNPs? (remove indels or polyalelic sites?)
# TBD


# outer loop: for every window
o <- rbindlist(foreach(windowStart = windowStarts) %do% {
    windowStop <- windowStart + windowSize - 1
    vcf.sub <- vcf[POS %between% c(windowStart, windowStop)]

    # inner loop: for every strain pair
    rbindlist(foreach(i = 1:nrow(strain_pairs)) %dopar% {

        # extract strain names from the strain_pairs table
        strainA = strain_pairs[i, strainA]
        strainB = strain_pairs[i, strainB]
        nDifferences <- nrow(vcf.sub[(get(strainA) == "0/0" & get(strainB) == "1/1") | (get(strainA) == "1/1" & get(strainB) == "0/0")])
        data.table("CHROM" = chromosome,
        "start" = windowStart,
        "stop" = windowStop, 
        "strainA" = strainA,
        "strainB" = strainB,
        "N" = nDifferences)
    })
})

fwrite(o, file=paste(chromosome, windowSize, windowStep, "tab", sep="."), quote=F, row.names=F, col.names=T, sep="\t")