#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)

# Number of strains (out of 133) with at least X many differences to all other 132 strains

thresholds <- 1:5

o <- foreach(threshold = thresholds, .combine = "rbind") %do% {
    tmp <- dat[N >= threshold][, .N, by=list(start, strainA)][N >= 100][, .N, by=list(start)]
    tmp[, "chromosome" := chromosome]
    tmp[, "threshold" := threshold]
    return(tmp)
}

res <- foreach(chromosome=1:16, .combine="rbind") %do% {
    dat <- fread(paste(chromosome, ".10000.10000.tab", sep=""))
    foreach(threshold = thresholds, .combine = "rbind") %do% {
        tmp <- dat[N >= threshold][, .N, by=list(start, strainA)][N == 132][, .N, by=list(start)]
        tmp[, "chromosome" := chromosome]
        tmp[, "threshold" := threshold]
        return(tmp)
    }
}
res[, stop := start + 10000 - 1]
setcolorder(res, c("chromosome", "start", "stop", "threshold", "N"))

res[, threshold := factor(threshold)]

res[threshold==5][order(-N)][1:10]
res[threshold==4][order(-N)][1:10]
res[threshold==3][order(-N)][1:10]
res[threshold==2][order(-N)][1:10]
res[threshold==1][order(-N)][1:10]



g <- ggplot(data=res, aes(x=start, y=N)) +
geom_point(aes(color=threshold)) +
facet_grid(chromosome~.) +
theme_few(12) +
labs(x="Position along chromosome", y="number of strains at least <threshold> sites away from all other 132")

ggsave(g, file="divergence.png", width=60, height=32, units="cm")