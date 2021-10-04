#!/usr/bin/env R

library(cluster)
library(data.table)
library(ggplot2)
library(foreach)

dat <- fread('1.500.500.tab')

mating_types <- fread('matingType.tsv')

geneInfo <- fread(cmd="zcat SacCer3_SGD_Genes.gz")
expressionInfo <- fread(cmd="zcat SacCer3_SGD_regulation.gz")


matA <- mating_types[Mating == "a", Standardized_name]

matAlpha <- mating_types[Mating == "b", Standardized_name]


window_start <- 60001



dat.sub <- dat[start == window_start][, c("strainA", "strainB", "N")]
dm <- dcast(dat.sub, strainA~strainB, value.var="N", fill=0)[, !1]


o <- foreach(window_start = unique(dat$start), .combine="rbind") %do% {

    dat.subA <- dat[strainA %in% matA & strainB %in% matA & start == window_start][, c("strainA", "strainB", "N")]
    dmA <- dcast(dat.subA, strainA~strainB, value.var="N", fill=0)[, !1]

    dat.subAlpha <- dat[strainA %in% matAlpha & strainB %in% matAlpha & start == window_start][, c("strainA", "strainB", "N")]
    dmAlpha <- dcast(dat.subAlpha, strainA~strainB, value.var="N", fill=0)[, !1]


    o.tmp <- foreach(k=2:50, .combine="rbind") %do% {
        aswA <- pam(dmA, k)$silinfo$avg.width
        aswAlpha <- pam(dmAlpha, k)$silinfo$avg.width
        out.A <- data.table("k"=k, "asw" = aswA, "mating_type" = "a", "window_start" = window_start)
        out.Alpha <- data.table("k"=k, "asw" = aswAlpha, "mating_type" = "alpha", "window_start" = window_start)
        rbindlist(list(out.A, out.Alpha))
    }
    o.tmp[, list("k"=which.max(asw)), by=list(mating_type, window_start)]
}

ggplot(o, aes(x=k, y=asw)) + geom_point() + facet_wrap(~window_start)


data.table(strain=dm$strainA, cluster=pam(dm, k=50, cluster.only=TRUE))

asw <- numeric(100)

for (k in 1:100)
    asw[k] <- pam(dm, k) $ silinfo $ avg.width

k.best <- which.max(asw)