library(ggplot2)
library(stringr)

source("~/setThemePoster.R")
args = commandArgs(trailingOnly = TRUE)

print("loading datasets")
niebs = read.delim(args[1],
                   header = FALSE,
                   col.names = c("chr", "start", "end", "name"))
niebs$len = niebs$end - niebs$start
allowed_chrs = unique(niebs$chr)

microsat = read.delim(args[2],
                      header = FALSE,
                      col.names = c("chr", "start", "end",
                                    "name", "score", "strand"))
microsat$len = microsat$end - microsat$start
microsat = microsat[microsat$chr %in% allowed_chrs, ]

genome = read.delim(args[3],
                    header = FALSE,
                    col.names = c("chr", "size"))
genome = genome[genome$chr %in% allowed_chrs, ]

niebs_annot = read.delim(args[4],
                         header = FALSE,
                         col.names = c("chr", "start", "end", "name"))
niebs_annot$len = niebs_annot$end - niebs_annot$start

niebXmicrosat = read.delim(args[5],
                           header = FALSE,
                           col.names = c("chr", "start", "end",
                                         "name", "score", "strand"))
niebXmicrosat$len = niebXmicrosat$end - niebXmicrosat$start

microsatxgap = read.delim(args[6],
                          header = FALSE,
                          col.names = c("chr", "start", "end",
                                        "name", "score", "strand"))
microsatxgap$len = microsatxgap$end - microsatxgap$start

print("calculating coverage")

agg_microsat = aggregate(microsat$len,
                         by = list(microsat$name),
                         FUN = sum)
colnames(agg_microsat) = c("name", "len") # bases covered by microsat category

agg_annot = aggregate(niebs_annot$len,
                      by = list(niebs_annot$name),
                      FUN = sum)
colnames(agg_annot) = c("type", "total") # bases covered by nieb category
agg_annot$cov = agg_annot$total / sum(genome$size)

microsat_by_cat_coverage =
  data.frame("name" = agg_microsat$name,
             "cov" = agg_microsat$len / sum(genome$size))


agg_niebXmicro = aggregate(niebXmicrosat$len,
                           by = list(niebXmicrosat$name),
                           FUN = sum)
colnames(agg_niebXmicro) = c("name", "len")

agg_niebXmicro$cov = agg_niebXmicro$len / sum(niebs$len)


agg_niebXmicro$genomecov =
  microsat_by_cat_coverage$cov
agg_niebXmicro = rbind(agg_niebXmicro,
                       data.frame(name = "total",
                                  len = sum(agg_niebXmicro$len),
                                  cov = sum(agg_niebXmicro$len) /
                                    sum(niebs$len),
                                  genomecov = sum(agg_microsat$len) /
                                    sum(genome$size)))
agg_niebXmicro$cov = agg_niebXmicro$cov * 100
agg_niebXmicro$genomecov = agg_niebXmicro$genomecov * 100
write.csv2(agg_niebXmicro,
           file = args[7])