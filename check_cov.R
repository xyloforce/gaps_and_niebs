args = commandArgs(trailingOnly = TRUE)

bed_colnames = c("chr", "start", "end", "name")

niebs_annot = read.delim(args[4],
                         header = FALSE,
                         col.names = bed_colnames)
niebs_annot$len = niebs_annot$end - niebs_annot$start

gaps = read.delim(args[1],
                  header = FALSE,
                  col.names = bed_colnames)
gaps = gaps[gaps$chr %in% niebs_annot$chr, ]
gaps$len = gaps$end - gaps$start

blocks = read.delim(args[2],
                    header = FALSE,
                    col.names = bed_colnames[0:3])
blocks = blocks[blocks$chr %in% niebs_annot$chr, ]
blocks$len = blocks$end - blocks$start

genome = read.delim(args[3],
                    header = FALSE,
                    col.names = c("chr", "size"))
genome = genome[genome$chr %in% niebs_annot$chr, ]
genome_size = sum(genome$size)

agg_annot = aggregate(niebs_annot$len,
                      by = list(niebs_annot$name),
                      FUN = sum)
colnames(agg_annot) = c("type", "total") # bases covered by nieb category
# agg_annot$cov = agg_annot$total / sum(genome$size)

cov_gaps = sum(gaps$len) / genome_size
cov_blocks = sum(blocks$len) / genome_size
cov_niebs = sum(niebs_annot$len) / genome_size
agg_annot$rel = agg_annot$total / sum(niebs_annot$len)
# agg_annot[agg_annot$type == "blocs",
#           "cov"] = agg_annot[agg_annot$type == "blocs", "total"] /
#   sum(blocks$len)
# agg_annot[agg_annot$type == "gaps",
#           "cov"] = agg_annot[agg_annot$type == "gaps", "total"] /
#   sum(gaps$len)

# agg_annot$rel = agg_annot$cov / cov_niebs
# agg_annot$cov = agg_annot$total / genome_size
print(agg_annot)