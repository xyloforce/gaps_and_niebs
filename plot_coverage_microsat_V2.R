library(ggplot2)
library(stringr)

source("~/setThemePoster.R")
args = commandArgs(trailingOnly = TRUE)
# args = c("data_human/niebs.bed",
#          "data_human/microsat_annotated.bed",
#          "data_human/genome.genome",
#          "data_human/niebs_annotated.bed",
#          "data_human/nieb_microsat.tsv",
#          "data_human/microsat_gap_niebs.bed")

print("loading datasets")
# we need to load niebs to select appropriate chromosomes
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

microsatxgap = read.delim(args[5],
                          header = FALSE,
                          col.names = c("chr", "start", "end",
                                        "name", "score", "strand"))
microsatxgap$len = microsatxgap$end - microsatxgap$start

uniq_overlap = read.delim(args[6],
                          header = TRUE,
                          col.names = c("chr", "start", "end", "name"))
uniq_overlap = uniq_overlap[uniq_overlap$chr %in% allowed_chrs, ]
uniq_overlap$len = uniq_overlap$end - uniq_overlap$start

microsat_uniq = read.delim(args[7],
                           header = TRUE,
                           col.names = c("chr",
                                         "start",
                                         "end",
                                         "name",
                                         "score",
                                         "strand"))
microsat_uniq$len = microsat_uniq$end - microsat_uniq$start

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

agg_microsat_uniq = aggregate(microsat_uniq$len,
                              by = list(microsat_uniq$name),
                              FUN = sum)

colnames(agg_microsat_uniq) = c("name", "total") # uniq bases covered by microsat category
cov_uniq_overlap = sum(uniq_overlap$len) /
  sum(genome$size)
agg_microsat_uniq$cov = agg_microsat_uniq$total /
  sum(uniq_overlap$len)

microsat_by_cat_coverage =
  data.frame("name" = agg_microsat$name,
             "cov" = agg_microsat$len / sum(genome$size))

microsatxgap[, c("name", "id")] = str_split(microsatxgap$name,
                                            "_",
                                            n = 2,
                                            simplify = TRUE)
microsatbytype = aggregate(microsatxgap$len,
                           by = list(microsatxgap$name,
                                     microsatxgap$id),
                           FUN = sum)
colnames(microsatbytype) = c("cover", "name", "total")

microsatbytype$cov = microsatbytype$total /
  agg_annot[match(microsatbytype$cover,
                  agg_annot$type), "total"]

print("merging dataframe")

final_df = microsatbytype
final_df = 
  rbind(final_df,
        data.frame("cover" = "unique part of overlap",
                   "name" = microsat_uniq$name,
                   "total" = microsat_uniq$len,
                   "cov" = agg_microsat_uniq[match(microsat_uniq$name,
                                                   agg_microsat_uniq$name),
                                             "cov"]))
final_df$rel = final_df$cov /
  microsat_by_cat_coverage[match(final_df$name,
                                 microsat_by_cat_coverage$name),
                           "cov"]

write.csv2(final_df, file = args[8])
print("plotting")
plot = ggplot(data = final_df[final_df$cover != "unique part of overlap", ],
              aes(x = name,
                  y = rel,
                  color = cover)) +
  geom_point(size = 3) +
  geom_point(data = final_df[final_df$cover == "unique part of overlap", ],
             shape = 21,
             fill = NA,
             size = 3) +
  coord_flip() +
  theme_poster +
  xlab("GC content") +
  ylab("Coverage") +
  scale_y_log10()

ggsave(args[9], width = 8, height = 10)
