library(ggplot2)

# args = c("data_human/blocks_total.bed",
#          "data_human/cleaned_uniq.bed")
source("~/setThemePoster.R")
args = commandArgs(trailingOnly = TRUE)

bed_colnames = c("chr", "start", "end")
blocks = read.delim(args[1],
                    header = FALSE,
                    col.names = bed_colnames)
gaps = read.delim(args[2],
                  header = FALSE,
                  col.names = c(bed_colnames,
                                "name"))

blocks$len = blocks$end - blocks$start
gaps$len = gaps$end - gaps$start

data = data.frame("len" = blocks$len,
                  "type" = "blocks")
data = rbind(data, data.frame("len" = gaps$len,
                              "type" = "gaps"))

stat_blocks = as.data.frame(quantile(blocks$len))
stat_blocks$level = row.names(stat_blocks)
row.names(stat_blocks) = seq_len(nrow(stat_blocks))
stat_blocks$type = "blocks"
stat_gaps = as.data.frame(quantile(gaps$len))
stat_gaps$level = row.names(stat_gaps)
row.names(stat_gaps) = seq_len(nrow(stat_gaps))
stat_gaps$type = "gaps"
head(stat_blocks)
head(stat_gaps)
colnames(stat_blocks) = c("value", "level", "type")
colnames(stat_gaps) = c("value", "level", "type")
stat_df = rbind(stat_blocks, stat_gaps)
write.csv2(stat_df[, c("type", "level", "value")], file = args[3])

plot = ggplot(data = data[data$len < 5000, ], aes(x = len)) +
  geom_histogram() +
  facet_wrap(~ type) +
  theme_poster +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Size") +
  ylab("Count")
ggsave(args[4], height = 10, width = 10)