library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
source("~/setThemePoster.R")
print("loading datasets")
microsat = read.delim(args[1],
                      header = FALSE,
                      col.names = c("chr",
                                    "start",
                                    "end",
                                    "name",
                                    "score",
                                    "strand"),
                      comment.char = "#")
microsat$gc_count = str_count(microsat$name, "[GC]")
microsat$gc_content = microsat$gc_count / nchar(microsat$name)

plot = ggplot(data = microsat,
              aes(x = gc_content)) +
  geom_histogram() +
  geom_vline(xintercept = quantile(microsat$gc_content),
             color = "red") +
  theme_poster +
  xlab("GC content") +
  ylab("Count")
ggsave(args[2])