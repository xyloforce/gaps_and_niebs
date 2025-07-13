library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

unique_part_of_overlap = read.delim(args[1],
                                    header = FALSE,
                                    col.names = c("chr",
                                                  "start",
                                                  "end",
                                                  "name"))

nieb_microsat_uniq = read.delim(args[2],
                                header = FALSE,
                                col.names = c("chr",
                                              "start",
                                              "end",
                                              "name"))

microsat_gap_niebs = read.delim(args[3],
                                header = FALSE,
                                col.names = c("chr",
                                              "start",
                                              "end",
                                              "name",
                                              "score",
                                              "strand"))

niebs_annotated = read.delim(args[4],
                             header = FALSE,
                             col.names = c("chr",
                                           "start",
                                           "end",
                                           "name"))

microsat_gap_niebs = microsat_gap_niebs[grep("overlap",
                                             microsat_gap_niebs$name), ]
niebs_annotated = niebs_annotated[niebs_annotated$name == "overlap", ]
microsat_gap_niebs$len = microsat_gap_niebs$end - microsat_gap_niebs$start
niebs_annotated$len = niebs_annotated$end - niebs_annotated$start
cov_microsat = sum(microsat_gap_niebs$len) / sum(niebs_annotated$len)

unique_part_of_overlap$len = unique_part_of_overlap$end - unique_part_of_overlap$start
nieb_microsat_uniq$len = nieb_microsat_uniq$end - nieb_microsat_uniq$start
cov_overlap = sum(nieb_microsat_uniq$len) / sum(unique_part_of_overlap$len)

results =  data.frame(names = c("cov_overlap",
                                "cov_microsat",
                                "cov_ratio"),
                      value = c(cov_overlap,
                                cov_microsat,
                                cov_overlap / cov_microsat))
write.csv2(results,
           file = args[5])