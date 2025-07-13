args = commandArgs(trailingOnly=TRUE)

bed_colnames = c("chr", "start", "end", "name")

gaps = read.delim(args[1],
                  header = FALSE,
                  col.names = bed_colnames)

blocs = read.delim(args[2],
                   header = FALSE,
                   col.names = bed_colnames)

overlap = read.delim(args[3],
                     header = FALSE,
                     col.names = bed_colnames)
gaps$name = "gaps"
blocs$name = "blocks"
overlap$name = "overlap"

gaps = rbind(gaps, blocs)
gaps = rbind(gaps, overlap)

write.table(gaps,
            file = args[4],
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
