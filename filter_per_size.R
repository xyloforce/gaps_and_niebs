args = commandArgs(trailingOnly = TRUE)

data = read.delim(args[1], header = F)
data$len = data$V3 - data$V2

write.table(data,
            file = args[2],
            quote = FALSE,
            sep = "\t",
            row.names = F,
            col.names = F)
