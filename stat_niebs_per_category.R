args = commandArgs(trailingOnly = TRUE)

data = read.delim(args[1],
                  header = FALSE,
                  col.names = c("chr", "start", "end", "name"))

total = nrow(data)
categories = aggregate(data$name, by = list(data$name), FUN = length)
colnames(categories) = c("category", "count")
print(total)

categories$rel = categories$count / total
write.csv2(categories, file = args[2])