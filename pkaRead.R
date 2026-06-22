txt <- readLines("Result_summ")  
txt <- gsub(" +", "\t", txt)  # Replace multiple spaces with a single tab
writeLines(txt, "cleaned_Result_summ.txt")  
a <- read.delim("cleaned_Result_summ.txt", header = TRUE, sep = "\t")
