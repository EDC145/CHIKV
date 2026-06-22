# Create a sequence of numbers from 1 to 394
numbersB <- 5:342
numbersF <- 1:393
# Create the combined labels for each number

result <- c(paste("B", numbersB), paste("F", numbersF))
# Print the result
writeLines(result, "Location.txt")

