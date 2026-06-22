library(openxlsx)
library(dplyr)

### modify asa file reading

oan <- read.table("full.asa", col.names = c("s","ResNum", "Restype", "chain", "AllatomsSum", "AllatomsPer.", "Non_P_side_Sum", "Non_P_side_Per.", "Polar_Side_Sum", "Polar_Side_Per.", "Total_Side_Sum", "Total_Side_Per.", "Main_Chain_Sum", "Main_Chain_Per."))
oan_A <- filter(oan, chain %in% c('A','B')) ###from diner
oan_B <- filter(oan, (chain %in% c('F'))) ### From dimmer



oanA <- read.table("e3e2.asa", col.names = c("s","ResNum", "Restype", "chain", "AllatomsSum", "AllatomsPer.", "Non_P_side_Sum", "Non_P_side_Per.", "Polar_Side_Sum", "Polar_Side_Per.", "Total_Side_Sum", "Total_Side_Per.", "Main_Chain_Sum", "Main_Chain_Per."))
### from monomer


oanB <- read.table("e1.asa", col.names = c("s","ResNum", "Restype", "chain", "AllatomsSum", "AllatomsPer.", "Non_P_side_Sum", "Non_P_side_Per.", "Polar_Side_Sum", "Polar_Side_Per.", "Total_Side_Sum", "Total_Side_Per.", "Main_Chain_Sum", "Main_Chain_Per."))
## from monomer



diffA <- oanA$AllatomsPer.-oan_A$AllatomsPer.
diffB <- oanB$AllatomsPer.-oan_B$AllatomsPer.



df <- data.frame(resNo = oanA$ResNum, resid = oanA$Restype,chain=oanA$chain, difference_SA= diffA, monomerASA= oanA$AllatomsPer., dimerASA = oan_A$AllatomsPer.)
df2 <- data.frame(resNo = oanB$ResNum, resid = oanB$Restype,chain=as.character(oanB$chain), difference_SA= diffB, monomerASA= oanB$AllatomsPer., dimerASA = oan_B$AllatomsPer.)


df3 <- rbind(df, df2)


for(i in 1:nrow(df3)){
  df3$perchangeMonoMinusDimerDivideMonox100[i] <- ((df3$monomerASA[i]-df3$dimerASA[i])/df3$monomerASA[i])*100 
}

levels(df3$chain) <- c(levels(df3$chain), "F")
df3$chain[df3$chain == "FALSE"] <- "F"

write.xlsx(df3, "ASADifference2.xlsx")






### 
