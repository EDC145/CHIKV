library(bio3d)
a <-read.pdb("/home/lab12/testROSETTA/chikv_rosetta/Structure/centroidAll/centroidfor_CHIKV.pdb")


print("copy paste the underpacked residues from my chosen seeds")


underpk_centroid <- read.pdb("/home/lab12/testROSETTA/chikv_rosetta/undpk/Voronoia/voronoia_3n42NoH/neighbors/centroid_chikv_undpk_noE3.pdb")

#seedUndB <- scan()
#seedUndF <- scan()

seedUndB <- c(7,8,35,50,60,70,99,111,123,131,139,141,153,156,160,162,167,169,171,172,190,201,225,227,231,234,248,255,256,257,258,259,260,262,264,266,268,269,270,271,285,287,295,317,330)
seedUndB<-paste(seedUndB,"B")
seedUndF <- c(18,28,29,38,40,46,48,49,51,60,100,104,106,131,133,137,161,163,171,203,204,205,221,261,287,329,331,368)
seedUndF <-paste(seedUndF,"F")

print("enter B chain seeds Intfc")
#seedsIntB <- scan()
seedsIntB <- c(18,29,38,134,154,165,176,226,238,240,243,300,306,42,36,298,168,261,152,138,170,177,301,40,39,308,136,151,334)
seedsIntB <- paste(seedsIntB, "B")

print("enter F chain seeds")
#seedsIntF <- scan()
seedsIntF <- c(50,57,105,115,117,181,241,253,256,88,54,112,254,116,231,228,230,260,92,93,259,249,113,95,52,89,244,245,247,183,252,223,224,225,226,227,229)
seedsIntF <- paste(seedsIntF, "F")

pdb <- data.frame(res = paste(a$atom$resno, a$atom$chain), x = a$atom$x, y= a$atom$y, z = a$atom$z)
undpdb <- data.frame(res = paste(underpk_centroid$atom$resno, underpk_centroid$atom$chain), x = underpk_centroid$atom$x, y= underpk_centroid$atom$y, z = underpk_centroid$atom$z)
d <- 0
check <- vector()

for(i in c(seedUndB, seedUndF)){
  for (j in undpdb$res) {
    d <- sqrt((pdb$x[pdb$res == i] - pdb$x[pdb$res == j])^2 + (pdb$y[pdb$res == i] - pdb$y[pdb$res == j])^2 + (pdb$z[pdb$res == i] - pdb$z[pdb$res == j])^2)
    if (d <= 7) {
      check <- unique(c(check, j,i))
    }
   }
}

### Check above  > WORKING


dB <- data.frame(matrix(ncol = length(c(check)), nrow = length(c(seedsIntB, seedsIntF))))
rownames(dB) <- c(seedsIntB, seedsIntF)
colnames(dB) <- c(check)

#for (i in c(seedsB,seedsF)) {
  #for(j in c(epitopeB,epitopeF)){
  #  d <- sqrt((pdb$x[pdb$res == i]^2 - pdb$x[pdb$res == j]^2)+(pdb$y[pdb$res == i]^2 - pdb$y[pdb$res == j]^2) + (pdb$z[pdb$res == i]^2 - pdb$z[pdb$res == j]^2))
#    dB[i,j] <- d
#  }
  
#}



for (i in c(seedsIntB, seedsIntF)) {
  for (j in c(check)) {
    # Subset coordinates based on residue and chain
    coord_i <- pdb[pdb$res == i, c("x", "y", "z")]
    coord_j <- pdb[pdb$res == j, c("x", "y", "z")]
    d <- sqrt((coord_i$x - coord_j$x)^2 + (coord_i$y - coord_j$y)^2 + (coord_i$z - coord_j$z)^2)
    dB[i, j] <- d
  }
}


openxlsx::write.xlsx(dB, "distancefromUnderpackedAndIntFcSeeds.xlsx",rowNames = TRUE)
