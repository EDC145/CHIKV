library(bio3d)
a <-read.pdb("/home/lab12/testROSETTA/chikv_rosetta/Structure/centroidAll/centroidfor_CHIKV.pdb")


print("copy paste the epitope residues")
#epitope <- scan()


epitopeB <- c(58 , 59  ,71 , 72  ,74 , 75  ,76 , 77, 118, 119 ,183, 184, 185, 186, 187, 189, 190, 191, 193, 196, 198, 199, 203, 204, 206, 207, 208, 209, 210, 211, 213, 214, 215 ,216, 217 ,218 ,219, 220, 221, 230, 231 ,232, 233 ,245, 246, 247, 248 ,250 ,251)
epitopeB<-paste(epitopeB,"B")
epitopeF <- c(61,63,64,65,67)
epitopeF <-paste(epitopeF,"F")

print("enter B chain seeds")
seedsB <- scan()
seedsB <- paste(seedsB, "B")

print("enter F chain seeds")
seedsF <- scan()
seedsF <- paste(seedsF, "F")

pdb <- data.frame(res = paste(a$atom$resno, a$atom$chain), x = a$atom$x, y= a$atom$y, z = a$atom$z)

dB <- data.frame(matrix(ncol = length(c(epitopeB,epitopeF)), nrow = length(c(seedsB, seedsF))))
rownames(dB) <- c(seedsB, seedsF)
colnames(dB) <- c(epitopeB,epitopeF)

#for (i in c(seedsB,seedsF)) {
  #for(j in c(epitopeB,epitopeF)){
  #  d <- sqrt((pdb$x[pdb$res == i]^2 - pdb$x[pdb$res == j]^2)+(pdb$y[pdb$res == i]^2 - pdb$y[pdb$res == j]^2) + (pdb$z[pdb$res == i]^2 - pdb$z[pdb$res == j]^2))
#    dB[i,j] <- d
#  }
  
#}



for (i in c(seedsB, seedsF)) {
  for (j in c(epitopeB, epitopeF)) {
    # Subset coordinates based on residue and chain
    coord_i <- pdb[pdb$res == i, c("x", "y", "z")]
    coord_j <- pdb[pdb$res == j, c("x", "y", "z")]
    d <- sqrt((coord_i$x - coord_j$x)^2 + (coord_i$y - coord_j$y)^2 + (coord_i$z - coord_j$z)^2)
    dB[i, j] <- d
  }
}


openxlsx::write.xlsx(dB, "distancefromEpitopeAndIntFcSeeds.xlsx",rowNames = TRUE)
