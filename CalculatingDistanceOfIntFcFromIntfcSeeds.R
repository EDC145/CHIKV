library(bio3d)
a <-read.pdb("/home/lab12/testROSETTA/chikv_rosetta/Structure/centroidAll/centroidfor_CHIKV.pdb")


print("copy paste the intfc residues")
#epitope <- scan()

intfc<- read.table("Location_of_IntFc")
intfc<- paste(intfc$V2, intfc$V1)

print("enter B chain seeds")
seedsB <- scan()
seedsB <- paste(seedsB, "B")

print("enter F chain seeds")
seedsF <- scan()
seedsF <- paste(seedsF, "F")

pdb <- data.frame(res = paste(a$atom$resno, a$atom$chain), x = a$atom$x, y= a$atom$y, z = a$atom$z)

dB <- data.frame(matrix(ncol = length(intfc), nrow = length(c(seedsB, seedsF))))
rownames(dB) <- c(seedsB, seedsF)
colnames(dB) <- intfc

#for (i in c(seedsB,seedsF)) {
  #for(j in c(epitopeB,epitopeF)){
  #  d <- sqrt((pdb$x[pdb$res == i]^2 - pdb$x[pdb$res == j]^2)+(pdb$y[pdb$res == i]^2 - pdb$y[pdb$res == j]^2) + (pdb$z[pdb$res == i]^2 - pdb$z[pdb$res == j]^2))
#    dB[i,j] <- d
#  }
  
#}



for (i in c(seedsB, seedsF)) {
  for (j in intfc) {
    # Subset coordinates based on residue and chain
    coord_i <- pdb[pdb$res == i, c("x", "y", "z")]
    coord_j <- pdb[pdb$res == j, c("x", "y", "z")]
    d <- sqrt((coord_i$x - coord_j$x)^2 + (coord_i$y - coord_j$y)^2 + (coord_i$z - coord_j$z)^2)
    dB[i, j] <- d
  }
}


openxlsx::write.xlsx(dB, "distancefromInterfaceAndIntFcSeeds.xlsx",rowNames = TRUE)
