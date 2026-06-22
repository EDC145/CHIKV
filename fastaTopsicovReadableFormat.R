a <- readLines("MSA.fasta")

substr(a[1], start = 1, stop = 1)


b<- a

c <-""
for (i in 1:length(a)){
  if(substr(a[i], start = 1, stop = 1)==">"){
    c <- paste(c,"\n", sep = "")
  }
  else{
    c <- paste(c, a[i], sep = "")
  }
  
}

write(c, "testAln")
