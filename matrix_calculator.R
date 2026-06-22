library(bio3d)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

a <- read.pdb("Relaxed_3n42_noH.pdb")
val <- read.table("Location_of_IntFc")

val <- paste(val$V1,":",val$V2,sep = "")

v <- val

v <- data.frame(v)
v[c("chain","resno")] <- str_split_fixed(v$v,":",2)
v<- v[c("chain","resno")]
v$resChain <- paste(v$resno, v$chain)

a1 <- a$atom
a1$resChain <- paste(a1$resno,a1$chain)
a4 <- a1

a1<- a1[!a1$elety %in% c("CA", "C", "N", "O"),]     #ADD THIS TO HAVE ONLY SIDE CHAINS


b <- data.frame(type = character(),
                eleno = numeric(),
                elety = character(),
                alt = character(),
                insert = character(),
               
                segid = character(),
                elesy = character(),
                charge = character(),
                resChain = character(),
                stringsAsFactors = FALSE)


x <- aggregate(x~resChain, a1, FUN=mean)
y <-  aggregate(y~resChain, a1, FUN=mean)
z <- aggregate(z~resChain, a1, FUN=mean)
bb <- aggregate(b~resChain, a1, FUN = mean)
b <- bind_rows(b, x)
b <- merge(b, y, by = "resChain")
b <- merge(b, z, by = "resChain")
b <- merge(b, bb, by = "resChain")
colnames(b)[which(names(b) == "bb")] <- "b"
a2 <- data_frame(chain= a1$chain, resno=a1$resno,resid= a1$resid, resChain=a1$resChain, o=a1$o)
a2 <- unique(a2[,c("chain","resno","resid",'resChain','o')])
b <- merge(b, a2, by = "resChain")
b$type <- "ATOM"
b$elety <- "CA"
b$alt <- ""
b$insert = ""
b$segid = ""
b$elesy <- ""
b$charge <- ""


a3 <- a4[a4$resid=="GLY",]
a3 <- a3[a3$elety=="CA",]
a3$alt<-""
a3$insert<- ""
a3$segid<- ""
a3$elesy<- ""
a3$charge<- ""

b2<- rbind(b, a3)
b <- b2

b<- b[with(b,order(chain, resno)),]
b$eleno <- 1:length(b$type)

#xyz <- rbind(b$x, b$y, b$z)

#write.pdb(file = "test.pdb", xyz = xyz, type = b$type, resno = b$resno,
 #         resid = b$resid, eleno = b$eleno, elety = b$elety, chain = b$chain, 
  #       o = b$o, b = b$b,  
   #       append = FALSE, verbose = FALSE, chainter = TRUE, end = TRUE, sse = FALSE, 
    #     print.segid = FALSE)


#for (i in 1:length(b$eleno)){
 #c <- b[b$eleno==i,]



###i have centroids for all atoms and residues
### now filtering for v which is our IntFc selected residues

bIntFc <- left_join(v, b, by="resChain")

bIntFc <- bIntFc[, !names(bIntFc) %in% "resno.x"]
bIntFc <- bIntFc[, !names(bIntFc) %in% "chain.x"]


colnames(bIntFc)[which(names(bIntFc) == "resno.y")] <- "resno"

colnames(bIntFc)[which(names(bIntFc) == "chain.y")] <- "chain"




bIntFc<- bIntFc[with(bIntFc,order(chain, resno)),]


#b <- b[, !names(b) %in% "resChain"]



 
#}
#class(b$eleno)








# Define the function to write the PDB file without using sprintf
write_pdb <- function(data, file) {
  pdb_lines <- apply(data, 1, function(row) {
    # Create a formatted string manually
    line <- paste0(
      format(row['type'], width = 6, justify = "left"),
      format(as.integer(row['eleno']), width = 5, justify = "right"),
      " ",
      format(row['elety'], width = 4, justify = "centre"),
      format(row['alt'], width = 1),
      format(row['resid'], width = 3, justify = "left"),
      " ",
      format(row['chain'], width = 1),
      format(as.integer(row['resno']), width = 4, justify = "right"),
      format(row['insert'], width = 1),
      "   ",
      formatC(as.numeric(row['x']), format = "f", digits = 3, width = 8),
      formatC(as.numeric(row['y']), format = "f", digits = 3, width = 8),
      formatC(as.numeric(row['z']), format = "f", digits = 3, width = 8),
      formatC(as.numeric(row['o']), format = "f", digits = 2, width = 6),
      formatC(as.numeric(row['b']), format = "f", digits = 2, width = 6),
      "           ",
      format(row['segid'], width = 2),
      format(row['charge'], width = 2)
    )
    return(line)
  })
  
  # Write lines to file
  writeLines(pdb_lines, file)
}
write_pdb(bIntFc, "centroidforIntfc_CHIKV.pdb")



p <- read.pdb("centroidforIntfc_CHIKV.pdb")


df <- data.frame(resid = as.character(p$atom$resid), resno = p$atom$resno, chain = as.character(p$atom$chain), x= p$atom$x, y = p$atom$y, z = p$atom$z)

df$chain <- as.character(df$chain)
df$resid <- as.character(df$resid)
df$reschain <- paste(df$resno,df$chain)


m <- matrix(0,nrow = length(df$resno), ncol =  length(df$resno))


df$x[1]-df$x[1]

cl <- data.frame(reschain = df$reschain)
cl$neig <-""


for(i in 1:length(df$resno)){
  for (j in 1:length(df$resno)){
    
    c = sqrt((df$x[i]-df$x[j])^2 + (df$y[i]-df$y[j])^2 + (df$z[i]-df$z[j])^2)
    m[i,j] <- c
    
    if((c<=6) && (df$reschain[i]!=df$reschain[j])){
      cl$neig[i] <- paste(cl$neig[i],df$reschain[j]," ,")
    }
    
    
  }
}

#rownames(m) <- paste(df$resno, df$chain)
#colnames(m) <- paste(df$resno, df$chain)


#write.table(m, "distanceMatrix.xlsx", row.names = TRUE, col.names = TRUE)


write.xlsx(cl,"6AngstormNeighborsIntfc_CHIKV_newResAdded.xlsx")



