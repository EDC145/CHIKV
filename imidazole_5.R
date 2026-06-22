library(bio3d)
library(tidyr)
library(dplyr)
pdb <- read.pdb("grepCleaned3n42_noe3.pdb")
pdb <- data.frame(resNo = pdb$atom$resno, resid = pdb$atom$resid, atom = pdb$atom$elety, atomid = pdb$atom$eleno, chain = pdb$atom$chain, x = pdb$atom$x, y = pdb$atom$y, z = pdb$atom$z)

pdb <- pdb[pdb$resid %in% c("ARG","HIS","LYS"),]

hispdb <- pdb[pdb$resid %in% c("HIS"),] 
hispdb <- hispdb[hispdb$atom %in% c("ND1","NE2","CE1","CG","CD2"),]
hispdb$rc <- paste(hispdb$resNo,hispdb$chain)

hiscentre <- data.frame(
  ResNo = integer(length(unique(hispdb$rc))),
  ResId = character(length(unique(hispdb$rc))),
  chain = character(length(unique(hispdb$rc))),
  x = numeric(length(unique(hispdb$rc))),
  y = numeric(length(unique(hispdb$rc))),
  z = numeric(length(unique(hispdb$rc))),
  rc = unique(hispdb$rc),
  stringsAsFactors = FALSE
)

# Compute centroid for each HIS residue
for (rc_value in unique(hispdb$rc)) {
  subset_hispdb <- hispdb[hispdb$rc == rc_value, ]
  hiscentre[hiscentre$rc == rc_value, "ResNo"] <- subset_hispdb$resNo[1]
  hiscentre[hiscentre$rc == rc_value, "ResId"] <- as.character(subset_hispdb$resid[1])
  hiscentre[hiscentre$rc == rc_value, "chain"] <- as.character(subset_hispdb$chain[1])
  hiscentre[hiscentre$rc == rc_value, "x"] <- mean(subset_hispdb$x)
  hiscentre[hiscentre$rc == rc_value, "y"] <- mean(subset_hispdb$y)
  hiscentre[hiscentre$rc == rc_value, "z"] <- mean(subset_hispdb$z)
}

cutoff <- 10   ###ChangeCutoffHere

wow <- data.frame(
  ResNo1 = integer(),
  ResId1 = character(),
  chain1 = character(),
  dis = numeric(),
  ResNo2 = integer(),
  ResId2 = character(),
  chain2 = character(),
  ATOM = character(),
  stringsAsFactors = FALSE
)


for (i in 1:nrow(hiscentre)){
  for (j in 1:nrow(pdb)) {
    
    d <- sqrt((hiscentre$x[i] - pdb$x[j])^2 + (hiscentre$y[i] - pdb$y[j])^2 + (hiscentre$z[i] - pdb$z[j])^2)
    #if (d <= cutoff & !(lyspdb$resNo[i] == pdb$resNo[j] & lyspdb$chain[i] == pdb$chain[j]) & !(lyspdb$atom[i] %in% c("C","CA","N","O")) & (pdb$atom[j] %in% c("C","CA","N","O"))){
    
    if (d <= cutoff & !(hiscentre$ResNo[i] == pdb$resNo[j] & hiscentre$chain[i] == pdb$chain[j]) & (pdb$atom[j] %in% c("NE", "NH1","NH2","HE","1HH1","2HH1","1HH2","2HH2","NZ","1HZ","2HZ","3HZ","ND1","NE2","HE2"))){ #(hiscentre$atom[i] %in% c("ND1","NE2","HE2","CE1","CG","CD2","HE1","HD2")) & 
      
      
      
      wow <- wow %>%
        add_row(ResNo1 = hiscentre$ResNo[i], ResId1 = hiscentre$ResId[i], chain1 = hiscentre$chain[i], dis = d, ResNo2 = pdb$resNo[j], ResId2 = pdb$resid[j], chain2 = pdb$chain[j], ATOM = pdb$atom[j])
    }
    
  }
}

wow <- arrange(wow, dis)


wow2 <- wow %>%
  rowwise() %>%
  mutate(
    ResNoA = min(ResNo1, ResNo2),
    ResIdA = ifelse(ResNo1 < ResNo2, ResId1, ResId2),
    ChainA = ifelse(ResNo1 < ResNo2, chain1, chain2),
    ResNoB = max(ResNo1, ResNo2),
    ResIdB = ifelse(ResNo1 > ResNo2, ResId1, ResId2),
    ChainB = ifelse(ResNo1 > ResNo2, chain1, chain2),
  ) %>%
  distinct(ResNoA, ResIdA, ChainA, ResNoB, ResIdB, ChainB, .keep_all = TRUE) %>%
  select(ResNo1, ResId1, chain1, dis, ResNo2, ResId2, chain2, ATOM)  # Keep original column names




c2c <- read.pdb("centroidfor_CHIKV.pdb")
c2c <- data.frame(resNo = c2c$atom$resno, resid = c2c$atom$resid, atom = c2c$atom$elety, atomid = c2c$atom$eleno, chain = c2c$atom$chain, x = c2c$atom$x, y = c2c$atom$y, z = c2c$atom$z)

c2c <- c2c[c2c$resid %in% c("ARG","HIS","LYS"),]

c2cd <- data.frame(
  ResNo1 = integer(),
  ResId1 = character(),
  chain1 = character(),
  disCentroid = numeric(),
  ResNo2 = integer(),
  ResId2 = character(),
  chain2 = character(),
  ATOM = character(),
  stringsAsFactors = FALSE
)

nrow(c2c)
c2c$x[112]

pairs <- expand.grid(i = 1:nrow(c2c), j = 1:nrow(c2c)) %>%
  filter(i != j) 

# Compute distances and filter by cutoff
c2cd <- pairs %>%
  rowwise() %>%
  mutate(d = sqrt((c2c$x[i] - c2c$x[j])^2 + (c2c$y[i] - c2c$y[j])^2 + (c2c$z[i] - c2c$z[j])^2)) %>%
  mutate(ResNo1 = c2c$resNo[i], ResId1 = c2c$resid[i], chain1 = c2c$chain[i],
         ResNo2 = c2c$resNo[j], ResId2 = c2c$resid[j], chain2 = c2c$chain[j],
         ATOM = c2c$atomid[j]) %>%
  select(ResNo1, ResId1, chain1, disCentroid = d, ResNo2, ResId2, chain2, ATOM)



c2cd <- arrange(c2cd, disCentroid)



wow3 <- wow2 %>%
  left_join(c2cd %>% select(ResNo1, ResId1, chain1, ResNo2, ResId2, chain2, disCentroid), 
            by = c("ResNo1", "ResId1", "chain1", "ResNo2", "ResId2", "chain2"))





int <- read.table("Location_of_IntFc")
int$com <- paste(int$V1, int$V2)
wow3$IsIntfc <- ""
for (i in 1:nrow(wow3)){
  if( (paste(wow3$chain1[i], wow3$ResNo1[i]) %in% int$com) & (paste(wow3$chain2[i], wow3$ResNo2[i]) %in% int$com)){
    wow3$IsIntfc[i] <- "yes,yes"
  }
  else if( (paste(wow3$chain1[i], wow3$ResNo1[i]) %in% int$com) & !(paste(wow3$chain2[i], wow3$ResNo2[i]) %in% int$com)){
    wow3$IsIntfc[i] <- "yes,no"
  }
  else if( !(paste(wow3$chain1[i], wow3$ResNo1[i]) %in% int$com) & (paste(wow3$chain2[i], wow3$ResNo2[i]) %in% int$com)){
    wow3$IsIntfc[i] <- "no,yes"
  }
  else{
    wow3$IsIntfc[i] <- "no"
  }
}

#propka <-  read.table("propka_DENV_summ")
wow3$dis <- round(wow3$dis, 4)
wow3$disCentroid <- round(wow3$disCentroid, 4)



pka <- read.table("pka_Result_summ", header = TRUE)
wow4 <- wow3 %>%
  mutate(com = paste(ResNo1, chain1, sep = "")) %>%
  left_join(pka %>% rename(com = RESIDUE) %>% select(com, pKa), by = "com")



asa <- read.table("full.asa", col.names = c("s","ResNum", "Restype", "chain", "AllatomsSum", "AllatomsPer.", "Non_P_side_Sum", "Non_P_side_Per.", "Polar_Side_Sum", "Polar_Side_Per.", "Total_Side_Sum", "Total_Side_Per.", "Main_Chain_Sum", "Main_Chain_Per."))
asa <- asa %>% select(ResNum,Restype,chain,AllatomsPer.) %>% rename(ResNo1 = ResNum, Resid1 = Restype, chain1 = chain, ASA_allAtomsPer = AllatomsPer.) %>% mutate(com = paste(ResNo1, chain1, sep = ""))
wow4 <- left_join(wow4, (asa %>% select(com, ASA_allAtomsPer)), by = "com")


depth <- read.table("DepthCHIKV-residue.depth", 
                    col.names = c("ch:no", "resid", "all-atom", "all-atom(stdev)", "MC-atom", 
                                  "MC-atom(stdev)", "SC-atom", "SC-atom(stdev)", "SC-polar-atom", 
                                  "SC-polar-atom(stdev)", "SC-nonpolar-atom", "SC-nonpolar-atom(stdev)"))

depth <- depth %>%
  separate(ch.no, c("chain", "resNo"), ":") %>%
  mutate(com = paste(resNo, chain, sep = "")) %>%
  rename(DepthValue = all.atom) %>%
  select(com, DepthValue)

wow4 <- wow4 %>%
  left_join(depth, by = "com")

c2c2 <- read.pdb("centroidfor_CHIKV.pdb")
c2c2 <- data.frame(
  resNo = c2c2$atom$resno,
  resid = c2c2$atom$resid,
  atom = c2c2$atom$elety,
  atomid = c2c2$atom$eleno,
  chain = c2c2$atom$chain,
  x = c2c2$atom$x,
  y = c2c2$atom$y,
  z = c2c2$atom$z
)

# Generate all pairs, excluding self-pairs
pairs <- expand.grid(i = 1:nrow(c2c2), j = 1:nrow(c2c2)) %>%
  filter(i < j)  # i < j avoids duplicate comparisons and self-pairs


cutoff_n = 7  ## put neighbor cutoff here

# Compute distances and filter by cutoff
c2cd2 <- pairs %>%
  rowwise() %>%
  mutate(
    d = sqrt((c2c2$x[i] - c2c2$x[j])^2 + 
               (c2c2$y[i] - c2c2$y[j])^2 + 
               (c2c2$z[i] - c2c2$z[j])^2),
    ResNo1 = c2c2$resNo[i], 
    ResId1 = c2c2$resid[i], 
    chain1 = c2c2$chain[i],
    ResNo2 = c2c2$resNo[j], 
    ResId2 = c2c2$resid[j], 
    chain2 = c2c2$chain[j],
    ATOM = c2c2$atomid[j]
  ) %>%
  ungroup() %>%
  filter(d <= cutoff_n) %>%
  select(ResNo1, ResId1, chain1, disCentroid = d, ResNo2, ResId2, chain2, ATOM)

# Add combined ResNo and chain column
c2cd2 <- c2cd2 %>%
  mutate(com = paste(ResNo1, chain1, sep = ""), com2 = paste(ResNo2, chain2, sep = ""))

wow5 <- wow4 %>%
  rowwise() %>%
  mutate(!!paste0("neighbors", cutoff_n) := paste(
    unique(c(
      paste(c2cd2$ResId2[c2cd2$com == com], c2cd2$com2[c2cd2$com == com]),
      paste(c2cd2$ResId1[c2cd2$com2 == com], c2cd2$com[c2cd2$com2 == com])
    )),
    collapse = ", "
  ))


writexl::write_xlsx(wow5, "distances_withPka_andNeighborsAndSASA.xlsx")


#writexl::write_xlsx(wow2, paste0("FromCentroidofImidazole",cutoff,".xlsx"))
#writexl::write_xlsx(wow, paste0("FromCentroidofImidazole_AlL",cutoff,".xlsx"))


plot(wow5$ASA_allAtomsPer, wow5$pKa)
plot(wow5$DepthValue, wow5$pKa)
plot(wow5$DepthValue, wow5$ASA_allAtomsPer)
