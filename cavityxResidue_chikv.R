## THIS FROM .vor

HOLENUMBER1 <-c(132, 137, 143 ,147 ,155 ,160 ,377 ,381 ,382, 383, 388 ,390 ,391 ,395, 399, 460, 462)
HOLENUMBER2 <-c(1746 ,1749, 1762, 1764, 1765 ,1767 ,1790, 2116 ,2122 ,2125, 2128 ,3794 ,3796 ,3799 ,3802 ,3803)
HOLENUMBER3 <-c(3347, 4097 ,4098, 4112, 4114, 4115 ,4201 ,4203, 4212 ,4216, 4345 ,4346 ,5211 ,5227)
HOLENUMBER4 <-c(3227, 3288 ,3305, 5590, 5591, 5592, 5600 ,5604, 5885 ,5886, 5887, 5889 ,5901)
HOLENUMBER5 <-c(3452 ,3456 ,3458 ,3469 ,3944 ,3945 ,4654 ,4655 ,4658 ,4666 ,4670)
HOLENUMBER6 <-c(5429, 5433, 5478, 5482 ,5483, 5794, 5798, 5803, 5808, 5814)
HOLENUMBER7<- c(3467 ,3471 ,3479, 3486, 4647 ,4648 ,4923 ,4928 ,4932)
HOLENUMBER8<- c(3488, 4613, 4645, 4647, 4762, 4765, 4767, 4909 ,4910)
HOLENUMBER9 <-c(255 ,257 ,1699 ,2302 ,2318 ,2319, 2320 ,2328)
HOLENUMBER10 <-c(383, 389 ,395 ,412 ,415, 416, 427, 462)
HOLENUMBER11 <-c(2272 ,3525 ,3527 ,3532 ,3536 ,3899 ,3904, 4871)
HOLENUMBER12 <-c(3446, 3450 ,3457, 3463, 4014 ,4022, 4026 ,4636)
HOLENUMBER13 <-c(3521 ,3912 ,4783 ,4798 ,4870 ,4883 ,4886)
HOLENUMBER14 <-c(3740, 3745, 3747, 3754, 3757 ,3803 ,3816 ,3820)
HOLENUMBER15 <-c(504 ,642 ,754 ,757 ,761 ,1140 ,2212)
HOLENUMBER16 <-c(1703, 2206, 2211, 2228, 2253, 2254, 3534)
HOLENUMBER17 <-c(1741 ,1748 ,3773 ,3776 ,3777 ,3781 ,3791)
HOLENUMBER18 <-c(2266, 2283, 2286, 3548 ,3550 ,3554, 3557)
HOLENUMBER19 <-c(2572 ,2574 ,2589, 2591, 2865, 2868, 2881)
HOLENUMBER20 <-c(3379, 3380 ,3397 ,3399 ,4076 ,4403 ,4418)
HOLENUMBER21 <-c(3897 ,3898 ,4796, 4803, 4870, 4875, 4878)
HOLENUMBER22 <-c(5433, 5447, 5450, 5462, 5478, 5803, 5808)
HOLENUMBER23 <-c(18, 24, 30, 234, 262, 265)
HOLENUMBER24 <-c(123, 131, 137, 142, 372, 391)
HOLENUMBER25 <-c(633 ,636 ,649 ,662 ,747, 752)
HOLENUMBER26 <-c(1059, 1061, 1227, 1230, 1237, 1242)
HOLENUMBER27 <-c(1657, 1662, 1665, 2374, 3959, 3962)
HOLENUMBER28 <-c(1716, 1722, 1727, 1734, 2149, 2281)
HOLENUMBER29 <-c(1789, 1856, 1921, 1943, 2113, 2118)
HOLENUMBER30 <-c(3232, 3244, 3248, 3281, 5290, 5292)
HOLENUMBER31 <-c(3462, 3463, 4627, 4637, 4928, 4929)
HOLENUMBER32 <-c(440, 450, 473, 474, 2356)
HOLENUMBER33 <-c(490, 496, 613, 772, 925)
HOLENUMBER34 <-c(1700, 1701, 1702, 2268, 2303)
HOLENUMBER35 <-c(1701 ,2268, 2274, 2278, 2303)
HOLENUMBER36 <-c(1703, 1714, 1726, 2203, 2206)
HOLENUMBER37 <-c(1709, 1716, 2152, 2163, 2281)
HOLENUMBER38 <-c(1956, 1961, 1969, 1976, 2052)
HOLENUMBER39 <-c(2504, 2590, 2607, 2673, 2960)
HOLENUMBER40 <-c(3312, 4132, 4148, 4331, 5240)
HOLENUMBER41 <-c(3380, 4401, 4403, 4410, 4418)
HOLENUMBER42 <-c(4493, 5105, 5108, 5111, 5120)
HOLENUMBER43 <-c(1039 ,1319 ,1323, 1335)
HOLENUMBER44 <-c(1566, 1569, 2417, 4489)
HOLENUMBER45 <-c(1703, 1726, 2206, 2254)
HOLENUMBER46 <-c(2266, 2286, 3557, 3815)
HOLENUMBER47 <-c(3462, 4637, 4648, 4928)
HOLENUMBER48 <-c(4115, 4203, 4327, 4346)
print(HOLENUMBER30)


library(bio3d)
pdb <- read.pdb("Relaxed_3n42_noH.vor.pdb")
a <- data.frame(
  ano = pdb$atom$eleno, 
  resNo = paste(pdb$atom$resno, pdb$atom$chain, sep = " "), 
  residue = pdb$atom$resid
)


hole <- list()

# Loop through the numbers 1 to 50 and dynamically access variables 1 to 48
for (i in 1:48) {       #### EDIT HERE
  var_name <- paste0("HOLENUMBER", i)  # Create variable name like "a1", "a2", etc.
  hole[[i]] <- get(var_name)  # Use get() to access the variable and store it in the list
}

which(a$ano %in% hole[[12]])

df <- data.frame(cav = 1:48)


for (i in 1:48) {
  # Get the positions where values in a$ano match those in hole[[i]]
  matching_positions <- which(a$ano %in% hole[[i]])
  
  # Convert the positions to a space-separated string and store in the 'values' column
  df$values[i] <- paste(a$resNo[matching_positions], collapse = ",")
}

for (i in 1:nrow(df)) {
  # Split the space-separated string into individual elements
  values_list <- unlist(strsplit(df$values[i], ","))
  
  # Remove duplicates
  unique_values <- unique(values_list)
  
  # Recombine the unique values into a space-separated string
  df$values[i] <- paste(unique_values, collapse = ",")
}

openxlsx::write.xlsx(df, "cavitiesAndResidues2_chikv.xlsx")
