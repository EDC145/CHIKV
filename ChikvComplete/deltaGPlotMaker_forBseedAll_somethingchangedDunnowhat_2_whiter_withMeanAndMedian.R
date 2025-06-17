library(dplyr)
library(openxlsx)
library(readODS)
library(ggplot2)
library(gridExtra)
library(cowplot)



#path1 <- list.files( pattern = ".ods", ignore.case = TRUE)

#path2 <- sub("\\.ods.*", "",path1)
#path2 <- unique(path2)

#for (pathi in path2){
  

pathi <- "CHIKV COmplete" 
  
 # fileToRead <- paste0(pathi,".ods")
 # d <- read.xlsx("allSeedCHIKV_B.ods")
d <- read.xlsx("CHIKV COmplete.xlsx")



d2 <- d %>% distinct(paste(PDB_File,seed), .keep_all = TRUE)


binw <- 2

min_val <- round(min(d2$`difference(control-new)`, na.rm = TRUE),3) 
max_val <- round(max(d2$`difference(control-new)`, na.rm = TRUE),3)

start <- -11  # Round down to ensure coverage
end <- 21

hist_data <- ggplot_build(
  ggplot(d2, aes(x = `difference(control-new)`)) + 
    geom_histogram(binwidth = binw, boundary = start)
)$data[[1]]

max_count <- max(hist_data$count, na.rm = TRUE)  # Get the highest bin count

meanx <- mean(d2$`difference(control-new)`)
median_x <- median(d2$`difference(control-new)`)
sd_X <- sd(d2$`difference(control-new)`)
meanAbsoluteDeviation_forMedian <- mad(d2$`difference(control-new)`)

#d2 <- d2 %>% filter(!is.na(`difference(control-new)`))
#ggplot(d2, aes(x = `difference(control-new)`)) + geom_freqpoly(bins = 1)

#hist(d2, xlab = "scores", col = "pink", border = "black")

histogram <-   ggplot(d2, aes(x = `difference(control-new)`)) + 
  geom_vline(xintercept = 3, color = "darkturquoise", linetype = "dashed", size = 0.7, alpha = 0.8)+
  geom_vline(xintercept = meanx, color = "darkorchid", linetype = "dashed", size = 0.7, alpha = 0.8)+
  geom_vline(xintercept = meanx+sd_X, color = "plum", linetype = "dotted", size = 0.5)+
  geom_vline(xintercept = meanx-sd_X, color = "plum", linetype = "dotted", size = 0.5)+
  geom_vline(xintercept = median_x, color = "green4", linetype = "dashed", size = 0.7, alpha = 0.8)+
  geom_vline(xintercept = median_x-meanAbsoluteDeviation_forMedian  , color = "olivedrab3", linetype = "dotted", size = 0.5)+
  geom_vline(xintercept = median_x+meanAbsoluteDeviation_forMedian  , color = "olivedrab3", linetype = "dotted", size = 0.5)+
#  geom_hline(yintercept = 0, color = "gray32", size = 0.5) +
#  geom_vline(xintercept = -21, color = "gray32", size = 0.5)+
  geom_histogram(binwidth = binw, color = "gray10", fill = NA, boundary = start) +
  geom_text(stat = "bin", binwidth = binw,   aes(label = ifelse(after_stat(count) > 30, after_stat(count), "")), 
            position = position_nudge(y = 50), color = "gray10", boundary = start, size =3.0) + # 
  scale_x_continuous(limits = c(start, end), breaks = seq(start, end, by = 2)) +
    scale_y_continuous(limits = c((0), (2400)), expand = c(0,0), breaks = seq(0, 2400, by = 300))+ 
  labs(title = paste0("CHIKV Complete"),x = "Δ Rosetta Score (WT - Mut)", y = "Frequency") + 
  theme_minimal() +
  theme(axis.title.x = element_text( vjust = 0, size =14),
        axis.title.y = element_text( vjust = 0, size =14),
        panel.grid.major = element_blank(),     # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),  # <-- Increased x-axis tick text
        axis.text.y = element_text(size = 11),  # <-- Increased y-axis tick text
       axis.line = element_line(color = "gray32"),  # Add axis lines
         axis.ticks = element_line(color = "gray32"),  # Add axis ticks# Remove minor grid lines)  # Add axis ticks
     plot.background = element_rect(fill = "white", color = NA),   # plot area background
     panel.background = element_rect(fill = "white", color = NA)   # panel background
     )  +annotate("text", x = 12 , y = 2100, label = paste("Mean:", round(meanx,3)), hjust = 0, color = "darkorchid") +
  annotate("text", x = 12 , y = 1950, label = paste("SD: ±", round(sd_X,3)), hjust = 0, color = "plum") +
  annotate("text", x = 12 , y = 1800, label = paste("Median:", round(median_x,3)), hjust = 0, color = "green4") +
  annotate("text", x = 12 , y = 1650, label = paste("Mean.Abs.Dev.: ±", round(meanAbsoluteDeviation_forMedian,3)), hjust = 0, color = "olivedrab3") +
  annotate("text", x = 12 , y = 1500, label = paste("Cutoff:", "3.0"), hjust = 0, color = "darkturquoise")

 # annotate("text", x = max_val/(4/3) , y = 1200, label = paste("Max:", max_val), hjust = 0, color = "darkgreen")

  ggsave("CHIKVpng_filename3_white.png", plot =histogram,
         width = 10, height = 5, dpi = 600)
#ggsave2(paste0(pathi, " graph.pdf"), histogram)
#}

#  d2 %>%    filter(is.na(`difference(control-new)`) |              `difference(control-new)` < (start-1) |              `difference(control-new)` > (end+1))

