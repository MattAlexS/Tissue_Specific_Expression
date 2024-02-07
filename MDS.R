library(MASS)
library(magrittr)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(rgl)

full <- read.csv(file = file.choose(),header = TRUE)
tissues <- read.csv(file = file.choose(),header = TRUE)
rem_dup <- unique(full)
lab <- tissues$Adrenal.Gland
rem_dup <- rem_dup[ -c(49,51,57) ]
rem_dup$V1 <- rem_dup$V1 + 1

data <- full[,2:60]
mds <- data %>%
  dist() %>%          
  isoMDS() %>%
  .$points %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2",
          label = factor(lab),
          title = "Enriched Tissue Centers",
          size = 1,
          repel = TRUE)


for (i in 1:32) {
  tester <- filter(rem_dup, V1 == i)
  t = paste("EnrichedTis", i, sep = "")
  tf = paste(t, ".png", sep = "")
  if(nrow(tester) > 5) {
    data = tester[,(2:60)]
    mds <- data %>%
      dist() %>%          
      isoMDS() %>%
      .$points %>%
      as_tibble()
    colnames(mds) <- c("Dim.1", "Dim.2")
    # Plot MDS
    dp <- ggscatter(mds, x = "Dim.1", y = "Dim.2",
              title = t,
              size = 1,
              repel = TRUE)
    ggarrange(dp) %>%
      ggexport(filename = tf,
      width = 480,
      height = 480,
      pointsize = 12,
      res = NA,
      verbose = TRUE
    )
  }
}

mds <- data %>%
  dist() %>%          
  cmdscale(k=3) %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2", "Dim.3")


plot3d(mds$Dim.1, mds$Dim.2, mds$Dim.3, xlab = "Dim 1", ylab = "Dim 2", zlab = "Dim 3")
text3d(mds$Dim.1, mds$Dim.2, mds$Dim.3,lab)
title3d(main = "Enriched Centers")

#ggscatter(mds, x = "Dim.1", y = "Dim.2", z= "Dim.3", 
#          size = 1,
#          repel = TRUE)

#Kruskal non-metric
mds <- data %>%
  dist() %>%          
  isoMDS() %>%
  .$points %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2",
          col = factor(lab),
          size = 1,
          repel = TRUE)

#sammon non-metric
mds <- data %>%
  dist() %>%          
  sammon() %>%
  .$points %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2",
          col = factor(lab),
          size = 1,
          repel = TRUE)