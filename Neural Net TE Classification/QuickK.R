#k-means clustering of tissue type
library(factoextra)
library(mclust)

data <- read.csv(file = file.choose(), header = FALSE)


cdb <- data[,2:ncol(data)]
tissues <- data[,1]

rownames(cdb) <- tissues


km.res <- kmeans(cdb, 32, iter.max = 10)

#adjustedRandIndex(tissues,km.res$cluster)

res.pca = prcomp(cdb, scale = FALSE)
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

print(km.res)

