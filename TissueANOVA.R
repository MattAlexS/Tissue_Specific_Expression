library(ggplot2)

data <- read.csv(file = file.choose())

tissues <- unique(data$Group)

x <- tissues
y <- tissues

table <- expand.grid(X=x, Y=y)

res.aov <- aov(APD ~ Group, data = data)

tuk <- TukeyHSD(res.aov)

summary(res.aov)

boncor <- choose(length(tissues),2)
z <- c()

for(i in tissues){
  for(j in tissues){
    comp <- paste(j, i, sep = "-")
    revcomp <- paste(i, j, sep = "-")
    if (i == j){
      z <- c(z, 0)
    }
    else if (comp %in% rownames(tuk$Group)){
      if (tuk$Group[comp, 4] < 0.05/boncor){
        if (tuk$Group[comp, 1] < 0){
          z <- c(z, -1)
        } 
        else {
          z <- c(z, 1)
        }
      }
      else {
        z <- c(z, 0)
      }
    }
    else if (revcomp %in% rownames(tuk$Group)){
      if (tuk$Group[revcomp, 4] < 0.05/boncor){
        if (tuk$Group[revcomp, 1] < 0){
          z <- c(z, 1)
        } 
        else {
          z <- c(z, -1)
        }
      }
      else {
        z <- c(z, 0)
      }
    }
  }
}
      
      

table$Z = z

ggplot(table, aes(X, Y, fill= Z)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu") +
  ggtitle("Tissue Compactness\nElevated Genes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 90)) +
  labs( y="Tissue 2", x="Tissue 1", fill="Relationship")
