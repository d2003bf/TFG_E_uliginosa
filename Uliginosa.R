#Script for the FST and PCA analysis of E. uliginosa developed by Daniel Bellido
setwd()
library(adegenet)
library(vcfR)

matrix <- read.vcfR("uliginosa_snps_target_filtered_thin.recode.vcf")
matrix1 <- vcfR2genlight(matrix)

#Use either pops or popsC depending on whether you want to filter based on the location or the country
pops<- c("P-SAN", "P-SAN", "P-SAN", "P-SAN", "P-GRA", "P-GRA", "P-GRA", "P-GRA", "P-ALC", "P-ALC", "P-ALC", "P-ALC", 
         "G-SIX","G-SIX","G-SIX","G-SIX","G-SIX","G-SIX", "G-VIL", "G-VIL","G-VIL","G-VIL","G-VIL","G-AND", "G-AND", "G-AND")
popsC<- c("P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "P", "G","G","G","G","G","G", "G", "G","G","G","G","G", "G", "G")
matrix1@pop <- as.factor(pops)
matrix1@pop <- as.factor(popsC)
popNames(matrix1)

pca1 <- glPca(matrix1)
#We select the first 3 PCA's
summary(pca1)

scatter(pca1)
#We use null to ensure there is nothing underneath the legend
scatter(pca1, posi = "null")
scatter(pca1, posi="bottomright")
mycolor <- c("red","green","blue","purple","orange","#FFD300")

s.class(pca1$scores,pop(matrix1), col=mycolor, clabel = 0.6, cpoint = 1, axesell = FALSE, cellipse = FALSE, cstar=0)
add.scatter.eig(pca1$eig , 2,1,2 , posi = "bottomright")

#FST
library(dartR.base)
library(poppr)
#We create the FST
fst <- gl.fst.pop(matrix1, nboots =1, percent = 95, nclusters = 1)
write.table(fst, "./fst_uliginosa.txt") #we save the fst into a table if need be
library(reshape2)
library(ggplot2)
#We can create a visual representation of the FST with a color gradient
melted <- melt(fst, na.rm = TRUE)
ggplot(data= melted, aes(Var1, Var2, fill = value)) + geom_tile(color= "white") + 
  scale_fill_gradient(low = "white", high = "red", name= "FST")


library(ape)
#Generamos un "arbol filogenetico" simplificado de nuestras especies
uli.tree <- nj(fst)
#for some reason the tree does not place the colors in the same order as the PCA
#We need to change the color to properly match our other figures
mycolor1 <- c("purple","orange","#FFD300","green","blue","red")
plot(uli.tree, type="unr", tip.col =mycolor1, font = 1)
annot <- round(uli.tree$edge.length, 2)
edgelabels(annot[annot>0], which(annot>0), frame= 'n')
#While we can include edgelabels for the distances it doesn't tend to look good in my 
#opinion so i opt for a scale bar
add.scale.bar()

