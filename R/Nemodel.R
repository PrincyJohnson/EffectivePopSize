# Effective Population size Estimation
# Vcf markers data - ld file
# Author: Josephine Johnson

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("gdsfmt")
# BiocManager::install("SNPRelate")
# library(gdsfmt)
# library(SNPRelate)
# library(dplyr)
# library(ggplot2)
# library(tidyr)

Nemodel<- function(bed,fam,bim,ld,cM, species_name){

  ### Build the gdsdb
  snpgdsBED2GDS(bed, fam, bim,out.gdsfn ="plink.gds")
  snpgdsSummary("plink.gds")
  genofile <- snpgdsOpen("plink.gds")
  ###  Getting the PCA object
  pca <- snpgdsPCA(genofile, num.thread=2,autosome.only=FALSE)
  ### Calculate the percentage of variance explained by the first 10 PC
  pc.percent <- pca$varprop*100

  table <- data.frame(sample.id = pca$sample.id,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],
                      EV4 = pca$eigenvect[,4],
                      EV5 = pca$eigenvect[,5],
                      stringsAsFactors = FALSE)
  colnames(table) <- c("ID", "EV1", "EV2", "EV3", "EV4", "EV5")

  write.table(table, file = "Eigenvalues.txt", row.names = F, col.names = F, sep = "\t", quote = F)


  ggplot(data = table) +
    geom_point(mapping = aes(x = EV1,y = EV2), color = 'darkblue', show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype="dotted") +
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(title = "PCA for SNPS in Depth 4",
         x = paste0("Principle Component 1 (",round(pc.percent[1])," %)"),
         y = paste0("Principle Component 2 (", round(pc.percent[2])," %)")) +
    theme_minimal()

  ggsave('PCA1.png',width=4, height=4,dpi=300)


  ggplot(data = table) +
    geom_point(mapping = aes(x = EV1,y = EV2), color = 'darkblue', show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype="dotted") +
    geom_vline(xintercept = 0, linetype="dotted") +
    labs(title = "PCA for SNPS in Depth 4",
         x = paste0("Principle Component 2 (",round(pc.percent[2])," %)"),
         y = paste0("Principle Component 3 (", round(pc.percent[3])," %)")) +
    theme_minimal()

  ggsave('PCA2.png',width=4, height=4,dpi=300)

  ### Read the ld file
  ld <- read.table(paste(ld,'.ld',sep=""),header=T,check.names=F,stringsAsFactors=F)

  ### calculation of basepairs using the position info from ld file

  bp<-(ld$BP_B-ld$BP_A)

  ### Calculation of distance(kb)

  Phy_dist<-bp/1000
  r2 <- ld$R2
  ld$Species<- species_name

  ### Plot of Linkage disequilibrium
  qplot(Phy_dist, r2, data=ld, ylab=bquote('LD ('*r^2*')'), xlab='distance (Kb)',
        group =Species, colour =Species, alpha=I(1/256), geom=c('smooth'), span=50000+
          theme(legend.justification=c(1,1), legend.position=c(1,1))) + theme_bw()
  ggsave('LDPlot.png',width=4, height=2,dpi=300)



  ### Calculate morgans
  morgans<- Phy_dist/(cM*100)
  Na_omit_data<- na.omit(data.frame(r2= r2, distance=morgans))  ###Removing Na's

  #Calculation of mean r2 for distance cM
  Meanr2 <- Na_omit_data %>% group_by(distance) %>% summarize(mean(r2)) ## average of r^2

  Meanr2$Species <- species_name

  ### Plot for mean r2
  qplot(distance, `mean(r2)`, data=Meanr2, ylab=bquote('Mean ('*r^2*')'), xlab='Recombination frequency (Morgans)',
        group =Species, colour =Species, alpha=I(1/256), geom=c('smooth'), span=50000+
          theme(legend.justification=c(1,1), legend.position=c(1,1))) + theme_bw()
  ggsave('Meanr2vsRecombination.png',width=4, height=2,dpi=300)

  ### To predict expected (r2), create a model with mean r2 and distance
  X<-as.matrix(cbind(1,Meanr2$distance))
  beta.hat<-solve(t(X)%*%X) %*% t(X) %*% Meanr2$`mean(r2)`
  mu.hat<- X %*% beta.hat  ###prediction

  ### Fit Sved formula using linear regression - remove intercept, method inspired from Juma et al., (2021) estimation of Ne in rice
  response<- (1/mu.hat)-1
  predictor<- Meanr2$distance*4

  beta_hat_Ne<-solve(t(predictor)%*%predictor) %*% t(predictor) %*% response

  ###Returning Ne estimate
  return(Ne=beta_hat_Ne)
}
