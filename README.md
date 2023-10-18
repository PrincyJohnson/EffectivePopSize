# EffectivePopSize
This package is designed to facilitate the estimation of effective population size (Ne) in crop species by analyzing linkage disequilibrium data stored in .ld files.


**Install Package:**

                if (!require("BiocManager", quietly = TRUE))
                       install.packages("BiocManager")
                       
                BiocManager::install("gdsfmt")
                
                BiocManager::install("SNPRelate")
                
                library(gdsfmt)
                
                library(SNPRelate)
                
                library(devtools)
                              
                devtools::install_github("PrincyJohnson/EffectivePopSize")
                              
                library(EffectivePopSize)



**How to use:**

              Nemodel(bed, bim, fam, ld, cM, species_name)

               bed = directory of the bed file eg. "/Documents/plink.bed"
                              
               bim = directory of the bim file eg. "/Documents/plink.bim"
                              
               fam = directory of the fam file eg. "/Documents/plink.fam"
                                
               ld = directory of the ld file eg. "Documents/ldfile"
                              
               cM = 1 cM value in your species (in Kb) eg. 200

               species_name =  eg. "Peas"

**This package offers a comprehensive suite of tools for genetic analysis, including Principal Component Analysis (PCA), Linkage Disequilibrium (LD) plot, recombination frequency plot, and precise estimation of Ne (effective population size).**
