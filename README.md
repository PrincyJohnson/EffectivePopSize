# EffectivePopSize
This Package helps to estimate effective population size (Ne) in crops using linkage disequilibrium file (.ld)


**Install Package:**

      Library(devtools)
      
      devtools::install_github("PrincyJohnson/EffectivePopSize")
      
      Library(EffectivePopSize)



**How to use:**

Nemodel(bed, bim, fam, ld, cM, species_name)

      bed = directory of the bed file eg. "/Documents/plink.bed"
      
      bim = directory of the bim file eg. "/Documents/plink.bim"
      
      fam = directory of the fam file eg. "/Documents/plink.fam"
      
      * Use prefix 'plink' for bed, bim, and fam files*
        
      ld = directory of the ld file eg. "Documents/ldfile"
      
      cM = 1 cM value in your species (in Kb)


**This package provides PCA, LD and recombination frequency plots and finally Ne estimation**
