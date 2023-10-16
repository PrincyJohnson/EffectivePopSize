\name{Nemodel}
\alias{Nemodel}
\title{Effective Population Size}
\usage{
Nemodel(bed, bim, fam, ld, cM, species_name)
}
\description{
This function provides PCA output, PCA plots, LD plots and Ne estimation
}
\examples{
dir.create('./usdald3') # create a folder for all results from GLM analysis
setwd('./usdald3')
bed <- "/Users/josephine.p.johnson/OneDrive - North Dakota University System/Documents/PhD/chapter1/LD_Ne/LD/vcf/plink.bed"
fam <- "/Users/josephine.p.johnson/OneDrive - North Dakota University System/Documents/PhD/chapter1/LD_Ne/LD/vcf/plink.fam"
bim <- "/Users/josephine.p.johnson/OneDrive - North Dakota University System/Documents/PhD/chapter1/LD_Ne/LD/vcf/plink.bim"
ld<-"C:/Users/josephine.p.johnson/OneDrive - North Dakota University System/Documents/PhD/chapter1/LD_Ne/LD/vcf/usdald"


Nemodel(bed=bed,fam=fam,bim=bim,ld=ld,cM=300,species_name="peas")
}
