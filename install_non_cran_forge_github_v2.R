.libPaths(c("/usr/local/lib/R/site-library", "/ricopili/dependencies/R_packages", .libPaths()))
#Install (if you haven't already) the devtools package via CRAN:
install.packages(c("devtools", "rstudioapi"))
install.packages("remotes")
library("devtools"); library("rstudioapi"); library("remotes"); 

Sys.setenv(GITHUB_PAT = "") #Insert token here, but NOT in GitHub repos

install.packages("GWAtoolbox", repos="http://R-Forge.R-project.org")

#Install postgwas from GitHub via devtools:
setRepositories(ind = 1:6)
install_github("merns/postgwas")

install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

install_github("ozancinar/poolR")
library(poolR)

devtools::install_github("MathiasHarrer/dmetar")

remotes::install_github("kassambara/easyGgplot2")


#OpenMx:
source('https://vipbg.vcu.edu/vipbg/OpenMx2/software/getOpenMx.R') #Development version - non-dev version is now also on CRAN
install.packages("umx") #Prerequisite for OpenMx

#ASSET:
if (!require("BiocManager", character.only = TRUE, quietly = TRUE)){
   install.packages("BiocManager", repos=NULL)
}
BiocManager::install("ASSET", version = "devel")
BiocManager::install("ASSET")

#LAVA: 
install.packages("BiocManager")
BiocManager::install("snpStats") #Prerequisite for LAVA
remotes::install_github("josefin-werme/LAVA")

#Also, download required reference files from Github into /notebooks/dms-gwas/opt/lava_refs

#LDpred-2: https://choishingwan.github.io/PRS-Tutorial/ldpred/#installing-ldpred-2
remotes::install_github("https://github.com/privefl/bigsnpr.git")

#Also, download required reference files from Github into /notebooks/dms-gwas/opt/ldpred_refs

#SAIGE:
install.packages(c("R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "SKAT")) #Prerequisite for SAIGE
devtools::install_github("leeshawn/MetaSKAT") #Prerequisite for SAIGE
devtools::install_github("weizhouUMICH/SAIGE")

#SBayesRC:
install.packages(c("Rcpp", "data.table", "stringi", "BH",  "RcppEigen")) #Prerequisite for SBayesRC
devtools::install_github("zhilizheng/SBayesRC")

