#Show where libraries are located:
.libPaths()
.libPaths( c("/usr/local/lib/R/site-library", "/ricopili/dependencies/R_packages", .libPaths() ) )
options('repos'="https://cran.rstudio.com/")
options('repos')
#wget https://cran.r-project.org/src/contrib/Archive/calibrate/calibrate_1.7.2.tar.gz
#install.packages("calibrate_1.7.2.tar.gz", repos=NULL, type="source")

#package ‘d3heatmap’ is not available (for R version 4.0.2)
#package ‘GWASTools’ is not available (for R version 4.0.2)
#"aplpack", "Rcmdr", "rgl", "tidyverse",
x <- c("devtools","calibrate", "car", "CARBayes", "data.table", "doBy", "ez","FactoMineR", 
    "gdtools", "ggmap", "ggpubr", "git2r", "googleAuthR", "haven", "heatmaply","Hmisc", "units", "xml2", 
    "labelled", "latticeExtra", "leaflet", "MBESS", "mediation", "memisc", "mitml", "mvoutlier", "plinkQC", "plotly", 
    "png", "psychmeta", "questionr", "randomcoloR", "RcmdrMisc", "RefManageR", "rgdal", "rgeos", "rgexf", 
    "RgoogleMaps", "RGraphics", "rio", "RMySQL", "robCompositions", "roxygen2", "sf", "spdep", "spelling", "splitstackshape", 
    "survminer", "survMisc", "svglite", "tcltk2", "ukbtools", "usethis", "V8", "vdiffr", "VIM", "writexl", "XML",
    "abind", "acepack", "apcluster", "ape", "arm", "askpass", "assertthat", "available", "backports", 
       "base64enc", "bayesplot", "bdsmatrix", "beeswarm", "BH", "bibtex", "bindr", "bindrcpp", "bit", 
       "bit64", "bitops", "blob", "boot", "brew", "broom", "BSDA", "calibrate", "callr", 
       "CARBayesdata", "carData", "caTools", "cellranger", "checkmate", "chron", "circlize", "classInt", "cli", 
       "clipr", "clisymbols", "cmprsk", "coda", "coin", "colorfulVennPlot", "colorspace", "colourpicker", 
       "commonmark", "compute.es", "corpcor", "corrgram", "corrplot", "covr", "cowplot", "coxme", "crayon", 
       "crosstalk", "cubature", "curl", "cvTools", "DBI", "deldir", "dendextend", 
       "DEoptimR", "desc", "DescTools", "diagram", "DiagrammeR", "diffobj", "digest", "diptest", 
       "doParallel", "dotCall64", "downloader", "dplyr", "DT", "dygraphs", "e1071", "effects", 
       "effsize", "ellipse", "ellipsis", "emmeans", "estimability", "evaluate", "exactRankTests", "expm",  
       "fansi", "farver", "fields", "flashClust", "flexmix", "fontBitstreamVera", 
       "fontLiberation", "fontquiver", "forcats", "foreach", "foreign", "forestplot", "formatR", "Formula", 
       "fpc", "freetypeharfbuzz", "fs", "futile.logger", "futile.options", "gdata", "gee", "geepack", 
       "generics", "geosphere", "getopt", "GetoptLong", "GGally", "ggbeeswarm", "ggcorrplot", "ggforce", 
       "ggplot2", "ggplot2movies", "ggrepel", "ggridges", "ggsci", "ggsignif", "ggstance", 
       "ggthemes", "gh", "GlobalOptions", "glue", "gmm", "gmodels", "GPArotation", 
       "gplots", "gridBase", "gridExtra", "gsubfn", "gtable", "gtools",  
       "hexbin", "highr", "hms", "htmlTable", "htmltools", "htmlwidgets", "httpuv", "httr", "hunspell", 
       "igraph", "igraph", "influenceR", "ini", "inline", "interactionTest", "interplot", "iotools", "irlba", 
       "iterators", "jomo", "jpeg", "jsonlite", "kernlab", "km.ci", "KMsurv", "knitr", "knitr", "labeling", 
       "laeken", "lambda.r", "later", "lavaan", "lazyeval", "leaps", 
       "LearnBayes", "libcoin", "lme4", "lmerTest", "lmtest", "loo", "lpSolve", "lubridate", "magrittr", 
       "manipulate", "manipulateWidget", "mapdata", "mapplots", "mapproj", "maps", "maptools", "markdown", 
       "MASS", "matlab", "Matrix", "matrixcalc", "MatrixModels", "matrixStats", "maxstat", "mclust", 
       "mcmc", "MCMCglmm", "MCMCpack", "memoise", "meta", "MetABEL", "metacart", 
       "metafor", "mi", "mice", "migest", "mime", "miniUI", "minqa", "mitools", "mnormt", "modeltools", 
       "multcomp", "multcompView", "MuMIn", "munsell", "mvtnorm", "nFactors", "nloptr", "nor1mix", 
       "norm", "nortest", "numDeriv", "openssl", "openxlsx", "optparse", "ordinal", "packrat", "pan", 
       "parmigene", "party", "pastecs", "pbivnorm", "pbkrtest", "pcaPP", "pheatmap", "pillar", "pkgbuild", 
       "pkgconfig", "pkgload", "pkgmaker", "plogr", "plotrix", "pls", "plyr",
        "polyclip", "polycor", "polynom", "prabclus", "praise", "prettyunits", "processx", "profvis", "progress", 
       "promises", "proto", "ps", "pscl", "psych", "purrr", "pvclust", "qap", "qqman", "qtl", 
       "quadprog", "quantreg", "R.methodsS3", "R.oo", "R.utils", "R2HTML", "R2WinBUGS", "R6", 
       "randomForest", "ranger", "raster", "RCircos", "rcmdcheck",  
       "RColorBrewer", "Rcpp", "RcppEigen", "RCurl", "readr", "readstata13", "readxl", "registry", 
       "relimp", "rematch", "remotes", "repr", "reshape", "reshape2", "rex",  
       "rjson", "rlang", "rmarkdown", "rmeta", "rngtools", "robustbase", "Rook", "rpart", "rpf", "rprojroot", "rrcov", "rrcovNA", 
       "rsconnect", "rsm", "RSQLite", "rstan", "rstanarm", "rstantools", "rstudioapi", "rtf", "Rtsne", 
       "rworldmap", "rworldxtra", "sandwich", "scales", "scatterplot3d", "selectr", "sem", "semTools", 
       "seriation", "sessioninfo", "sfsmisc", "sgeostat", "shape", "shapefiles", "shiny", "shinyjs", 
       "shinystan", "shinythemes", "snow", "snowfall", "sourcetools", "sp", "spam", "SparseM", "spData", 
       "sqldf", "sROC", "StanHeaders", "stringi", "stringr", "strucchange", "survey", "survival", "sys",  
       "TeachingDemos", "tensorA", "testthat", "TH.data", "threejs", "tibble", "tidyr", "tidyselect", 
       "tinytex", "tmvtnorm", "trimcluster", "truncnorm", "TSP", "tweenr", "ucminf",  
       "UpSetR", "utf8", "vcd", "vctrs", "VennDiagram", "vipor", 
       "viridis", "viridisLite", "visNetwork", "vroom", "webshot", "whisker", "withr", "wordcloud",  
       "WriteXLS", "xfun", "xopen", "xtable", "xts", "yaml", "zCompositions", "zeallot", "zip", "zoo")
#install.packages(x, repos="https://cran.rstudio.com/")
#for (pkg in x) {
#    writeLines(pkg)
#    if (!require(pkg,character.only = TRUE)) {
#        install.packages(pkg, dependencies=T, repos="https://cran.rstudio.com/")
#        suppressPackageStartupMessages(library(pkg))
#        #if(!require(pkg,character.only = TRUE)) stop("Package not found")
#    }
#}

#lapply(x, library, character.only = T, warn.conflicts = F, quietly = T)
#lapply(x, library, character.only = T, logical.return = FALSE, warn.conflicts = F, quietly = T, verbose=F)

# Hmisc loads "lattice" as well
# install.packages("available")
# #also installing the dependencies ‘ISOcodes’, ‘tokenizers’, ‘janeaustenr’, ‘stopwords’, ‘tidytext’, ‘SnowballC’, ‘stringdist’, ‘udapi’, ‘yesno’
# library(available)
# available(x)

library(tools) #base package so no need to install, but still need to load
library(graphics) #base package so no need to install, but still need to load
library(stats) #base package so no need to install, but still need to load
library(utils) #base package so no need to install, but still need to load
library(methods) #base package so no need to install, but still need to load
sort((.packages())) #which packages are loaded
#"snpRF" deprecated
#"xlsx","xlsxjars" doesn't work anymore because installation require rJava
#"metacor": required package "gsl" not installed properly, so metacor not installed
# package 'ggplot2' dependency: "vctrs"
# package 'biomaRt' required by "mglR" could not be found
# "e1071" --> not properly installed -  "VIM" and "mvoutlier" depend on this
# "fifer" --> not properly installed - downloaded source from cran archive
#install.packages("Rcmdr") #--> opens a pop-up window --> needs Xquartz and Xcode
#, "XLConnect", "rJava", "xlsx"
#x <- c("mvtnorm", "corpcor") # --> one of these packages deactivates the 'trim' function
#x <- c("NLP", "tm") # --> Package NLP also has an annotate function, and because you loaded that after ggplot2, NLP::annotate ended up masking ggplot2::annotate
#Package ‘mglR’ was removed from the CRAN repository.



# # Loading required package: class
# # Failed with error:  ‘lazy-load database '/Library/Frameworks/R.framework/Versions/3.6/Resources/library/class/R/class.rdb' is corrupt’
# #"class", "cluster", "fifer", "gclus", "gsl", "mgcv", "BRugs", "rjags", "R2jags", "miscF", "NADA", "nlme", "onyxR", , "rJava",  

# #GitHub packages:
# #packages ‘TwoSampleMR’, ‘poolR’ are not available in CRAN
# x <- c("TwoSampleMR", "poolR", "OpenMx", "umx", "GWAtoolbox")

# for (pkg in x) {
# pkgTest(pkg)
# }

# #Bioconductor packages
# x <- c("AnnotationDbi", "AnnotationFilter", "Biobase", "BiocGenerics", "BiocManager", "BiocParallel", "BiocVersion", "biomaRt", "Biostrings", "ComplexHeatmap", "DelayedArray", "EnsDb.Hsapiens.v75", "ensembldb", "GenomeInfoDb", "GenomeInfoDbData", "GenomicAlignments", "GenomicFeatures", "GenomicRanges", "IRanges", "limma", "mygene", "org.Hs.eg.db", "org.Mm.eg.db", "ProtGenerics", "Rsamtools", "rtracklayer", "S4Vectors", "SummarizedExperiment", "XVector", "zlibbioc")

# for (pkg in x) {
# pkgTest(pkg)
# }



#update.packages(x)
# remove.packages("rJava")
# remove.packages("xlsx")
# install.packages("xlsx") #also installs rJava

#list of loaded packages:
print(sort((.packages())))
#sort(tolower((.packages())))

asNumeric <- function(x) as.numeric(as.character(x))
asCharacter <- function(x) as.character(x)
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))
factorsCharacter <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asCharacter))
se <- function(x) sqrt(var(x)/length(x))
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
neg <- function(x) -x
count_nas <- function(x) { return(sum(is.na(x)))  } 
count_not_nas <- function(x) { return(sum(!is.na(x)))  }
count_ones <- function(x) { return(sum(x==1))  }
# You can use the pskill function in the R "tools" package to interrupt the current process and return to the console. Concretely, I have the following function defined in a startup file that I source at the beginning of each script. You can also copy it directly at the start of your code, however. Then insert halt() at any point in your code to stop script execution on the fly. This function works well on GNU/Linux and judging from the R documentation, it should also work on Windows (but I didn't check).
# halt: interrupts the current R process; a short idle time prevents R from
# outputting further results before the SIGINT (= Ctrl-C) signal is received 
halt <- function(hint = "Process stopped.\n") {
    writeLines(hint)
    require(tools, quietly = TRUE)
    processId <- Sys.getpid() 
    pskill(processId, SIGINT)
    idleTime <- 1.00
    Sys.sleep(idleTime)
}


#lsmeans is being deprecated and replaced with emmeans

#detach("package:memisc", character.only=TRUE) #memisc loads package "MASS" as well; one of these packages makes TukeyHSD not work
#detach("package:MASS", character.only=TRUE)

###FIFER package has been archived on cran, so need to isntall from github
###has several dependencies that need to be installed first
# install.packages("tidyverse")
# #also installing the dependencies ‘clipr’, ‘fs’, ‘modelr’, ‘reprex’
# install.packages("xtable")
# library(devtools)
# install_github("dustinfife/fifer")

# y <- c("gtools", "ggplot2")
# install.packages(y)
# #also installing the dependencies ‘stringi’, ‘magrittr’, ‘colorspace’, ‘Rcpp’, ‘stringr’, ‘RColorBrewer’, ‘dichromat’, ‘munsell’, ‘labeling’, ‘digest’, ‘gtable’, ‘plyr’, ‘reshape2’, ‘scales’
# y <- c("Hmisc", "psych", "foreign", "memisc", "xlsx")
# install.packages(y)
# #also installing the dependencies ‘chron’, ‘Formula’, ‘latticeExtra’, ‘acepack’, ‘gridExtra’, ‘data.table’, ‘mnormt’, ‘rJava’, ‘xlsxjars’
# y <- c("doBy", "VIM", "gdata", "car", "XLConnect", "xtable")
# install.packages(y)
# #also installing the dependencies ‘zoo’, ‘minqa’, ‘nloptr’, ‘RcppEigen’, ‘DEoptimR’, ‘lmtest’, ‘lme4’, ‘SparseM’, ‘MatrixModels’, ‘robustbase’, ‘sp’, ‘vcd’, ‘e1071’, ‘laeken’, ‘pbkrtest’, ‘quantreg’, ‘XLConnectJars’
# y <- c("diagram", "rrcovNA", "nFactors", "mvoutlier", "matlab")
# install.packages(y)
# #also installing the dependencies ‘reshape’, ‘mvtnorm’, ‘pcaPP’, ‘pls’, ‘cvTools’, ‘GGally’, ‘sROC’, ‘shape’, ‘rrcov’, ‘norm’, ‘sgeostat’, ‘robCompositions’
# y <- c("GPArotation", "mi", "mitools", "rmeta", "meta", "metafor")
# install.packages(y)
# #also installing the dependencies ‘abind’, ‘coda’, ‘arm’
# y <- c("ez", "boot", "questionr", "pastecs", "FactoMineR")
# install.packages(y)
# #also installing the dependencies ‘httpuv’, ‘mime’, ‘jsonlite’, ‘R6’, ‘assertthat’, ‘tibble’, ‘lazyeval’, ‘DBI’, ‘BH’, ‘evaluate’, ‘formatR’, ‘markdown’, ‘yaml’, ‘shiny’, ‘miniUI’, ‘rstudioapi’, ‘highr’, ‘classInt’, ‘htmltools’, ‘ellipse’, ‘flashClust’, ‘leaps’, ‘scatterplot3d’, ‘dplyr’, ‘knitr’
# y <- c("RGraphics", "devtools", "RCircos", "cowplot", "MBESS")
# install.packages(y)
# #also installing the dependencies ‘curl’, ‘openssl’, ‘pbivnorm’, ‘quadprog’, ‘StanHeaders’, ‘rpf’, ‘matrixcalc’, ‘httr’, ‘memoise’, ‘whisker’, ‘git2r’, ‘withr’, ‘gsl’, ‘lavaan’, ‘OpenMx’, ‘sem’, ‘semTools’
#also installing the dependencies ‘showtextdb’, ‘grImport2’, ‘rsvg’, ‘showtext’, ‘sysfonts’, ‘tikzDevice’
# y <- c("lubridate", "colorfulVennPlot", "gplots", "VennDiagram")
# install.packages(y)
# #also installing the dependencies ‘bitops’, ‘lambda.r’, ‘futile.options’, ‘caTools’, ‘futile.logger’
# y <- c("pvclust", "pheatmap", "mclust", "fifer")
# install.packages(y)

#also installing the dependencies ‘gargle’, ‘XML’
#also installing the dependency ‘systemfonts’

# install.packages("psychmeta")
# also installing the dependencies ‘gmm’, ‘R2Cuba’, ‘nleqslv’, ‘bibtex’, ‘rprojroot’, ‘tinytex’, ‘tmvtnorm’, ‘nor1mix’, ‘fungible’, ‘RCurl’, ‘xml2’, ‘RefManageR’, ‘rmarkdown’
# install.packages("MCMCglmm")
# also installing the dependencies ‘ape’, ‘corpcor’, ‘tensorA’, ‘cubature’
# install.packages("lsmeans")
# also installing the dependencies ‘TH.data’, ‘estimability’, ‘multcomp’
# install.packages("mediation")
# also installing the dependency ‘lpSolve’
# install.packages("ordinal")
# also installing the dependency ‘ucminf’
# install.packages("memisc")
# also installing the dependency ‘repr’
# install.packages("rworldmap")
# also installing the dependencies ‘dotCall64’, ‘spam’, ‘fields’
# install.packages("randomcoloR")
# also installing the dependencies ‘V8’, ‘Rtsne’
# install.packages("RMySQL")
# also installing the dependency ‘DBI’
# install.packages("RSQLite")
# also installing the dependencies ‘bit’, ‘bit64’, ‘blob’
# install.packages("circlize")
# also installing the dependency ‘GlobalOptions’
# install.packages("party")
# also installing the dependencies ‘strucchange’, ‘coin’
# install.packages("colourpicker")
# also installing the dependency ‘shinyjs’
# install.packages("Rcmdr")
# also installing the dependencies ‘readstata13’, ‘nortest’, ‘survey’, ‘RcmdrMisc’, ‘effects’, ‘tcltk2’, ‘relimp’
# install.packages("CARBayes")
# also installing the dependencies ‘shapefiles’, ‘LearnBayes’, ‘deldir’, ‘CARBayesdata’, ‘MCMCpack’, ‘spdep’
# install.packages("dygraphs")
# also installing the dependency ‘xts’
# install.packages("ggmap")
# also installing the dependencies ‘jpeg’, ‘geosphere’
# install.packages("pkgmaker")
# also installing the dependency ‘registry’
# install.packages("roxygen2")
# also installing the dependency ‘commonmark’
# install.packages("rstanarm")
# also installing the dependencies ‘crosstalk’, ‘PKI’, ‘RJSONIO’, ‘packrat’, ‘DT’, ‘rsconnect’, ‘shinythemes’, ‘threejs’, ‘loo’, ‘shinystan’
# install.packages("rgl")
# also installing the dependencies ‘praise’, ‘testthat’, ‘processx’, ‘webshot’, ‘manipulateWidget’, ‘magick’, ‘alphashape3d’
# install.packages("metacart")
# also installing the dependency ‘rpart.plot’ (## Loading required package: rpart ## Loading required package: rpart.plot ## Loading required package: ggplot2 ## Loading required package: gridExtra)
# install.packages("ggpubr")
# also installing the dependency ‘polynom’
#library(grid) # is a base package in newer R version, so no need to install or load
#"tools" is a base package in newer R version, so no need to install  --> contains smartbind so do need to load
#"lattice" is a base package in newer R version, so no need to install 
#package 'tools' is a base package, and should not be updated 
#package 'tcltk' is a base package, and should not be updated 

# dev_mode(on=T)
# install_github("baptiste/gtable")
#remove.packages(x, "/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

#phenotypicForest includes polarHistogram: plots a large number of histogram on a wheel, to save space.
#devtools::install_github("chrislad/phenotypicForest")
#

#http://sohu.io/questions/1416896/change-text-color-for-cells-using-tablegrob-in-r
#source("/Users/gabriella/Documents/GENUS/neuropsych/code/drawDetails.table.R")
#source("/Users/gabriella/Documents/GENUS/neuropsych/code/find_aes.R")

# getOption("defaultPackages") #automatically loaded
#[1] "datasets"  "utils"     "grDevices" "graphics"  "stats"     "methods"

#sessionInfo()

# # Install from CRAN
# install.packages("ukbtools")
# # Install latest development version
# devtools::install_github("kenhanscombe/ukbtools", build_vignettes = TRUE, dependencies = TRUE)
