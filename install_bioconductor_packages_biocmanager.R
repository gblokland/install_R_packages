.libPaths(c("/usr/local/lib/R/site-library", "/ricopili/dependencies/R_packages", .libPaths()))
#install.packages("BiocManager")
library(BiocManager)
#Installing package(s) ‘Biobase’, ‘IRanges’, ‘AnnotationDbi’
#also installing the dependencies ‘BiocGenerics’, ‘S4Vectors’

# Install specific packages, e.g., “GenomicFeatures” and “AnnotationDbi”, with
BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "GO.db"))
#also installing the dependencies ‘matrixStats’, ‘snow’, ‘DelayedArray’, ‘BiocParallel’, ‘SummarizedExperiment’, ‘GenomeInfoDbData’, ‘zlibbioc’, ‘XML’, ‘Rsamtools’, ‘GenomicAlignments’, ‘GenomeInfoDb’, ‘GenomicRanges’, ‘XVector’, ‘Biostrings’, ‘rtracklayer’, ‘biomaRt’
BiocManager::install(c("SpidermiR"))
#also installing the dependencies ‘ggsci’, ‘ggsignif’, ‘exactRankTests’, ‘KMsurv’, ‘km.ci’, ‘hwriter’, ‘geneplotter’, ‘rjson’, ‘ggpubr’, ‘maxstat’, ‘survMisc’, ‘cmprsk’, ‘ShortRead’, ‘DESeq’, ‘aroma.light’, ‘locfit’, ‘GetoptLong’, ‘R.oo’, ‘R.methodsS3’, ‘annotate’, ‘ALL’, ‘gsubfn’, ‘proto’, ‘downloader’, ‘survminer’, ‘EDASeq’, ‘edgeR’, ‘rvest’, ‘ComplexHeatmap’, ‘R.utils’, ‘genefilter’, ‘ConsensusClusterPlus’, ‘doParallel’, ‘sva’, ‘limma’, ‘ggrepel’, ‘sqldf’, ‘networkD3’, ‘TCGAbiolinks’, ‘miRNAtap’, ‘miRNAtap.db’, ‘org.Hs.eg.db’
BiocManager::install(c("STRINGdb"))
#also installing the dependencies ‘hash’, ‘plotrix’
BiocManager::install(c("hierGWAS"))
#also installing the dependencies ‘fastcluster’, ‘glmnet’, ‘fmsb’
BiocManager::install(c("UniProt.ws"))
#also installing the dependencies ‘dbplyr’, ‘BiocFileCache’, ‘rappdirs’
BiocManager::install(c("KEGGREST"))
BiocManager::install(c("KEGGgraph"))
#also installing the dependency ‘graph’
BiocManager::install(c("clusterProfiler"))
#also installing the dependencies ‘triebeard’, ‘RcppArmadillo’, ‘fastmatch’, ‘urltools’, ‘gridGraphics’, ‘tidygraph’, ‘graphlayouts’, ‘DO.db’, ‘fgsea’, ‘europepmc’, ‘ggplotify’, ‘ggraph’, ‘DOSE’, ‘enrichplot’, ‘GOSemSim’, ‘qvalue’, ‘rvcheck’
BiocManager::install(c("org.Mm.eg.db"))
BiocManager::install(c("org.Rn.eg.db"))
BiocManager::install(c("org.Hs.eg.db"))
BiocManager::install(c("Rgraphviz"))
BiocManager::install(c("affy"))
#also installing the dependencies ‘affyio’, ‘preprocessCore’
BiocManager::install(c("minet"))
#also installing the dependency ‘infotheo’
BiocManager::install(c("pathview"))
BiocManager::install(c("DelayedArray"))
BiocManager::install(c("RBGL"))
BiocManager::install(c("supraHex"))
#also installing the dependency ‘hexbin’
BiocManager::install(c("c3net"))
BiocManager::install(c("dnet"))
BiocManager::install(c("mygene"))
###biocLite("fifer")
BiocManager::install(c("NMF"))
#also installing the dependencies ‘aplpack’, ‘stringi’
BiocManager::install(c("mglR"))
BiocManager::install(c("grasp2db"))
BiocManager::install(c("AnnotationFilter"))
BiocManager::install(c("EnsDb.Hsapiens.v75", "ensembldb"))
#also installing the dependency ‘ProtGenerics’
BiocManager::install(c("GeneticsDesign"))

BiocManager::install(version = '3.10') #R3.6.0
BiocManager::install("ASSET", version = "devel")
#also installing the dependency ‘msm’
BiocManager::install(c("genefu"))
#also installing lots of dependencies

setRepositories() #select Bioconductor repos

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = '3.5') #R3.4.0
BiocManager::install("genefu")

BiocManager::install(c("GWASTools"))
# Bioconductor version 3.10 (BiocManager 1.30.10), R 3.6.2 (2019-12-12)
# Installing package(s) 'GWASTools'
# also installing the dependencies ‘gdsfmt’, ‘GWASExactHW’, ‘DNAcopy’, ‘logistf’, ‘quantsmooth’
