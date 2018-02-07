
cat('###################################################\n')
cat('###################################################\n')
cat('Connecting to bioconductor.org...\n')
source("https://bioconductor.org/biocLite.R")
cat('###################################################\n')
cat('###################################################\n')
cat('Installing the packages from Bioconductor...\n')
biocLite(c("Biobase", "lambda.r", "futile.options", "SummarizedExperiment", "futile.logger", "snow", "RCurl", "XML", "zlibbioc","GenomicAlignments", "bitops", "BiocParallel", "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "Biostrings", "rtracklayer", "XVector", "Rsamtools", "BSgenome"))
biocLite("BSgenome.Hsapiens.UCSC.hg19")

pkgs <- setdiff(c("ggplot2","plotly","stringr","shiny","shinyjs","reshape2","rhandsontable","foreach","doParallel"),installed.packages())
cat('###################################################\n')
cat('###################################################\n')
cat('Installing the packages from CRAN...\n')
install.packages(pkgs,repos="http://cran.us.r-project.org")



