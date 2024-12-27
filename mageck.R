
#3.6.1版本的R中
BiocManager::install("MAGeCKFlute")
library(MAGeCKFlute)
library(ggplot2)

#mageck，核心算法：rra
#mageck-vispr，核心算法：mle


#rra结果的后续分析流程
#gene_summary file
file1 = file.path(system.file("extdata", package = "MAGeCKFlute"),
                  "testdata/rra.gene_summary.txt")
countsummary = read.delim(file1, check.names = FALSE)
head(countsummary)
View(countsummary)

#sgrna_summary file
file2 = file.path(system.file("extdata", package = "MAGeCKFlute"),
                  "testdata/rra.sgrna_summary.txt")

FluteRRA(file1, file2, prefix ="Test", organism="hsa",
          outdir = "./")


#mle结果的后续分析流程，基于beta-score的分析与可视化
file3 = file.path(system.file("extdata", package = "MAGeCKFlute"),
                  "testdata/mle.gene_summary.txt")
countsummary = read.delim(file3, check.names = FALSE)
head(countsummary)
View(countsummary)
FluteMLE(file3, treatname="plx", ctrlname="dmso", prefix="Test", organism="hsa")




#step by step analysis
#input data
file4 = file.path(system.file("extdata", package = "MAGeCKFlute"),
                  "testdata/countsummary.txt")
countsummary = read.delim(file4, check.names = FALSE)
head(countsummary)

#gini index of sgRNAs
source("BarView.R")
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads")
#zeroCounts of sgRNAs
countsummary$Missed = log10(countsummary$Zerocounts)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")

# Read mapping
MapRatesView(countsummary)


#开始对分析结果进行下游分析



#对mageck-vispr的mle test分析结果进行可视化（其实就是mle.gene_summary.txt)
file3 = file.path("mle.gene_summary.txt")

FluteMLE(file3, treatname="plx", ctrlname="dmso", prefix="Test0521", organism="hsa")


