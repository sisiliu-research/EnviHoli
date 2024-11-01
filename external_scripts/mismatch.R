#!/usr/bin/env Rscript

#==== what? ====
# 1. combine all *.mismatches.txt.gz
# save to RData

#==== load packages ====
load.packs=c("readr", "dplyr")
sapply(load.packs, require, character=TRUE)

#== show input args
args = commandArgs(trailingOnly=TRUE)
message("\nInput arguments\n")
print(args)

FILES=args[1]
FILE_END=args[2]

SAMPLES=args[3]
TAXLI=args[4]

OUT_TAB=args[5]
OUT_FILE_NAME=args[6]

#==== file list ====
FILE_LIST=list.files(path = FILES, pattern = paste0("*", FILE_END), full.names = T, recursive = T)
FILE_LIST
#==== samples infor. ====
samples=read.csv(SAMPLES)

#==== load taxli ====
load(TAXLI)

#=== loop for each gz file ====
mismatch_freq <- function(x) c(x[1:4]/sum(x[1:4]),x[5:8]/sum(x[5:8]),x[9:12]/sum(x[9:12]),x[13:16]/sum(x[13:16]))
# calculate frequency for all paired nuclitides substitituion and merged all libs into one RData file
if(TRUE) {
  all_freq_li=NULL
  for (i in FILE_LIST) {
    print(basename(i))
    dfi=read.table(gzfile(i))
    names(dfi)=c("tax_id",	"direction", "position",	"AA",	"AC",	"AG",	"AT",	"CA",	"CC",
                 "CG",	"CT",	"GA",	"GC",	"GG",	"GT",	"TA",	"TC",	"TG",	"TT")
    
    nucleotide_columns=dfi[-c(1:3)]
    freq_rows=apply(nucleotide_columns, 1, mismatch_freq)
    freq_rows=as.data.frame(t(freq_rows))
    freq_rows=cbind(dfi[, c(1:3)], freq_rows)
   
    # attach sample name and infor.
    freq_rows$sample=gsub(FILE_END, "", basename(i))
    freq_rows=left_join(freq_rows, samples, by = "sample")

    # attach taxlineage
    freq_rows_li=left_join(freq_rows, taxlineage, by = "tax_id")
    # merge
    all_freq_li=rbind(all_freq_li, freq_rows_li)
  }
  save(all_freq_li, file = paste0(OUT_TAB, "/", OUT_FILE_NAME, "_MergedData.RData"))
}

