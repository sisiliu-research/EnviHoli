#!/usr/bin/env Rscript

#==== what? ====
# 1. combine all _metaDMGout.csv_taxlineage.csv
# save to RData

#==== load packages ====
load.packs=c("readr", "dplyr")
sapply(load.packs, require, character=TRUE)

#==== current working dir ====
getwd()

#==== input args ====
args = commandArgs(trailingOnly=TRUE)
print(args)

FILES=args[1]
FILE_END=args[2]
END_IN_FILE=args[3]

SAMPLES=args[4]
TAXLI=args[5]
OUT_TAB=args[6]
OUT_FILE_NAME=args[7]


#==== list files ====
files0=list.files(path = FILES, pattern = paste0("*", FILE_END), full.names = T, recursive = T)
files0
#==== merge ====
MergedData <- do.call(rbind, lapply(files0, read.csv))

#==== remove END_IN_FILE ====
MergedData$sample=gsub(END_IN_FILE, "", MergedData$sample)

#==== attach samples infor. ====
samples=read.csv(SAMPLES)
MergedData=left_join(MergedData, samples, by = "sample")

#==== attach taxlineage
load(TAXLI)
MergedData=left_join(MergedData, taxlineage, by = "tax_id")

#==== save ====
save(MergedData, file = paste0(OUT_TAB, "/", OUT_FILE_NAME, "_MergedData.RData"))


