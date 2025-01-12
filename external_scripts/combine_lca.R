#!/usr/bin/env Rscript

#==== what? ====
# 1. combine all .lca.taxlineage.tsv
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

SAMPLES=args[3]

OUT_TAB=args[4]
OUT_FILE_NAME=args[5]


#==== list files ====
files0=list.files(path = FILES, pattern = paste0("*", FILE_END, "$"), full.names = T, recursive = T)
files0

#==== merge ====
MergedData <- do.call(rbind, lapply(files0, read_tsv))

#==== attach samples information
samples=read.csv(SAMPLES)
MergedData=left_join(MergedData, samples, by = "sample")

#==== save ====
save(MergedData, file = paste0(OUT_TAB, "/", OUT_FILE_NAME, "_MergedData.RData"))

