#!/usr/bin/env Rscript

#==== what? ====
# 1. read *.lca.txt.gz (metadmg)
# 2. attach key ranks and full lineage
# 3. save to *.txt (a tab-separated ( sep = “\t”) file)

#==== load packages ====
load.packs=c("readr", "dplyr", "stringr")
sapply(load.packs, require, character=TRUE)

#==== input args ====
args = commandArgs(trailingOnly=TRUE)
message("\nInput arguments\n")
print(args)

#==== Funs ====
get_nodes <- function(taxdir) {
  # We only read the first five tab-separated fields
  n <- scan(paste(taxdir,sep="/"),what=as.list(character(5)),
            quote="",flush=TRUE,sep="\t")
  # The information we want is in fields 1, 3 and 5
  # (fields 2 and 4 contain the vertical bar)
  tax_id <- as.numeric(n[[1]])
  ncbi_rank <- n[[5]]
  nodes <- data.frame(tax_id, ncbi_rank)
  return(nodes)
}

# Define a function to split the cell values
split_cell <- function(cell_value) {
  parts <- str_split(cell_value, ":", simplify = TRUE)
  
  if (length(parts) == 2) {
    # Only one col, set rank to NA
    tax_id <- as.numeric(parts[1])
    taxon <- parts[2]
    rank <- NA
  } else {
    # Extract specific parts according to your rules
    tax_id <- as.numeric(parts[1])
    rank <- parts[length(parts)]
    taxon <- paste(parts[2:(length(parts) - 1)], collapse = ":")
  }
  return(data.frame(tax_id, taxon, rank))
}

#==== combine lineages ====
# rank
nodes=as.data.frame(get_nodes(args[6]))

#full lineage
fulllineage=as.data.frame(read_tsv(args[4],
                                   col_names = c("tax_id", "ncbi_name", "path"),
                                   col_types = ("n-c-c-")))
fulllineage$path=paste0(fulllineage$path, " ", fulllineage$ncbi_name)

# change NA root to root, NA cellular organisms to cellular organisms
fulllineage$path=ifelse(fulllineage$tax_id == 1, "root", fulllineage$path)
fulllineage$path=ifelse(fulllineage$tax_id == 131567, "cellular organisms", fulllineage$path)

# The ranks are species, genus, family, order, class, phylum, kingdom, and superkingdom.
rankedlineage=as.data.frame(read_tsv(args[5],
                                     col_names = c("tax_id", "ncbi_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "clade", "superkingdom"),
                                     col_types = ("n-c-c-c-c-c-c-c-c-c-c-"))) %>% select(-ncbi_name)
# combine
taxlineage=left_join(nodes, fulllineage, by = "tax_id")
taxlineage=left_join(taxlineage, rankedlineage, by = "tax_id")

# add infor. at 7 ranks
taxlineage$species=ifelse(taxlineage$ncbi_rank == "species", taxlineage$ncbi_name, taxlineage$species)
taxlineage$genus=ifelse(taxlineage$ncbi_rank == "genus", taxlineage$ncbi_name, taxlineage$genus)
taxlineage$family=ifelse(taxlineage$ncbi_rank == "family", taxlineage$ncbi_name, taxlineage$family)
taxlineage$order=ifelse(taxlineage$ncbi_rank == "order", taxlineage$ncbi_name, taxlineage$order)
taxlineage$class=ifelse(taxlineage$ncbi_rank == "class", taxlineage$ncbi_name, taxlineage$class)
taxlineage$phylum=ifelse(taxlineage$ncbi_rank == "phylum", taxlineage$ncbi_name, taxlineage$phylum)
taxlineage$kingdom=ifelse(taxlineage$ncbi_rank == "kingdom", taxlineage$ncbi_name, taxlineage$kingdom)
taxlineage$clade=ifelse(taxlineage$ncbi_rank == "clade", taxlineage$ncbi_name, taxlineage$clade)
taxlineage$superkingdom=ifelse(taxlineage$ncbi_rank == "superkingdom", taxlineage$ncbi_name, taxlineage$superkingdom)

#==== read lca ====
li=paste0(args[1], "/", args[2], args[3])
print(paste0("input lca.txt.gz file: ", li))

# aggregate 
DF1 = data.frame(taxon=character(), stringsAsFactors=F)
DF2.1 = read.csv(file = gzfile(li), header=F, sep="\t", stringsAsFactors=F, fill=T,
                  col.names = paste0("V",seq_len(60)), comment.char = "#")
DF2.2 = data.frame(taxa=DF2.1[,2],
                   count=rep(1, dim(DF2.1)[1]),
                   stringsAsFactors = F)
DF2.3 = aggregate(DF2.2[,2]~DF2.2[,1], data=DF2.2, FUN=sum)
colnames(DF2.3) = c("taxon", "taxonReads")
DF1 = merge(DF1, DF2.3, by="taxon", all=T)

# split
sldf=do.call(rbind, lapply(DF1$taxon, split_cell))
sldf=cbind(sldf, DF1[, 2])
names(sldf)[4]="taxonReads"

# attach taxlineage
tempdf=sldf
tempdf=left_join(tempdf, taxlineage, by = "tax_id")

# NA in path
na_in_path=tempdf[is.na(tempdf$path), ]

if (isTRUE(nrow(tempdf) == nrow(DF1)) & nrow(na_in_path) == 0) {
  # add sample
  tempdf$sample=args[2]
  # save
  write_tsv(tempdf, file = paste0(args[1], "/", args[2], gsub(".txt.gz", "", args[3]), ".taxlineage.tsv"))

} else if (isTRUE(nrow(tempdf) == nrow(DF1)) & nrow(na_in_path) > 0) {
  # add sample
  tempdf$sample=args[2]
  # save
  write_tsv(tempdf, file = paste0(args[1], "/", args[2], gsub(".txt.gz", "", args[3]), ".taxlineage.tsv"))
  # save
  write_tsv(na_in_path, file = paste0(args[1], "/", args[2], gsub(".txt.gz", "", args[3]), ".no.taxlineage.tsv"))

} else {
  message("\nerror occurs: left_join(taxlineage)\n")
}








