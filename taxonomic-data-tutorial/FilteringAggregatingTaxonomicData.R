## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Load r packages---------------------------------------------------------------------------------------------------------------------------------
library(readr)
library(dplyr)



## ----Load lca RData, generated from EnviHoli pipeline------------------------------------------------------------------------------------------------
load("path/to/btoko_lca.v1_MergedData.RData")

# > names(MergedData)
# [1] "tax_id"       "taxon"        "rank"         "taxonReads"   "ncbi_rank"    "ncbi_name"    "path"         "species"      "genus"       
# [10] "family"       "order"        "class"        "phylum"       "kingdom"      "clade"        "superkingdom" "sample"       "lib_id"      
# [19] "lib_id_order" "years"        "ka"



## ----Filtering by superkingdom-----------------------------------------------------------------------------------------------------------------------
Viruses <- MergedData %>%
  filter(superkingdom %in% "Viruses")

Eukaryota <- MergedData %>%
  filter(superkingdom %in% "Eukaryota")



## ----Filtering by superkingdom and taxonomic rank----------------------------------------------------------------------------------------------------
Viruses_family <- MergedData %>%
  filter(superkingdom %in% "Viruses") %>%
  filter(!is.na(family))

Eukaryota_family <- MergedData %>%
  filter(superkingdom %in% "Eukaryota") %>%
  filter(!is.na(family))


## ----Filtering by superkingdom and taxonomic rank and sediments--------------------------------------------------------------------------------------
Viruses_family_sed <- MergedData %>%
  filter(superkingdom %in% "Viruses") %>%
  filter(!is.na(family)) %>%
  filter(years > -999)

Eukaryota_family_sed <- MergedData %>%
  filter(superkingdom %in% "Eukaryota") %>%
  filter(!is.na(family)) %>%
  filter(years > -999)



## ----Filtering by superkingdom and taxonomic rank and controls---------------------------------------------------------------------------------------
Viruses_family_contr <- MergedData %>%
  filter(superkingdom %in% "Viruses") %>%
  filter(!is.na(family)) %>%
  filter(years < -999)

Eukaryota_family_contr <- MergedData %>%
  filter(superkingdom %in% "Eukaryota") %>%
  filter(!is.na(family)) %>%
  filter(years < -999)



## ----Filtering by taxa name--------------------------------------------------------------------------------------------------------------------------
taxa_name <- c("Mammuthus", "Salix alba", "Cyanobacteriota")

lca_mytaxa <- MergedData %>%
  filter(taxon %in% taxa_name)


## ----Aggregating total counts------------------------------------------------------------------------------------------------------------------------
# total count of Viruses
sum(Viruses$taxonReads)

# total count of Viruses per libraries
Viruses_lib_count <- aggregate(.~sample, Viruses[, c("sample", "taxonReads")], FUN = sum)

# total count of Viruses per sediments
Viruses_sed_count <- aggregate(.~sample, Viruses[Viruses$years > -999, c("sample", "taxonReads")], FUN = sum)



