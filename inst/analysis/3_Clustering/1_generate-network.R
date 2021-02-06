#!/usr/bin/env Rscript

# author: twab
# title: SwipProteomics
# description: generate protein co-variation (correlation) network and perform
#   network enhancement

## ---- Input:
root <- "~/projects/SwipProteomics"

## ---- Output:

# * adjm.rda
# * ne_adjm.rda

# * ne_adjm.csv --> for leidenalg clustering!

# NOTE: large ouput files (>100mb) saved in root/rdata bc too big to be tracked by git


## ---- prepare the working environment

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load data in root/data
data(swip_tmt)
data(swip_gene_map)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(neten) # twesleyb/neten
  library(igraph)
  library(data.table)
})


## ---- create covariation network

message("Generating covariation network...")

# no median summarization of bioreplicates
# network is constructed from log2(Intensity) ~ Abundance
dm <- swip_tmt %>%
	reshape2::dcast(Protein ~ Mixture + Condition, value.var = "Abundance") %>%
	as.data.table() %>%
	as.matrix(rownames="Protein")


# there are a small number proteins with some missing vals
# e.g. Q9QUN7 = low abundance only quantified in 4/7 fractions
# remove these proteins
idx <- apply(dm,1,function(x) any(is.na(x)))
warning(sum(idx)," proteins with any missing values are removed.")
filt_dm <- dm[!idx,]

stopifnot(!any(filt_dm<0))

## calculate coorrelation matrix
adjm <- cor(t(filt_dm), method="pearson",use="complete.obs")


## ---- network enhancement

# Wang et al., 2018 (Nature Communications; PMID:30082777)

message("Performing network enhancement...")

# FIXME: is there any room for speed improvements here? Probably...
ne_adjm <- neten(adjm) # result is robust to neten parameters


## ---- save networks as csv in rdata

# coerce to data.table and save adjm.csv
adjm_dt <- as.data.table(adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","adjm.csv")
data.table::fwrite(adjm_dt, myfile)
message("saved: ", myfile)

# coerce to data.table and save ne_adjm.csv
ne_adjm_dt <- as.data.table(ne_adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","ne_adjm.csv")
data.table::fwrite(ne_adjm_dt, myfile)
message("saved: ", myfile)


## ---- save data as rda

# NOTE: data are saved in root/rdata
# adjm could be make smaller by melting to edge list
# but it is still too big at ~ 150 mb
# ne_adjm is smaller bc it is sparse, but still to large to be easily tracked by
# git

# adjm
myfile <- file.path(root,"rdata","adjm.rda")
save(adjm, file=myfile,version=2)
message("saved: ", myfile)

# ne adjm
myfile <- file.path(root,"rdata","ne_adjm.rda")
save(ne_adjm, file=myfile,version=2)
message("saved: ", myfile)
