#!/usr/bin/env Rscript 

# title: SwipProteomics
# description: generate protein co-variation (correlation) network
# author: twab

## ---- Input:
root <- "~/projects/SwipProteomics"
input_adjm <- file.path(root,"rdata","adjm.rda")


# species from which to collect PPIs for PPI graph
os_keep = c(9606, 10116, 10090)


## ---- Output:
# ppi_adjm.rda

# NOTE: large ouput files saved in root/rdata bc too big to be tracked by git

stopifnot(file.exists(input_adjm))


## ---- prepare the working environment

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load data in root/data
data(gene_map)

# load data in root/rdata
load(input_adjm)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(getPPIs) # twesleyb/getPPIs
  library(data.table)
})


# load mouse PPIs compiled from HitPredict
data(musInteractome) 


## ---- create ppi network

# NOTE: PPIs are NOT used to identify communities

# map uniprot to entrez
uniprot <- colnames(adjm)
idx <- match(uniprot,gene_map$uniprot)
entrez <- gene_map$entrez[idx]

# given entrez, collect ppis from musInteractome
ppi_df <- musInteractome %>% 
	subset(osEntrezA %in% entrez & osEntrezB %in% entrez) %>% 
	subset(Interactor_A_Taxonomy %in% os_keep) %>%
	subset(Interactor_B_Taxonomy %in% os_keep) %>%
	dplyr::select(osEntrezA, osEntrezB)

# map back to uniprot and cast to matrix for igraph
idx <- match(ppi_df$osEntrezA,gene_map$entrez)
idy <- match(ppi_df$osEntrezB,gene_map$entrez)
ppi_dm <- ppi_df %>% 
	dplyr::mutate(ProtA = gene_map$uniprot[idx], ProtB = gene_map$uniprot[idy]) %>%
	dplyr::select(ProtA,ProtB) %>% 
	as.matrix()

# create igraph graph
g <- igraph::graph_from_edgelist(ppi_dm, directed=FALSE)

# simplify (weight=0,1) and get the adjacency matrix
ppi_adjm <- as.matrix(igraph::as_adjacency_matrix(igraph::simplify(g)))

# collect proteins that are missing
missing_prots <- colnames(adjm)[colnames(adjm) %notin% colnames(ppi_adjm)]
# these proteins are unconnected^, but we include them in the ppi_adjm so the
# networks have matching vertex sets

# add missing cols
tmp_cols <- matrix(0, nrow=nrow(ppi_adjm),ncol=length(missing_prots))
colnames(tmp_cols) <- missing_prots
rownames(tmp_cols) <- rownames(ppi_adjm)
tmp_dm <- cbind(ppi_adjm,tmp_cols)

# add missing rows
tmp_rows <- matrix(0, nrow=length(missing_prots),ncol=ncol(tmp_dm))
colnames(tmp_rows) <- colnames(tmp_dm)
rownames(tmp_rows) <- missing_prots

# full(ppi)_adjm
full_adjm <- rbind(tmp_dm,tmp_rows)

# sort rows and cols to match adjm
ppi_adjm <- full_adjm[colnames(adjm),rownames(adjm)]


## ---- save networks

# coerce to data.table and write to file
ppi_dt <- as.data.table(ppi_adjm,keep.rownames="Protein")
myfile <- file.path(root,"rdata","ppi_adjm.csv")
data.table::fwrite(ppi_dt, myfile)
message("saved: ", myfile)


## ---- save as rda

# ppi adjm
myfile <- file.path(root,"rdata","ppi_adjm.rda")
save(ppi_adjm, file=myfile,version=2)
message("saved: ", myfile)
