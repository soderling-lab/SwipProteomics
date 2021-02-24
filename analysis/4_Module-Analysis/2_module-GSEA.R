#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: analysis of modules for GSE
# author: Tyler W Bradshaw

## ---- Input parameters

BF_alpha <- 0.05 
FE_threshold <- 2
save_results <- TRUE


## ---- Set-up the workspace 

# project root dir
root <- "~/projects/SwipProteomics"

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(geneLists) # for gene lists (pathways) and hyperTest function
})

# load functions in root/R and data in root/data
devtools::load_all(root, quiet = TRUE)

# load the data from root/data
data(sig_modules) 
data(swip_gene_map) # gene_map
data(swip_partition) # partition
data(wash_interactome) # == wash iBioID sig_prots
data(msstats_sig_prots) # sig_prots

# Project Directories
datadir <- file.path(root, "data")
rdatdir <- file.path(root, "rdata")
tabsdir <- file.path(root, "tables")

# Load the gene lists from twesleyb/geneLists
# FIXME: geneLists should be associated with Pubmed IDs!
data(list = "ePSD") # Uezu et al., 2016 [4]
data(list = "iPSD") # Uezu et al., 2016 [5]
data(list = "corum") # CORUM protein complexes [1]
data(list = "takamori2006SV") # presynaptic proteome from Takamori et al. [3]
data(list = "uniprotSubcell") # uniprot subcellular anno for network prots
data(list = "lopitDCpredictions") # protein predicted subcellular locations [2]
data(list = "itzhak2017") # protein predicted subcellular locations


## ---- Do work 

# Add retriever complex
# Retriever complex from McNally et al., 2017. [6]
retriever <- c("Vps35l", "Vps26c", "Vps29")

# get entrez ids for retriever prots
names(retriever) <- gene_map$entrez[match(retriever, gene_map$symbol)]

# Collect list of modules, map Uniprot accession to Entrez.
modules <- split(names(partition), partition)[-1]
module_entrez <- lapply(modules, function(x) {
  gene_map$entrez[match(x, gene_map$uniprot)]
})
names(module_entrez) <- paste0("M", names(module_entrez))
all_entrez <- unlist(module_entrez, use.names = FALSE)

# Collect WASH BioID genes
wash_prots <- wash_interactome
wash_genes <- na.omit(gene_map$entrez[match(wash_prots, gene_map$uniprot)])

# Clean-up lopit dc predictions
names(lopitDCpredictions) <- paste("LopitDC:", names(lopitDCpredictions))

# Clean-up corum names
names(corum) <- paste("CORUM:", names(corum))

# Clean-up Takamori et al. pathway names
names(takamori2006SV) <- paste("Takamori et al., 2006:", names(takamori2006SV))

# Clean-up iPSD names
names(iPSD) <- paste("Uezu et al., 2016:", names(iPSD))

# Clean-up ePSD names
names(ePSD) <- paste("Uezu et al., 2016:", names(ePSD))

# clean-up itzhak names
names(itzhak2017) <- paste("Itzhak et al., 2017:", names(itzhak2017))

# Collect list of entrez ids for pathways of interest
gene_lists <- c(
  list("WASH-iBioID" = wash_genes), # 1
  corum, # 2
  lopitDCpredictions, # 3 and 4_
  list("McNally et al., 2017: Retriever Complex" = names(retriever)),
  takamori2006SV, # 5
  iPSD, # 6
  ePSD, # 7
  uniprotSubcell,
  list("SigProts" = sig_prots),
  itzhak2017
)

# Remove lists with less than 3 proteins
part_entrez <- setNames(partition,
			nm=mapID(names(partition),"uniprot","entrez"))
idx <- which(sapply(gene_lists, function(x) sum(x %in% names(part_entrez))<3))
gene_lists <- gene_lists[-idx]

# these are the pathways
#sapply(gene_lists,length) %>% knitr::kable()

# Loop to perform GSE for each pathway
message("\nPerforming GSE analysis for all modules:")
results <- list()
pbar <- txtProgressBar(max = length(gene_lists), style = 3)
for (experiment in names(gene_lists)) {
  # Get pathway specific genes
  pathway_genes <- gene_lists[[experiment]]
  # Background is union of all network genes and pathway genes
  background <- unique(c(all_entrez, pathway_genes))
  # Loop to perform hypergeometric test for enrichment
  results_list <- list()
  for (i in c(1:length(module_entrez))) {
    results_list[[i]] <- hyperTest(
      pathway_genes,
      module_entrez[[i]],
      background
    )
  }
  names(results_list) <- paste0("M", names(modules))
  # Collect results in a data.table
  hyper_dt <- as.data.table(do.call(rbind, results_list),
    keep.rownames = "Module"
  )
  # Adjust p-values
  hyper_dt$FDR <- p.adjust(hyper_dt$"P-value", method = "BH")
  hyper_dt$Padjust <- p.adjust(hyper_dt$"P-value", method = "bonferroni")
  # Add module size annotation
  sizes <- sapply(module_entrez, length)
  hyper_dt <- tibble::add_column(hyper_dt,
    "Module Size" = sizes,
    .after = "Module"
  )
  # Add pathway annotation
  hyper_dt <- tibble::add_column(hyper_dt,
    Pathway = experiment,
    .after = "Module Size"
  )
  # Add total number of pathway genes
  hyper_dt <- tibble::add_column(hyper_dt,
    "Total Pathway Genes" = length(pathway_genes),
    .after = "Pathway"
  )
  hyper_dt <- tibble::add_column(hyper_dt,
    "N Pathway Genes" = sum(pathway_genes %in% all_entrez),
    .after = "Total Pathway Genes"
  )
  # Number of pathway genes in each module
  n <- sapply(module_entrez, function(x) sum(pathway_genes %in% x))
  hyper_dt <- tibble::add_column(hyper_dt,
    "n Pathway Genes in Module" = n,
    .after = "N Pathway Genes"
  )
  # Sort by fold enrichment
  hyper_dt <- hyper_dt %>% arrange(desc(`Fold enrichment`))
  Entrez <- sapply(module_entrez[hyper_dt$Module], function(x) x[x %in% as.numeric(gene_lists[[experiment]])])
  hyper_dt$Genes <- lapply(Entrez, function(x) {
    paste(gene_map$symbol[match(x, gene_map$entrez)], collapse = "; ")
  })
  # Return the results
  results[[experiment]] <- hyper_dt
  setTxtProgressBar(pbar, value = match(experiment, names(gene_lists)))
}
close(pbar)


## ---- Collect the results in a single data.table

# collect results
dt <- dplyr::bind_rows(results)

# only sig + enriched results:
sig_dt <- dt %>% filter(Padjust < BF_alpha) %>% 
	filter(`Fold enrichment` > FE_threshold) 


## ---- status

m <- length(unique(sig_dt$Module))
M <- length(modules)
message(m, " of ", M, " modules exhibit some significant GSE.")

# modules with sig lopitDC enrichment
message("Modules enriched for LopitDC subcellular compartments:")
idx <- grepl("LopitDC", sig_dt$Pathway)
sig_dt %>% filter(idx) %>% 
	select(Module, Pathway, Padjust, `Fold enrichment`) %>% 
	knitr::kable()

# modules with uniprot enrichment
message("UniProt enriched: ")
sig_dt %>% filter(grepl("Uniprot",Pathway)) %>%
	select(Module, Pathway, Padjust, `Fold enrichment`) %>% 
	knitr::kable()

message("Sig Modules: ")
sig_dt %>% 
	subset(Module %in% sig_modules) %>%
	select(Module, Pathway, Padjust, `Fold enrichment`) %>% 
	arrange(as.numeric(gsub("M","",Module))) %>%
	knitr::kable()

message("NS Modules: ")
sig_dt %>% 
	subset(Module %notin% sig_modules) %>%
	select(Module, Pathway, Padjust, `Fold enrichment`) %>% 
	arrange(as.numeric(gsub("M","",Module))) %>%
	knitr::kable()


## ---- save results

# save gene lists
myfile <- file.path(root,"data","gene_lists.rda")
save(gene_lists, file=myfile,version=2)
message("saved: ", myfile)

# save as rda
module_gsea <- sig_dt
myfile <- file.path(root, "data", "module_gsea.rda")
save(module_gsea, file = myfile, version = 2)
message("saved: ", myfile)

idx <- order(as.numeric(gsub("M","",sig_dt$Module)))
sig_dt <- sig_dt[idx,] # sort by module
tmp_df <- data.table(Pathway=names(gene_lists),
	  Entrez = sapply(gene_lists,paste,collapse=";"))
tmp_list <- list("Module GSEA" = sig_dt,"Pathways" = tmp_df)

# save as excel
myfile <- file.path(root,"tables", "SWIP-TMT-Module-GSEA.xlsx")
write_excel(tmp_list,myfile)
message("saved :", myfile)
