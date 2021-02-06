#!/usr/bin/env Rscript

# title: SwipProteomics
# description: MSstatsTMT analysis of SWIP-TMT Proteomics
# author: twab

analyze_subset = FALSE

## MSstatsTMT options
MBimpute <- TRUE # impute missing val within an MS run?
rm_single <- TRUE # rm prots identified by single feature?
global_norm <- TRUE # global peptide-level sample normalization?
reference_norm <- TRUE # normalize to norm condition (SPQC) samples?
remove_norm_channel <- TRUE # rm norm condition (SPQC) data?

## threshold for protein significance
FDR_alpha = 0.05 # FDR = Benjamini Hochberg FDR

## ---- prepare the working R environment

root <- "~/projects/SwipProteomics"
renv::load(root)

#library(SwipProteomics)
devtools::load_all(root)

# load data in root/data
data(pd_psm)
data(swip_gene_map)
data(pd_annotation)
data(mut_vs_control) # 'Mutant-Control' comparison
data(msstats_contrasts) # 'intra-BioFraction' comparisons

# NOTE: msstats_contrasts is a matrix specifying pairwise contrasts between all
# 'BioFraction.Control' and 'BioFraction.Mutant' Conditions.

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  # my MSstats forks:
  library(MSstats) # twesleyb/MSstats
  library(MSstatsTMT) # twesleyb/MSstats
})

## NOTE: my fork attempts to remove much of MSstats's verbosity.
## You may find some debug.log files on your computer still, their source
## has not been found yet.
## twesleyb/MSstatsTMT includes access to the internal functions used by
## MSstatsTMT to fit protein-wise models and perform statistical comparisons.
## MSstatsTMT is a wrapper around MSstats. My fork allows you to pass arguments
## for parallel processing to MSstats::proteinSummarization to speed things up.


## ---- analyze a subset of the data


if (analyze_subset) {
  nprot = 23
  proteins <- sample(unique(pd_psm$"Master.Protein.Accessions"),nprot)
  message("Analyzing a subset of the data, ", nprot, " proteins.")
  pd_psm <- pd_psm %>% subset(`Master.Protein.Accessions` %in% proteins)
}



## ---- convert PSM data to MSstatsTMT format
# proteins with a single feature (i.e. peptide) are removed if rm_single

message("\nConverting PD PSM-level data into MSstatsTMT format.")

t0 <- Sys.time()

# aprox 7 minutes for all proteins
suppressMessages({ # verbosity
  msstats_psm <- PDtoMSstatsTMTFormat(pd_psm,
    pd_annotation, # see 0_PD-data-preprocess.R
    which.proteinid = "Master.Protein.Accessions",
    rmProtein_with1Feature = rm_single
  )
})

message("\nTime to reformat PSM data: ",
	round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")


## ---- protein-level summarization and normalization
# perform protein summarization for each run with MSstatsTMT

# NOTE: my fork allows you to pass additional args to underlying MSstats
# dataProcess function  -- speed things up by specifying the number of cores to
# be used for parallel processing.

message("\nPerforming normalization and protein summarization using MSstatsTMT.")

n_cores <- parallel::detectCores() - 1

t0 <- Sys.time()

suppressMessages({ # verbosity
  msstats_prot <- proteinSummarization(msstats_psm,
    method = "msstats",
    remove_norm_channel = remove_norm_channel,
    global_norm = global_norm, # perform global norm using 'norm' condition
    MBimpute = MBimpute,
    reference_norm = reference_norm,
    clusters = n_cores
  )
})

# This takes about 11 minutes for 8.5 k proteins with 23 cores
# FIXME: fix warnings messages about closing clusters.
proteins <- unique(as.character(msstats_prot$Protein))
message(
  "\nTime to summarize ", length(proteins), " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes."
)


## ---- perform statistical comparisons with MSstatsTMT

# NOTE: for the pairwise contrasts, MSstats fits the lmer model:
# fx <- formula(Abundance ~ 1 + Condition + (1|Mixture)) # lmerTest::lmer

# We specify Condition as Genotype.BioFraction for all intra-fraction
# comparisons. T-statistics are moderated using ebayes methods in limma.

# We perform the two analyses seperately so we can specify moderated = TRUE for
# intra-BioFraction comparisons and moderated = FALSE for overall
# 'Mutant-Control' comparison.

message("\nAssessing protein-level comparisons with MSstatsTMT.")

t0 <- Sys.time()

## 1. 'intra-BioFraction' comparisons
suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity FIXME:
    results1 <- groupComparisonTMT(msstats_prot,
      msstats_contrasts,
      moderated = TRUE
    )
  })
})

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform 8 'intra-BioFraction' comparisons for ",
  formatC(length(proteins),big.mark=","), " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")


## 2. 'Mutant-Control' comparison
# repeat for overall comparision

t0 <- Sys.time()

# contrast.matrix should be a matrix!
suppressWarnings({ # about closing clusters FIXME:
  suppressMessages({ # verbosity FIXME:
    results2 <- groupComparisonTMT(msstats_prot,
				   contrast.matrix = mut_vs_control,
				   moderated = FALSE) })
})

# This takes about 21 minutes for 8.5 k proteins
message(
  "\nTime to perform 'Mutant-Control' comparison for ",
  formatC(length(proteins),big.mark=","), " proteins: ",
  round(difftime(Sys.time(), t0, units = "min"), 3), " minutes.")


# combine intra-BioFraction and Mutant-Control results
msstats_results <- rbind(results1,results2)


## ---- format msstats_prot for downstream analysis

# clean-up the data
msstats_prot$Run <- NULL
msstats_prot$TechRepMixture <- NULL
msstats_prot$Channel <- as.character(msstats_prot$Channel)
msstats_prot$BioReplicate <- as.character(msstats_prot$BioReplicate)
msstats_prot$Condition <- as.character(msstats_prot$Condition)
msstats_prot$Mixture <- as.character(msstats_prot$Mixture)
msstats_prot$Genotype <- sapply(strsplit(msstats_prot$Condition,"\\."),"[", 1)
msstats_prot$BioFraction <- sapply(strsplit(msstats_prot$Condition,"\\."),"[", 2)
msstats_prot <- msstats_prot %>% subset(!is.na(Abundance))


# map uniprot to gene symbols and entrez ids
proteins <- unique(as.character(msstats_prot$Protein))
idx <- match(msstats_prot$Protein,gene_map$uniprot)
msstats_prot <- msstats_prot %>%
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")


## ---- format msstats_results for downstream analysis

# drop NA pvals and remove SingleMeasurePerCondition
msstats_results <- msstats_results %>%
	filter(!is.na(adj.pvalue)) %>% filter(is.na(issue)) %>% select(-issue)
colnames(msstats_results)[colnames(msstats_results) == "Label"] <- "Contrast"
colnames(msstats_results)[colnames(msstats_results) == "pvalue"] <- "Pvalue"
colnames(msstats_results)[colnames(msstats_results) == "adj.pvalue"] <- "FDR"

# annotate with gene Symbols and Entrez ids
idx <- match(msstats_results$Protein,gene_map$uniprot)
msstats_results <- msstats_results %>%
	tibble::add_column(Symbol = gene_map$symbol[idx],.after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx],.after="Symbol")

# bonferroni padjust
msstats_results <- msstats_results %>% group_by(Contrast) %>%
	mutate(Padjust=p.adjust(Pvalue,method="bonferroni"))


## ---- proteins with overall sig change

df <- subset(msstats_results,Contrast=='Mutant-Control')
sig_prots <- unique(as.character(df$Protein)[df$FDR < FDR_alpha])


## ---- save msstats_results as an excel document

# format the data for saving
results_list <- msstats_results %>% group_by(Contrast) %>% group_split()
names(results_list) <- sapply(results_list,function(x) unique(x$Contrast))
namen <- names(results_list)
new_names <- gsub("Mutant.F[0-9]{1,2}-Control\\.","",namen)
names(results_list) <- new_names

# sort by pvalue
results_list <- lapply(results_list,arrange,Pvalue)

# sort results list
idx <- c("F4","F5","F6","F7","F8","F9","F10","Mutant-Control")
results_list <- results_list[idx]

class(results_list) <- "list"

# summarize results
lapply(results_list, summarize, sum(FDR<0.05)) %>%
	bind_rows(.id="Contrast") %>% knitr::kable()


## ---- save results

# save results as excel document
myfile <- file.path(root,"tables","SWIP-MSstatsTMT-Results.xlsx")
write_excel(results_list,myfile)
message("saved: ", myfile)


## save data in root/data

# save msstats_prot -- the normalized protein data
myfile <- file.path(root, "data", "msstats_prot.rda")
save(msstats_prot, file = myfile, version = 2)
message("saved: ", myfile)

# save msstats_results
myfile <- file.path(root, "data", "msstats_results.rda")
save(msstats_results, file = myfile, version = 2)
message("saved: ", myfile)

# save sig_prots
myfile <- file.path(root, "data", "msstats_sig_prots.rda")
save(sig_prots, file = myfile, version = 2)
message("saved: ", myfile)
