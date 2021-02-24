#!/usr/bin/env Rscript

# title: SwipProteomics
# description: module-level analysis with mixed models
# author: Tyler W Bradshaw

## ---- Inputs

# Input data in root/data/
root = "~/projects/SwipProteomics"


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# load the data
data(swip_tmt)
data(swip_gene_map)
data(swip_partition)  # partition

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
	library(doParallel)
})


## ---- Function

# fit mixed-model to log2 relative (scaled to sum) Intensity
fx <-  log2(rel_Intensity) ~ 0 + Condition + (1|Protein)

fitModule <- function(prots, tidy_prot, fx) {
  # build list of input args for lmerTest
  lmer_args <- list()
  lmer_args[["formula"]] <- fx
  lmer_args[["data"]] <- tidy_prot %>% subset(Protein %in% prots)
  # fit the model with some lmer control
  lmer_args[["control"]] <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- do.call(lmerTest::lmer, lmer_args)
  # assess overall contrast and collect results
  LT <- getContrast(fm,"Mutant","Control")
  result <- lmerTestContrast(fm,LT) %>%
	  mutate(Contrast='Mutant-Control') %>% unique() %>%
	  mutate('nProts'=length(prots))
  return(result)
} #EOF


## ---- main

# loop to fit module-level models and assess overall contrast
modules <- split(names(partition),partition)[-1]
names(modules) <- paste0("M",names(modules))

message("k Modules: ", length(modules))


## ---- loop to fit module-level models and assess contrast

# scale Intensity
tidy_prot <- swip_tmt %>%
	group_by(Protein) %>%
	mutate(rel_Intensity=Intensity/sum(Intensity))

# examine fit of wash complex proteins
washc_prots <- gene_map$uniprot[grepl("Washc*", gene_map$symbol)]
prots = washc_prots[washc_prots %in% tidy_prot$Protein]

fm = lmerTest::lmer(fx, tidy_prot %>% subset(Protein %in% prots))
LT = getContrast(fm,"Mutant","Control")

lmerTestContrast(fm,LT) %>%
	mutate(Contrast = 'Mutant-Control') %>%
	mutate(Pvalue = formatC(Pvalue)) %>%
	mutate(nProts = length(prots)) %>%
	unique() %>% knitr::kable()


# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() -1)

# loop to do module-level analysis
results_list <- foreach(module = names(modules)) %dopar% {
  fitModule(modules[[module]], tidy_prot, fx)
} # EOL
names(results_list) <- names(modules)

## collect results
results_df <- bind_rows(results_list, .id="Module") %>%
	mutate(FDR = p.adjust(Pvalue,method="BH")) %>%
	mutate(Padjust = p.adjust(Pvalue,method="bonferroni")) %>%
	arrange(Pvalue)


## ---- save results

# drop singular col
results_df$isSingular <- NULL

# re-arrange column order
results_df <- results_df %>%
dplyr::select(Module, nProts, Contrast, log2FC,
		percentControl, SE, Tstatistic,
		Pvalue, FDR, Padjust, DF, S2)

# annotate candidate sig modules (Bonferroni padjust < 0.05)
results_df <- results_df %>%
	mutate(candidate = Padjust < 0.05) %>%
	arrange(desc(candidate))

# summary
message("n Sig modules: ", sum(results_df$candidate))

# list of results

# data.frame describing network partition
idx <- match(names(partition),gene_map$uniprot)
df <-  data.table(UniProt = names(partition),
	 Entrez = gene_map$entrez[idx],
	 Symbol = gene_map$symbol[idx],
	 Membership = partition)

# results list:
results_list <- list()
results_list[["Partition"]] <- df %>% arrange(Membership)
results_list[["Module Results"]] <- results_df

# save in root/tables
myfile <- file.path(root,"tables","SWIP-TMT-Module-Results.xlsx")
write_excel(results_list, myfile)
message("saved :", myfile)

# save results as rda in root/data
module_results <- results_df
myfile <- file.path(root,"data", "module_results.rda")
save(module_results, file=myfile, version=2)
message("saved :", myfile)

# save sig modules
sig_modules <- module_results$Module[module_results$candidate]
myfile <- file.path(root,"data", "sig_modules.rda")
save(sig_modules, file=myfile, version=2)
message("saved :", myfile)
