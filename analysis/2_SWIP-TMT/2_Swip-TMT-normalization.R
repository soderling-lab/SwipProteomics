#!/usr/bin/env Rscript

# title: Swip TMT Proteomics
# description: preprocessing and statistical analysis of Swip TMT proteomics
# experiment performed by JC
# author: Tyler W A Bradshaw

## ---- input

# project's root directory
root = "~/projects/SwipProteomics"

# INPUT data is zipped up in root/data/
zip_file = "TMT.zip"
input_meta = "TMT-samples.csv"
input_data = "TMT-raw-peptide.csv"


# ---- functions

## NOTE: The functions used in this script are not robust. They were written
## to work with the input arguments that they are provided, and in many cases
## will not perform as expected if passed different arguments. I attempted to
## keep the data in a tidy-ish format throughout. This decision makes some
## operations like plotting easier, but makes other operations like
## normalization more cumbersome and computationally costly.


## ---- prepare the R workspace

# prepare the R workspace for the analysis

# load required packages and functions
suppressPackageStartupMessages({
	library(dplyr) # for manipulating data
	library(getPPIs) # twesleyb/getPPIs for mouse PPIs
	library(geneLists) # twesleyb/geneLists for gene mapping fun
	library(data.table) # for working with tables
	library(doParallel) # for parallel processing
})

# load project specific functions and data
devtools::load_all(root, quiet=TRUE)

# project directories:
datadir <- file.path(root, "data") # Key pieces of data saved as rda
rdatdir <- file.path(root, "rdata") # Temporary data files
tabsdir <- file.path(root, "tables") # Output tables saved as excel files
downdir <- file.path(root, "downloads") # Misc downloads/temporary files

# create project output directories if necessary
if (!dir.exists(datadir)) { dir.create(datadir) }
if (!dir.exists(downdir)) { dir.create(downdir) }
if (!dir.exists(rdatdir)) { dir.create(rdatdir) }
if (!dir.exists(tabsdir)) { dir.create(tabsdir) }


## ---- load the raw data and sample info

# datadir [1] ~/projects/SwipProteomics/data"

# extract the raw TMT data from zipped file
myfile <- file.path(datadir, zip_file)
unzip(myfile, exdir=downdir) # unzip into root/downloads/

# read TMT.csv data into R with data.table::fread
myfile <- file.path(downdir, tools::file_path_sans_ext(zip_file), input_data)
peptides <- data.table::fread(myfile)

# load sample information
myfile <- file.path(downdir,tools::file_path_sans_ext(zip_file),input_meta)
samples <- data.table::fread(myfile)

# format cfg force column -- this is the Force in g's used to obtain
# the subcellular fraction.
samples$"Cfg Force (xg)" <- formatC(samples$"Cfg Force (xg)",big.mark=",")


## ---- map all Uniprot accession numbers to stable entrez IDs

message("\nCreating gene identifier map.")

# first, remove any non-mouse proteins from the data
peptides <- peptides %>% filter(grepl("OS=Mus musculus",Description))

# remove these Immunoglobin proteins:
ig_prots <- c("P01631","P01646","P01665","P01680","P01746","P01750",
	      "P01786","P01864","P01878","P03975","P06330","P03987")
peptides <- peptides %>% filter(Accession %notin% ig_prots)

# collect all Uniprot IDs
uniprot <- unique(peptides$Accession)

# map Uniprot IDs to Entrez using online MGI batch query function
entrez <- geneLists::queryMGI(uniprot)
names(entrez) <- uniprot

# map any remaining missing IDs by hand
message("Mapping missing IDs by hand.\n")
missing <- entrez[is.na(entrez)]
mapped_by_hand <- c(P05214=22144, P0CG14=214987, P84244=15078)
entrez[names(mapped_by_hand)] <- mapped_by_hand

# check: Have we successfully mapped all Uniprot IDs to Entrez?
check <- sum(is.na(entrez)) == 0
if (!check) { stop("Unable to map all UniprotIDs to Entrez.") }

# map entrez ids to gene symbols using twesleyb/getPPIs
# NOTE: getIDs is just an easier-to-use wrapper around AnnotationDbi::mapIDs
# NOTE: you need the org.Mm.eg.db package from Bioconductor for mapping mouse genes
gene_symbols <- geneLists::getIDs(entrez,from="entrez",to="symbol",species="mouse")

# create gene identifier mapping data.table
gene_map <- data.table(uniprot = names(entrez),
                       entrez = entrez,
	               symbol = gene_symbols)
gene_map$id <- paste(gene_map$symbol,gene_map$uniprot,sep="|")


## ---- tidy-up the input data from Proteome Discover

# convert PD df into tidy data.frame
message("\nLoading raw data from Proteome Discover (PD2.2).")
cols <- colnames(peptides)[!grepl("Abundance",colnames(peptides))]
tidy_peptide <- tidyProt(peptides, id.vars=cols)

# annotate tidy data with additional meta data from samples
tidy_peptide <- left_join(tidy_peptide,samples,by="Sample")

# summary of peptide/protein quantification:
message(paste("\nSummary of initial peptide/protein quantification",
	      "after removing contaiminants:"))
n_samples <- length(unique(tidy_peptide$Sample))
n_proteins <- length(unique(tidy_peptide$Accession))
n_peptides <- length(unique(tidy_peptide$Sequence))
df <- data.frame("Samples"=as.character(n_samples),
		 "Proteins"=formatC(n_proteins,big.mark=","),
		 "Peptides"=formatC(n_peptides,big.mark=","))
knitr::kable(df)


## ---- perform sample loading normalization

# perform sample normalization
# normalization is done for each mixture independently

# NOTE: Grouping by Experiment, Channel does't work because
# Sample TMT channels (e.g. 126N were used in different expirements)


message("\nPerforming sample loading normalization.")

sl_peptide <- normSL(tidy_peptide, groupBy=c("Experiment","Sample"))
save(sl_peptide,file="sl_peptide.rda",version=2)
stop()


## ---- impute peptide-level missingness

# Impute missing peptide values with k-nearest neighbors (KNN) algorithm.
# * Missing QC values will not be imputed.
# * Peptides (rows) with more than 50% missingness will not be imputed.
# Values in these rows are masked (replaced with NA).
# NOTE: KNN imputing is done for each experimental group seperately.

message("\nImputing a small number of missing peptide values.")

imputed_peptide <- imputeKNNpep(sl_peptide, groupBy="Experiment",
				samples_to_ignore="SPQC", quiet=FALSE)


## ---- examine reproducibility of QC measurements

# assess reproducibility of QC measurements and remove QC peptide measurments
# that are irreproducible

# this strategy was adapted from Ping et al., 2018 (pmid: 29533394)
# For each experiment, the ratio of SPQC tech replicates is calculated
# These ratios are then binned based on average Intensity into
# 5 bins. For each bin, ratios that are outside
# +/- 4x standard deviation from the bin's mean (should be centered at 0) are removed

message("\nRemoving peptides with irreproducible QC measurements.")

filt_peptide <- filtQC(imputed_peptide,controls="SPQC", quiet=FALSE)


## ---- summarize to protein level by Sum

message("\nSummarizing proteins as the sum of their peptides.")

proteins <- sumProt(filt_peptide)


## ---- SL normalization

message("\nPerforming sample loading normalization between experiments.")

sl_protein <- normSL(proteins, groupBy="Sample")


## ---- IRS normalization

# equalize QC measurements between experiments. Adjusts protein
# measurements of biological replicates simultaneously

# accounts for protein quantification by different peptides in
# each experiment

message(paste("\nStandardizing protein measurements between",
	      "experiments by IRS normalization."))

irs_protein <- normIRS(sl_protein,controls="SPQC",robust=TRUE)

# check, protein-wise QC measurements are now equal in a random protein
# generate list of QC proteins, and check the average of the three Exps
qc_proteins <- irs_protein %>% filter(Treatment == "SPQC") %>%
	group_by(Accession,Experiment) %>%
	dplyr::summarize(Treatment = unique(Treatment),
		  "Mean(Intensity)" = mean(Intensity,na.rm=TRUE),
		  .groups = "drop") %>% group_by(Accession) %>% group_split()
message("\nIntra-experimental QC means are equal after IRS normalization:")
knitr::kable(sample(qc_proteins,1))


## ---- protein level filtering

# remove proteins that:
# * are identified by a single peptide
# * contain too many (>50%) missing values
# * contain any missing QC values
# * have outlier protein measurements:
#       mean(logRatio(replicates)) outside +/- nSD * mean(binned logRatio))

#FIXME: filtProt is slow!

message(paste("\nFiltering proteins, this may take several minutes."))

filt_protein <- filtProt(irs_protein, controls="SPQC", nbins=5, nSD=4, summary=TRUE)

# at this point there are no remaining missing values
check <- sum(is.na(filt_protein$Intensity)) == 0
if (!check) { stop("Why are there missing values!") }


## ---- combine the final normalized data and sample meta data

cols <- intersect(colnames(samples),colnames(filt_protein))

swip_tmt <- filt_protein %>%
	# combine with sample meta data
	left_join(samples,  by = cols) %>%
	# drop QC
	filter(Treatment != "SPQC") %>%
	# other munge
	mutate(Genotype = Treatment) %>%
	mutate(Protein = Accession) %>%
	mutate(Mixture = gsub("Exp","M", Experiment)) %>%
	mutate(BioFraction = Fraction) %>%
	mutate(Abundance = log2(Intensity)) %>%
	mutate(Condition = interaction(Genotype,BioFraction)) %>%
	mutate(Subject = as.numeric(interaction(Mixture,Genotype))) %>%
	dplyr::select(Protein,Mixture,Genotype,BioFraction,Condition,Subject, Intensity,Abundance) %>%
	# calculate relative_Intensity (sum normalization)
	group_by(Protein) %>%
	mutate(rel_Intensity = Intensity / sum(Intensity))


## ----  save key results

# final normalized protein in tidy format as rda object
myfile <- file.path(datadir,"swip_tmt.rda")
save(swip_tmt,file=myfile,version=2)
message("saved: ", myfile)

# save gene_map
myfile <- file.path(datadir,"swip_gene_map.rda")
save(gene_map,file=myfile,version=2)
message("saved: ", myfile)
