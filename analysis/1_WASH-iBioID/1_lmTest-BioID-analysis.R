#!/usr/bin/env Rscript

# title: WASH iBioID Proteomics Analysis
# description: preprocessing and statistical analysis of WASH1 (Washc1) iBioID
# experiment performed by JC
# author: Tyler W A Bradshaw


## ---- inputs

FDR_alpha = 0.05 # FDR significance threshold for protein enrichment
enrichment_threshold = log2(4.0) # enrichment threshold

## Input data in root/data/BioID.zip/
zipfile = "BioID.zip"
datafile = "BioID_raw_protein.csv" 

## Output in root/tables:
# * WASH_BioID_Results.xlsx

## Output in root/data:
# * wash_interactome.rda -- the WASH iBioID proteome (sig enriched prots)

## ---- renv

root <- "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)


## ---- Functions 

check_group_CV <- function(tidy_prot) {
  # Check CV of WASH, Control, and SPQC groups.
  # Calculate the Coefficient of Covarition (CV).
	cv <- function(x) { sd(x)/mean(x) }
	# Annotate with QC and biological replicates groups.
	tmp_dt <- tidy_prot %>% as.data.table()
	tmp_dt$Group <- sapply(strsplit(tmp_dt$Sample," "),"[",2)
	# Calculate protein-wise CV.
	tmp_res <- tmp_dt %>% group_by(Accession, Group) %>% 
		dplyr::summarize(n = length(Intensity), 
				 CV=cv(Intensity),.groups="drop")
	# NOTE: NA can arise if protein was not quantified in all replicates.
	cv_summary <- tmp_res %>% dplyr::filter(!is.na(CV)) %>% 
		group_by(Group) %>% 
		summarize('Mean(CV)' = mean(CV,na.rm=TRUE),
			  'SD(CV)' = sd(CV, na.rm=TRUE),
			  'n' = length(CV),.groups="drop")
	return(cv_summary)
} #EOF


tidyProt <- function(raw_data,id.vars,species=NULL,
		     samples=NULL,summary=FALSE){
	# description: a function to tidy proteomics data.
#	suppressPackageStartupMessages({
#		library(data.table)
#		library(tibble)
#		library(dplyr)
#	})
	dt <- as.data.table(raw_data) %>%
		melt(id.vars = id.vars,
		     variable.name="Sample",
		     value.name="Intensity",
		     variable.factor=FALSE) # Don't coerce to factor.
	# Remove proteins that do not coorespond to species of interest.
	idx <- grepl(paste0("OS=",species),dt$Description)
	if (!all(idx)) {
		n_out <- length(unique(dt$Accession[!idx]))
		msg <- paste(n_out,"proteins are not from", lquote(species),
			     "and will be removed.")
		warning(msg,call.=FALSE)
		dt <- dt %>% filter(grepl(paste0("OS=",species),Description))
	}
	# Insure zeros are NA.
	dt$Intensity[dt$Intensity==0] <- NA
	# Return tidy df.
	return(as.data.table(dt))
} #EOF


normSL <- function(tp, groupBy=NULL){
#	suppressPackageStartupMessages({
#		library(dplyr)
#		library(data.table)
#	})
	# Create an expression that will be evaluated to dynamically 
	# group data based on user provided grouping factors. Default's
	# to Sample -- every replicate.
	tp <- ungroup(as.data.frame(tp))
	cmd <- paste0("group_by(tp, ",paste(groupBy,collapse=", "),")")
	# Calculate column sums for grouped samples.
	# FIXME: if .groups is set to 'drop' to avoid Warning msg, then error?
	data_list <- eval(parse(text=cmd)) %>% 
		summarize(Total=sum(Intensity,na.rm=TRUE),.groups="drop") %>%
		group_split()
	# Calculate normalization factors.
	data_SL <- lapply(data_list,function(x) {
				     x$"Mean(Total)" <- mean(x$Total,na.rm=TRUE)
				     x$NormFactor <- x$"Mean(Total)"/x$Total
				     x$Norm <- x$Total*x$NormFactor
			      return(x) }) %>% bind_rows()
	# Collect normalization factors as named vector.
	SL_factors <- data_SL$NormFactor
	names(SL_factors) <- as.character(data_SL[["Sample"]])
	# Normalize measurements by sample factors.
	tp$Intensity <- tp$Intensity * SL_factors[as.character(tp[["Sample"]])]
	return(as.data.table(tp))
} #EOF


normSP <- function(tp, pool){
	# perform normalization to pooled QC samples
	# store a copy of the data
	tp <- ungroup(tp)
	tp_copy <- tp
	# group pooled together
	tp$Group <- as.numeric(grepl(paste(pool,collapse="|"),tp$Sample))
	tp_list <- tp %>% group_by(Accession,Group) %>% 
		dplyr::summarize(Mean_Intensity=mean(Intensity,na.rm=TRUE),
			  n = length(Intensity), .groups="drop") %>%
	as.data.table() %>% 
	arrange(Accession,Group) %>% 
	group_by(Accession) %>% group_split()
	# loop to calculate normalization factors
        new_list <- list()
	for (i in 1:length(tp_list)){
		x <- tp_list[[i]]
		x$NormFactor <- c(x$Mean_Intensity[2]/x$Mean_Intensity[1],1)
		x$Norm_Mean_Intensity <- x$Mean_Intensity * x$NormFactor
		new_list[[i]] <- x
	}
	tp_list <- new_list
	# collect in a df
	df <- do.call(rbind,tp_list) %>% 
		dplyr::select(Accession,Group,NormFactor)
	# merge with input data
	tp_norm <- left_join(tp,df,by=c("Accession","Group"))
	# Perform normalization step.
	tp_norm$Intensity <- tp_norm$Intensity * tp_norm$NormFactor
	tp_norm <- tp_norm %>% dplyr::select(colnames(tp_copy))
	tp_norm <- as.data.table(tp_norm)
	# return the normalized data
	return(tp_norm)
} #EOF


imputeKNNprot <- function(tidy_prot,ignore="QC",k=10,rowmax=0.5,colmax=0.8,
			  quiet=TRUE){
	# Determine how many missing values each protein has,
	# and then determine the permisible number of missing values 
	# for any given protein.
#	suppressPackageStartupMessages({
#		library(impute)
#		library(dplyr)
#		library(tibble)
#		library(data.table)
#	})
	# Store a copy of the input data.
	tp <- tp_in <- tidy_prot
	# How many missing values are there?
	tp$Intensity[tp$Intensity==0] <- NA
	N_missing <- sum(is.na(tp$Intensity))
	if (N_missing ==0) { 
		message(paste("There are no missing values.",
			   "Returning untransformed data."))
		return(tp_in)
	}
	# Separate data to be imputed.
	if (!is.null(ignore)) {
		# Do not impute:
		tp_ignore <- tp %>% filter(grepl(ignore,Sample)) %>% 
			dplyr::select(Accession,Sample,Intensity) %>% 
			as.data.table
		# Data to be imputed:
		tp_impute <- tp %>% filter(!grepl(ignore,Sample)) %>% 
			as.data.table
	} else {
		tp_ignore <- NULL
	}
	# Cast the data into a matrix.
	dm <- tp_impute %>%
		dcast(Accession ~ Sample, value.var="Intensity") %>% 
		as.matrix(rownames=TRUE)
	# Don't impute rows (proteins) with too many missing values.
	n_missing <- apply(dm,1,function(x) sum(is.na(x)))
	limit <- ncol(dm) * rowmax
	rows_to_ignore <- n_missing > limit
	n_ignore <- sum(rows_to_ignore)
	if (n_ignore > 0) {
		msg <- paste(n_ignore,"proteins have more than",limit,
			     "missing values, these values cannot be imputed,",
			     "and will be ignored.")
		warning(msg,call.=FALSE)
	}
	# Total number of missing values to be imputed.
	n_imputed <- sum(is.na(dm[!rows_to_ignore,]))
	if (!quiet){
		message(paste("There are",n_imputed, "missing",
		      "values that will be replaced by imputing."))
	}
	# Perform KNN imputing.
		# 
	if (quiet) {
		# suppress output from impute.knn
		silence({
			data_knn <- impute.knn(log2(dm[!rows_to_ignore,]),
					       k=k,colmax=colmax,rowmax=rowmax)
		})
	} else {
		data_knn <- impute.knn(log2(dm[!rows_to_ignore,]),
				       k=k,colmax=colmax,rowmax=rowmax)
	}
	# Collect the imputed data
	dm_knn <- dm
	dm_knn[!rows_to_ignore,] <- 2^data_knn$data
	# Melt into tidy df
	dt_knn <- as.data.table(dm_knn,keep.rownames="Accession")
	tp_imputed <- melt(dt_knn,id.vars="Accession",
			  variable.name="Sample",value.name="Intensity")
	# combine with any samples that were ignored
	tp_imputed <- rbind(tp_ignore,tp_imputed)
	# combine with input meta data
	tp_in$Intensity <- NULL
	tp_in$Sample <- as.character(tp_in$Sample)
	tp_imputed$Sample <- as.character(tp_imputed$Sample)
	tp_out <- left_join(tp_in,tp_imputed,by=c("Sample","Accession")) %>%
		as.data.table
	return(tp_out)
} #EOF


## ---- Prepare the workspace 

# imports
suppressPackageStartupMessages({
	library(dplyr) # For manipulating the data
	library(getPPIs) # For mapping gene identifiers
	library(geneLists) # For a list of mito proteins
	library(data.table) # For working with data.tables
})


# project directories:
datadir <- file.path(root,"data") # key datasets
rdatdir <- file.path(root,"rdata") # temp data files
tabsdir <- file.path(root,"tables") # final xlsx tables
downdir <- file.path(root,"downloads") # misc/temp files

# create dirs if they dont exist
if (!dir.exists(datadir)){ dir.create(datadir) }
if (!dir.exists(rdatdir)){ dir.create(rdatdir) }
if (!dir.exists(tabsdir)){ dir.create(tabsdir) }
if (!dir.exists(downdir)){ dir.create(downdir) }


## ---- Load the raw proteomics data 

# extract the raw data from zipped file
myfile <- file.path(datadir, zipfile)
unzip(myfile) # unzip 

# read into R with data.table::fread
myfile <- file.path(getwd(), tools::file_path_sans_ext(zipfile), datafile)
raw_prot <- fread(myfile)

# clean-up
myfile <- file.path(downdir, tools::file_path_sans_ext(zipfile))
unlink(myfile,recursive=TRUE)
unlink("./BioID", recursive=TRUE)

# tidy-up the data
tidy_prot <- tidyProt(raw_prot,species="Mus musculus",
		      id.vars=c("Accession","Description","Peptides"))

## insure that keratins have been removed

# collect vector or keratin proteins
keratins <- tidy_prot %>% 
	filter(grepl("Keratin|keratin",Description)) %>%
	dplyr::select(Accession) %>% 
	unlist() %>% unique()

# remove keratin prots
tidy_prot <- tidy_prot %>% 
	dplyr::filter(Accession %notin% keratins)

## remove mitochondrial contaiminants

# load mitochondrial protein list from twesleyb/geneLists

data(list="mitocarta2", package = "geneLists")

# collect mitocarta entrez ids
mito_entrez <- unlist(mitocarta2, use.names=FALSE)

# rm these additional mito prots
hand_anno_mito <- c("Mtres1", "Gcdh", "Clpx", "Pdk3", "Mrps36","Hscb",
	  "Aldh2","Shmt2","Ciapin1","Ssbp1","Bckdha","Dap3")
hand_anno_entrez <- getPPIs::getIDs(hand_anno_mito, 'symbol', 'entrez', 'mouse')

# to be removed if in mito
mito <- c(mito_entrez,hand_anno_entrez) 

# annotate data with entrez ids
tidy_prot <- tidy_prot %>% 
	mutate(Entrez = getPPIs::getIDs(Accession,'uniprot','entrez','mouse'))

# total number of mito prots to be removed
nMito <- sum(mito %in% tidy_prot$Entrez)

# rm mito prots
tidy_prot <- tidy_prot %>% dplyr::filter(Entrez %notin% mito)

# cast the tidy raw protein data into a matrix -- we will save this as part of
# an excel workbook with the results
tidy_dm <- tidy_prot %>% as.data.table() %>%
	dcast(Accession + Description + Peptides ~ Sample, 
	      value.var = "Intensity")

#message("\nSummary of group CVs:")
#knitr::kable(check_group_CV(tidy_prot))


## ---- create gene map 

# map uniprot to entrez
uniprot <- unique(tidy_prot$Accession)
entrez <- mgi_batch_query(uniprot)

# fix missing
entrez[is.na(entrez)]
missing_entrez <- c("P10853" = 319180)
entrez[names(missing_entrez)] <- missing_entrez

# map entrez to symbols
symbol <- getPPIs::getIDs(entrez, from="Entrez", to="Symbol", species="mouse")

# there should be no missing entrez
stopifnot(!any(is.na(entrez)))

# there should be no missing symbols
stopifnot(!any(is.na(symbol)))

# collect ids as a data.table
gene_map <- data.table(uniprot, entrez, symbol)


## ---- sample loading normalization 

# munge to add Condition, Run, BioReplicate, Subject cols to data
tidy_prot <- tidy_prot %>% 
	mutate(Condition = sapply(strsplit(Sample,"\\ "),"[",2)) %>% 
	mutate(Run = sapply(strsplit(Sample,"\\ "),"[",1)) %>%
	mutate(BioReplicate = sapply(strsplit(Sample,"\\ "),"[",3)) %>%
	mutate(Subject = interaction(Condition, BioReplicate))

# sl normalization
SL_prot <- normSL(tidy_prot, groupBy="Sample")

# Check, column sums (Run total intensity) should now be equal:

#df <- SL_prot %>% group_by(Sample) %>% 
#	summarize("Total Intensity"=sum(Intensity,na.rm=TRUE),.groups="drop")
#knitr::kable(df)


## ---- protein level filtering 

# Remove one hit wonders -- proteins identified by a single peptide
one_hit_wonders <- unique(SL_prot$Accession[SL_prot$Peptides == 1])
n_ohw <- length(one_hit_wonders)

message(paste("... Number of one-hit-wonders:", n_ohw))

# remove proteins with any missing QC vals
dm <- SL_prot %>% 
	reshape2::dcast(Accession ~ Sample, value.var="Intensity") %>% 
	as.data.table() %>% 
	as.matrix(rownames="Accession")

# collect prots with missing QC
idy <- grep("QC",colnames(dm))
idx <- apply(dm[,idy], 1, function(x) any(is.na(x)))

# number of missing:
is_missing <- rownames(dm)[idx]
n_miss <- length(is_missing[is_missing %notin% n_ohw])

# Remove proteins with more than 50% missingness as these cannot be imputed
df <- SL_prot %>% 
	dplyr::filter(!grepl("QC",Sample)) %>% 
	group_by(Accession) %>%
	summarize(N=length(Intensity),
		  n_missing=sum(is.na(Intensity)),.groups="drop")

# number of sparsly quantified proteins
is_sparse <- unique(df$Accession[df$n_missing>0.5*df$N])
n_sparse <- length(is_sparse)

# subset the data
filt_prot <- SL_prot %>% 
	dplyr::filter(Accession %notin% n_ohw) %>% 
	filter(Accession %notin% is_missing) %>%
	filter(Accession %notin% is_sparse)

# final number of quantified proteins
prots <- unique(filt_prot$Accession)
n_prot <- length(prots)


## ---- sample pool normalization 
# perform normalization to QC samples

SPN_prot <- normSP(filt_prot, pool="QC")


## ---- impute missing values 
# Proteins with missing values are less abundant than those without. 
# This is evidence that missing values are MNAR and can be imputed with
# the KNN algorithm.

message("\nImputing missing protein values using the KNN algorithm (k=10).")

imp_prot <- imputeKNNprot(SPN_prot, k=10, rowmax=0.5, colmax=0.8, quiet=FALSE)

## Status:
knitr::kable(check_group_CV(imp_prot))


## ---- statistical testing for SWIP

wash_bioid <- imp_prot

swip <- gene_map$uniprot[match("Washc4",gene_map$symbol)]

# simple linear model:
fx <- log2(Intensity) ~ Condition

# the data:
df <- wash_bioid %>% subset(Accession == swip) %>% filter(Condition != "QC")

# fit the model
fm <- lm(fx, data = df)

# create a contrast:
LT <- coef(fm) # model coefficients
LT[] <- 0
LT['ConditionWASH'] <- 1

#lmTestContrast(fm, LT) %>% mutate(Pvalue=formatC(Pvalue)) %>% knitr::kable()


## ---- loop to perform test for all proteins

# empty lists for models, sigma2, and degrees of freedom for each fit
fm_list <- list()
s2_list <- list()
df_list <- list()

# fit the models
proteins <- unique(wash_bioid$Accession)
for (prot in proteins) {
	df <- wash_bioid %>% 
		subset(Accession == prot) %>% 
		filter(Condition != "QC")
	fm <- lm(fx, data = df)
	fm_list[[prot]] <- fm
	s2_list[[prot]] <- sigma(fm)^2 # calculate sigma2 for moderation
	df_list[[prot]] <- fm$df.residual # get df for moderation
}

# perform t-statistic moderation
eb_var <- limma::squeezeVar(unlist(s2_list), unlist(df_list))
df_prior <- eb_var$df.prior
s2_prior <- eb_var$s2.prior
if (is.null(s2_prior)) { s2_prior <- 0 }

# examine SWIP's result with moderation
df <- wash_bioid %>% subset(Accession == swip) %>% filter(Condition != "QC")
fm <- lm(fx, data=df)

#lmTestContrast(fm, LT, s2_prior, df_prior) %>% mutate(Pvalue=formatC(Pvalue)) %>% 
#	knitr::kable()

# loop to moderate comparisons for all proteins
result_list <- list()

for (prot in proteins) {
	fm <- fm_list[[prot]]
	result_list[[prot]] <- lmTestContrast(fm,LT,s2_prior,df_prior)
}

# collect results
results_df <- dplyr::bind_rows(result_list,.id="Protein")

# calc FDR
results_df <- results_df %>% mutate(FDR = p.adjust(Pvalue, method="BH"))

# anno with sig and up
results_df <- results_df %>% 
	mutate(up = log2FC > enrichment_threshold) %>% 
	mutate(sig = FDR < FDR_alpha) %>%
	mutate(candidate = up & sig)

# annotate with gene Symbols and Entrez IDs
idx <- match(results_df$Protein,gene_map$uniprot)
results_df <- results_df %>% 
	tibble::add_column(Symbol = gene_map$symbol[idx], .after="Protein") %>%
	tibble::add_column(Entrez = gene_map$entrez[idx], .after="Symbol")

# final sort
results_df <- results_df %>% 
	arrange(desc(sig), desc(up), Pvalue, desc(log2FC)) %>%
	dplyr::select(-sig,-up)

#data.table("nSig"=sum(results_df$candidate)) %>% knitr::kable()


## ---- save results

# save statistical results as excel
bioid_results <- results_df
myfile <- file.path(root,"tables","WASH-iBioID-Results.xlsx")
write_excel(list("WASH-BioID"=bioid_results), myfile)
message("saved: ", myfile)

# save final normalized protein as rda
myfile <- file.path(root,"data","wash_bioid.rda")
save(wash_bioid, file=myfile, version=2)
message("saved: ", myfile)

# save bioid results as rda
myfile <- file.path(root,"data","bioid_results.rda")
save(bioid_results, file = myfile, version = 2)
message("saved: ", myfile)

# collect sig enriched prots 
wash_interactome <- bioid_results %>% 
	filter(candidate == TRUE) %>% 
	dplyr::select(Protein) %>% 
	unlist(use.names=FALSE) %>% 
	unique()

# save wash_interactome as rda
myfile <- file.path(root,"data","wash_interactome.rda")
save(wash_interactome, file=myfile, version=2)
message("saved: ", myfile)
