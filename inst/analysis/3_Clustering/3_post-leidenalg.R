#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: clean-up results from from running leidenalg executable

## project root dir
root <- "~/projects/SwipProteomics"

## input in root/rdata
input_part <- "swip_partition.csv"

## output is [input_part].rda
output_part <- paste0(tools::file_path_sans_ext(input_part),".rda")


## ---- prepare the R env

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})


## ---- load partition

myfile <- file.path(root,"rdata",input_part)
df <- data.table::fread(myfile,drop=1)

partition_list <- unlist(apply(df,1,function(x) list(x)),
			 recursive=FALSE, use.names=FALSE)

# add 1 bc python is 0-based
part <- partition_list[[1]] + 1

# set small modules to 0
modules <- split(names(part),part)
too_small <- as.numeric(names(which(sapply(modules,length) < 5)))
part[part %in% too_small] <- 0


## status

message("Total number of modules: ", length(unique(part))-1)


## ---- save as rda

# save partition
partition <- part
myfile <- file.path(root,"data", output_part)
save(partition, file = myfile, version = 2)
message("saved: ", myfile)
