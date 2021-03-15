#!/usr/bin/env Rscript

# title: Swip Proteomics Plotting
# description: generate module colors
# authors: Tyler W Bradshaw

## ---- INPUT

# input data in root/data
root = "~/projects/SwipProteomics"

# not clustered == "gray"
NC_color = "#BEBEBE" 

## ---- OUTPUT 
# * color assignments for every module in graph partition


## ---- FUNCTIONS 

str_to_vec <- function(response) {
  # Parse the python dictionary returned as a string from 
  # system(random_color.py)
	vec <- gsub("'","",gsub("\\[|\\]","",
				trimws(unlist(strsplit(response,",")))))
	return(vec)
}


## ---- Prepare the workspace 

# Load renv
renv::load(root, quiet=TRUE)

# Global Imports
suppressPackageStartupMessages({
	library(dplyr)
	library(data.table)
})

# Local Imports 
devtools::load_all(root, quiet=TRUE)

# Load TMT data and partition
data(swip)
data(mut_color)
data(swip_partition)


## ---- Generate Module colors 

# The number of colors we need
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))
n_colors <- length(modules) # NOTE: M0 will be gray. M(Swip) will be purple.

# Path to python script which is a simple script that uses the python 
# port of randomcolors to generate random colors.
script <- file.path(root,"Py","random_color.py")

# Generate n random colors
cmd <- paste(script,"--count", n_colors, "--luminosity", "bright")
response <- system(cmd, intern = TRUE)

#  Parse the response
colors <- toupper(str_to_vec(response))

if (mut_color %in% colors) { stop("Duplicate colors.") }

# Module color assignments
# Initialize a vector for the module colors
module_colors <- rep(NA,length(modules))
names(module_colors) <- names(modules)

# Insure that M0 is gray and WASH community/module is #B86FAD
module_colors["M0"] <- NC_color

wash_module <- names(which(sapply(modules, function(x) swip %in% x)))

module_colors[wash_module] <- mut_color

# The reamining colors are random
idx <- is.na(module_colors)
module_colors[idx] <- sample(colors,sum(idx))


## ---- Save the data 

message(paste("\nSaving colors."))

# Save updated module colors
namen <- paste0(gsub("partition","colors",part_file),".rda")
myfile <- file.path(root,"data", namen)
save(module_colors,file=myfile,version=2)


################################################################################
## save colors for ne_surprise2 partition
################################################################################


## ---- load the data

data(ne_surprise2_partition); subpart <- partition

data(ne_surprise_colors) # module_colors
data(ne_surprise_partition) # partition


## ---- color mapping

# prepare module color vector
part <- setNames(paste0("M",partition), nm = names(partition))
prot_colors <- setNames(module_colors[part], nm = names(part))
color_list <- split(prot_colors, subpart)
new_colors <- sapply(color_list, function(x) (unique(x)))
names(new_colors) <- paste0("M",names(new_colors))
new_colors[["M0"]] <- "#BEBEBE" 
module_colors <- unlist(new_colors)


## ---- save

myfile <- file.path(root,"data","ne_surprise2_colors.rda")
save(module_colors, file=myfile, version=2)
