#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: do something

# INPUT:
root = "~/projects/SwipProteomics"
fig_width = 5
fig_height = 5

## OUTPUT 
# * pdf of module protein pca plot in figs/Modules

## ---- Prepare the workspace 

# Load renv
renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)

# global imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})


# project directories:
figsdir <- file.path(root,"figs","Samples")

# if necessary, create figsdir
if (!dir.exists(figsdir)) {
	dir.create(figsdir)
}

## ---- Prepare the data for ploting 

# load the proteomics data
data(swip_tmt)

# cast into a matrix
dm <- swip_tmt %>%
	reshape2::dcast(Protein ~ interaction(Mixture, Genotype, BioFraction), 
			value.var= "Abundance") %>%
		as.data.table() %>% as.matrix(rownames="Protein")

# there should be no missing or negative values
is_na <- apply(dm,1,function(x) any(is.na(x)))
subdm <- dm[!is_na,]
is_neg <- apply(subdm,1,function(x) any(x<0))


# PCA
pca <- prcomp(t(subdm[!is_neg,]))
pca_summary <- as.data.frame(t(summary(pca)$importance))
idx <- order(pca_summary[["Proportion of Variance"]],decreasing=TRUE)
pca_summary <- pca_summary[idx,]
top2_pc <- head(pca_summary[["Proportion of Variance"]],2)
names(top2_pc) <- head(rownames(pca_summary),2)

# Plot axis labels:
x_label <- paste0(names(top2_pc)[1],
		  " (PVE: ",round(100*top2_pc[1],2)," %)")
y_label <- paste0(names(top2_pc)[2],
		  " (PVE: ",round(100*top2_pc[2],2)," %)")

# Collect data for plotting.
df <- as.data.frame(pca$x[,names(top2_pc)])
colnames(df) <- c("x","y")

# annotate with group info
geno <- sapply(strsplit(rownames(df),"\\."),"[",2) # <- 
frac <- sapply(strsplit(rownames(df),"\\."),"[",3) # <-
df$group <- interaction(geno,frac)


## ---- generate the plot

plot <- ggplot(df, aes(x,y,color=group)) + geom_point(size=4)
plot <- plot + xlab(x_label)
plot <- plot + ylab(y_label)
plot <- plot + theme(axis.title.x = element_text(color = "black")) 
plot <- plot + theme(axis.title.x = element_text(size = 11))
plot <- plot + theme(axis.title.x = element_text(face = "bold"))
plot <- plot + theme(axis.title.y = element_text(color = "black")) 
plot <- plot + theme(axis.title.y = element_text(size = 11))
plot <- plot + theme(axis.title.y = element_text(face = "bold"))
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(fill=NA))
plot <- plot + scale_x_continuous(expand = c(0,0))
plot <- plot + scale_y_continuous(expand = c(0,0))


## ---- Save to file

myfile <- file.path(figsdir,"sample_PCA.pdf")
ggsave(myfile, plot, width=fig_width,height=fig_height)
message("saved: ", myfile)
