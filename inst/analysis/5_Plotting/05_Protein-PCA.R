#!/usr/bin/env Rscript

# title: SwipProteomics
# description:
# author: twab

## ---- INPUTs 


## ---- OPTIONs 
fig_width = 5
fig_height = 5


## ---- OUTPUTs 
# * pdf of protein pca plot with module colors


## ---- Prepare the workspace 

# Load renv
root <- "~/projects/SwipProteomics"
renv::load(root,quiet=TRUE)

# Global Imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# Local Imports
devtools::load_all(root, quiet=TRUE)

# Project directories:
figsdir <- file.path(root,"figs","Proteins")

if (!dir.exists(figsdir)) { 
  dir.create(figsdir) 
  message("mkdir ", figsdir) 
}


## ---- Prepare the data for ploting 

# load the data
data(swip_tmt)

# load partition
data(swip_partition)

# load module_colors
data(swip_colors)

# all modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# coerce tidy data to a matrix
dm <- swip_tmt %>% 
	filter(Protein %in% names(partition)) %>%
	group_by(Protein) %>%
	mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	reshape2::dcast(Protein ~ Mixture + Condition, value.var= "rel_Intensity") %>%
	as.data.table() %>% as.matrix(rownames="Protein")

# Drop un-clustered proteins
idx <- rownames(dm) %in% modules[["M0"]]
filt_dm <- dm[!idx,]

# insure there are no missing vals
idx <- apply(dm, 1, function(x) any(is.na(x)))
filt_dm <- dm[!idx,]

# do pca
pca <- prcomp(log2(filt_dm))
pca_summary <- as.data.frame(t(summary(pca)$importance))

# get top 2 components
idx <- order(pca_summary[["Proportion of Variance"]],decreasing=TRUE)
pca_summary <- pca_summary[idx,]
top2_pc <- head(pca_summary[["Proportion of Variance"]],2)
names(top2_pc) <- head(rownames(pca_summary),2)

# Plot axis labels:
x_label <- paste0(names(top2_pc)[1],
		  " (PVE: ",round(100*top2_pc[1],2)," %)")
y_label <- paste0(names(top2_pc)[2],
		  " (PVE: ",round(100*top2_pc[2],2)," %)")

# Collect data for plotting
df <- as.data.frame(pca$x[,names(top2_pc)])
colnames(df) <- c("x","y")

# Generate the plot
plot <- ggplot(df)
plot <- plot + aes(x, y)
plot <- plot + geom_point()
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

# Get plot's data, annotate with protein module membership
df <- plot$data
df$module <- paste0("M",partition[rownames(df)])

# Annotate with color assignment
df$color <- module_colors[df$module]

# Add color to plot
plot <- plot + geom_point(data = df, aes(colour=factor(module)))
plot <- plot + scale_colour_manual(values=df$color)

# Drop the legend
plot <- plot + theme(legend.position = "none")

# Save to file
myfile <- file.path(figsdir,"Protein_PCA.pdf")
ggsave(myfile, plot, width=fig_width,height=fig_height)
