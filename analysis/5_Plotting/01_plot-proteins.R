#!/usr/bin/env Rscript

# title: SwipProteomics
# description: plot protein abundance
# authors: Tyler W Bradshaw

## ---- Inputs 


## ---- Output 
# * a single pdf with plots of all proteins


## ---- Set-up the workspace 

# load renv
root <- "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)

# Load additional functions in root/R/
devtools::load_all(root, quiet=TRUE)

# load the data
data(swip)
data(swip_tmt) 
data(swip_colors) # module_colors
data(swip_gene_map) # gene_map
data(swip_partition)  # partition

# we annotate the plots with MSstatsTMT stats
data(msstats_results) 
data(msstats_sig_prots) # sig_prots


# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})


# project dirs:
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Proteins")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# set theme for the plots:
ggtheme()
set_font("Arial", font_path=fontdir)


## ---- Functions 

plotProtein <- function(protein, prot_df, gene_map,
			 sig_prots, legend=FALSE) {
  # a function that generates the plot

  # annotate title in red if protein has sig change in 'Mutant-Control' contrast
  title_colors <- c("darkred"=TRUE,"black"=FALSE)
  colors <- c("#000000","#303030","#5E5E5E", # WT Blacks
  	    "#942192","#B847B4","#DC6AD7") # Swip Purples

  # subset the data
  gene <- gene_map$symbol[match(protein,gene_map$uniprot)]
  title_color <- names(which(title_colors == (protein %in% sig_prots)))
  df <- subset(prot_df,Protein == protein)

  # insure factor levels are set correctly
  df$BioFraction <- factor(df$BioFraction,
  			      levels=c("F4","F5","F6","F7","F8","F9","F10"))
  df$Mixture <- factor(df$Mixture, levels=c("M1","M2","M3"))
  df$Genotype <- factor(df$Genotype, levels=c("Control","Mutant"))

  # collect FDR stats
  stats_df <- df %>% group_by(Genotype,BioFraction) %>% 
  		summarize(`Max Abundance` = max(Abundance),
  			  FDR = unique(FDR),.groups="drop")
  stats_df$ypos <- 1.02 * max(stats_df$`Max Abundance`)
  stats_df <- stats_df %>% filter(Genotype == "Control")
  stats_df$symbol <- ""
  stats_df$symbol[stats_df$FDR<0.1] <- "."
  stats_df$symbol[stats_df$FDR<0.05] <- "*"
  stats_df$symbol[stats_df$FDR<0.005] <- "**"
  stats_df$symbol[stats_df$FDR<0.0005] <- "***"

  # generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction, y = Abundance)
  plot <- plot + aes(group = interaction(Mixture,Genotype))
  plot <- plot + aes(colour = interaction(Mixture,Genotype))
  plot <- plot + aes(shape=Mixture)
  plot <- plot + aes(group = interaction(Mixture,Genotype))
  plot <- plot + aes(colour = interaction(Mixture,Genotype))
  plot <- plot + geom_point(size=2)
  plot <- plot + geom_line()
  plot <- plot + ggtitle(paste(gene,"|",protein))
  plot <- plot + theme(plot.title=element_text(color=title_color))
  plot <- plot + ylab("log2(Protein Intensity)")

  # annotate with significance stars
  check <- all(is.na(stats_df$FDR))
  any_sig <- any(stats_df$FDR<0.1)
  if (!check & any_sig) { 
    plot <- plot + annotate("text", 
			    x=stats_df$BioFraction, 
			    y=max(stats_df$ypos), 
			    label=stats_df$symbol,size=7)
  }

  # add custom colors and modify legend title and labels
  mylabs <- paste(c(rep('Control',3),rep('Mutant',3)),c(1,2,3))
  plot <- plot + scale_colour_manual(name="Subject", values=colors,labels=mylabs) 
  plot <- plot + ggtitle(paste(gene,protein,sep=" | "))
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black",
						  size=11, angle = 0, 
						  hjust = 1, family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black",
						  size=11, angle = 0, 
						  hjust = 1, family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())

  # remove legend
  if (!legend) {
	  plot <- plot + theme(legend.position = "none")
  }
  return(plot)
} #EOF


## ---- main 

## combine protein data and statistical results

# we will use stats to annotate plots with stars
results <- msstats_results %>% filter(Contrast != 'Mutant-Control') %>% 
	mutate(BioFraction = sapply(strsplit(Contrast,"\\."),"[",3))
shared_cols <- intersect(colnames(swip_tmt),colnames(results))

prot_df <- left_join(swip_tmt, results, by=shared_cols)

# annotate with module membership
prot_dt <- prot_df %>% filter(Protein %in% names(partition))
prot_df$Module <- paste0("M",partition[prot_df$Protein])

# sort proteins by module membership, drop M0
sorted_prots <- as.character(unlist(split(names(partition),partition)[-1]))

stopifnot(all(sorted_prots %in% prot_df$Protein))

message("\nGenerating plots for ", 
	formatC(length(sorted_prots),big.mark=","), " proteins.")


## ---- loop to generate plots

plot_list <- list()

pbar <- txtProgressBar(max=length(sorted_prots),style=3)

for (protein in sorted_prots) {

	# generate a proteins plot
        plot <- plotProtein(protein, prot_df, gene_map, sig_prots)

	# annotate with module assignment
	plot_label <- paste("Module:", partition[protein])
	yrange <- plot$data %>% dplyr::filter(Protein == protein) %>% 
		select(Abundance) %>% range()
	ypos <- yrange[1] - 0.1* diff(yrange)
	plot <-  plot + annotate(geom="label", x=7, y=ypos, label=plot_label)

	# store in list
	plot_list[[protein]] <- plot

	# update pbar
	setTxtProgressBar(pbar, value=match(protein, sorted_prots))
} #EOL
close(pbar)


# Generate a plot with a legend
plot <- plotProtein(swip, prot_df, gene_map, sig_prots, legend = TRUE)
plot_legend <- cowplot::get_legend(plot) 


## ---- save as pdf

# legend
myfile <- file.path(figsdir,"SWIP-Protein-Abundance-legend.pdf")
ggsave(plot_legend, file=myfile, width=4.5, height=4.5)
message("saved: ", myfile)

# plot list
message("\nSaving plots as a single pdf, this will take several minutes.")
myfile <- file.path(figsdir, "SWIP-Protein-Abundance.pdf")
ggsavePDF(plot_list, myfile)
message("saved: ", myfile)
