#!/usr/bin/env Rscript

# title: SwipProteomics
# author: Tyler W Bradshaw
# description: plot protein abundance

# plot median Abundance

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"


## ---- Prepare environment 

renv::load(root,quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)


# load the data
data(swip)
data(swip_tmt) 
data(swip_gene_map) # gene_map
data(swip_partition) # partition

# we annotate the plots with stats from MSstatsTMT
data(msstats_results) 
data(msstats_sig_prots) # sig_prots


# other imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs","Proteins")

if (! dir.exists(figsdir)) {
	dir.create(figsdir, recursive=TRUE)
}

# set plotting theme
ggtheme()
set_font("Arial",font_path=fontdir)


## ---- prepare the data

# all protein modules
modules <- split(names(partition),partition)
names(modules) <- paste0("M",names(modules))

# drop 'Mutant-Control' contrasts from results and then
# combine msstats_results with msstats_prot
# we use statistical results from intra-BioFraction comparisons to annotate
# the plots
filt_results <- msstats_results %>% filter(Contrast != "Mutant-Control") %>%
	mutate(BioFraction = sapply(strsplit(Contrast, "\\."),"[",3))
shared_cols <- intersect(colnames(swip_tmt),colnames(filt_results))
prot_df <- left_join(swip_tmt,filt_results,by=shared_cols)

# annotate with module membership
prot_df <- prot_df %>% filter(Protein %in% names(partition)) %>%
	mutate(Module = paste0("M",partition[Protein])) 


## ---- Function 

# a function to generate protein profile plot:

protein = sample(unique(swip_tmt$Protein),1)

plotProfile <- function(protein, gene_map, swip_tmt, msstats_results) {

  # colors for plot
  wt_color <- "#47b2a4"
  mut_color <- "#b671af"

  # map uniprot to gene name
  gene <- gene_map$symbol[match(protein, gene_map$uniprot)]

  # proteins with overall sig change
  sig_prots <- msstats_results %>% 
	filter(Contrast == 'Mutant-Control' & FDR<0.05) %>%  ungroup() %>%
	select(Protein) %>% unlist() %>% as.character() %>% unique()

  # title = darkred if significant difference for overall comparison 
  title_colors <- c("darkred"=TRUE,"black"=FALSE)
  title_color <- names(which(title_colors==protein %in% sig_prots))

  # prepare stats for labeling
  stats_df <- prot_df %>% filter(Protein == protein) %>% 
  	select(Protein,BioFraction,FDR) %>% unique()
  stats_df$stars <- ""
  stats_df$stars[stats_df$FDR < 0.05] <- "*"
  stats_df$stars[stats_df$FDR < 0.005] <- "**"
  stats_df$stars[stats_df$FDR < 0.0005] <- "***"

  # prepare the data for plotting
  # * Abundance = log2(Intensity)
  # * calc scaled intensity
  # * summarize median of 3x mix
  df <- prot_df %>% 
	  filter(Protein == protein) %>% 
	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	  mutate(scale_Intensity = scale01(log2(rel_Intensity))) %>%
	  group_by(Protein,Condition) %>% 
	  summarize(med_Intensity = median(scale_Intensity),
		    ymin = min(scale_Intensity),
		    ymax = max(scale_Intensity), .groups="drop")

    # munge to annotate with Genotype and BioFraction
    condition <- as.character(df$Condition)
    df$Genotype <- factor(sapply(strsplit(condition,"\\."),"[",1),
		        levels=c("Control","Mutant"))
    df$BioFraction <- factor(sapply(strsplit(condition,"\\."),"[",2),
		        levels=c("F4","F5","F6","F7","F8","F9","F10"))

    # generate the plot
    plot <- ggplot(df)
    plot <- plot + aes(x = BioFraction)
    plot <- plot + aes(y = med_Intensity)
    plot <- plot + aes(group = Genotype)
    plot <- plot + aes(colour = Genotype)
    plot <- plot + aes(shape = Genotype)
    plot <- plot + aes(fill = Genotype)
    plot <- plot + aes(shade = Genotype)
    plot <- plot + aes(ymin=ymin)
    plot <- plot + aes(ymax=ymax)
    plot <- plot + geom_line()
    plot <- plot + geom_ribbon(alpha=0.1, linetype="blank")
    plot <- plot + geom_point(size=2)
    plot <- plot + ggtitle(paste(gene," | ",protein))
    plot <- plot + ylab("Scaled Intensity")
    plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
    plot <- plot + theme(axis.text.x = element_text(color="black", size=11))
    plot <- plot + theme(axis.text.x = element_text(angle = 0, hjust = 1)) 
    plot <- plot + theme(axis.text.x = element_text(family = "Arial"))
    plot <- plot + theme(axis.text.y = element_text(color="black", size=11))
    plot <- plot + theme(axis.text.y = element_text(angle = 0, hjust = 1)) 
    plot <- plot + theme(axis.text.y = element_text(family = "Arial"))
    plot <- plot + theme(panel.background = element_blank())
    plot <- plot + theme(axis.line.x=element_line())
    plot <- plot + theme(plot.title = element_text(color=title_color))
    plot <- plot + theme(axis.line.y=element_line())
    plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))
    plot <- plot + scale_fill_manual(values=c(wt_color,mut_color))
    plot <- plot + theme(legend.position = "none")

    # annotate with module membership and significance stars
    yrange <- setNames(range(df$med_Intensity),nm=c("min","max"))
    ypos1 <- yrange["min"]  + 0.05 * diff(yrange) # Module#
    ypos2 <- yrange["max"] + 0.05 * diff(yrange) # stars
    module_label <- paste0("M",partition[protein])
    plot <- plot + annotate(geom="label",x=7.25,y=ypos1,label=module_label)
    plot <- plot + annotate("text",x=stats_df$BioFraction,y=ypos2,
		  label=stats_df$stars,size=7)

  return(plot)
} #EOF


## ---- generate plots 

# protein names sorted by module membership
sorted_prots <- as.character(unlist(split(names(partition),partition)[-1]))

message("\nGenerating plots for ",
	formatC(length(sorted_prots),big.mark=",")," proteins.")

# loop for all proteins
plots <- list()
pbar <- txtProgressBar(max=length(sorted_prots),style=3)
for (protein in sorted_prots) {
	plots[[protein]] <- plotProfile(protein, gene_map, swip_tmt, msstats_results)
        setTxtProgressBar(pbar,value=match(protein,sorted_prots))
} # EOL for proteins
close(pbar)


## ---- save results 

message("\nSaving plots as a single pdf, this will take several minutes.")
myfile = file.path(figsdir,"SWIP-Protein-Profiles.pdf")
ggsavePDF(plots, file=myfile)
