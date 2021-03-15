#!/usr/bin/env Rscript

# title: SwipProteomics
# description: plot WASH complex proteins
# author: Tyler W Bradshaw

## ---- Inputs

# input data in root/data/
root = "~/projects/SwipProteomics"


## ---- Prepare the R environment

renv::load(root, quiet=TRUE)

devtools::load_all(root, quiet=TRUE)

# load the data
data(swip_tmt)
data(washc_prots)

# imports
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(data.table)
	library(doParallel)
})

# project dirs
fontdir <- file.path(root, "fonts")
figsdir <- file.path(root, "figs", "Modules")

if (! dir.exists(figsdir)) {
	dir.create(figsdir,recursive=TRUE)
}

# set plotting theme and font
ggtheme()
set_font("Arial", font_path=fontdir)


## ---- function 

plotWASHC <- function(swip_tmt, prots=washc_prots) {

  # color for Control condition
  wt_color = "#47b2a4"
  mut_color = "#b671af"

  # subset
  subdat <- swip_tmt %>% subset(Protein %in% prots)

  # number of proteins in module
  nprots <- length(unique(subdat$Protein))

  # set factor order (levels)
  subdat$Genotype <- factor(subdat$Genotype,levels= c("Control","Mutant"))
  subdat$BioFraction <- factor(subdat$BioFraction,
			 levels=c("F4","F5","F6","F7","F8","F9","F10"))

  # prepare the data
  df <- subdat %>% 
	  group_by(Protein) %>%
	  mutate(rel_Intensity = Intensity/sum(Intensity)) %>%
	  group_by(Protein, Genotype, BioFraction) %>% 
	  summarize(med_Intensity = median(rel_Intensity), .groups="drop") %>% 
	  mutate(scale_Intensity = scale01(log2(med_Intensity)))

  # get module fitted data by fitting linear model to scaled Intensity
  fx <- scale_Intensity ~ 0 + Genotype:BioFraction + (1|Protein)
  fm <- lmerTest::lmer(fx, df)

  # collect coefficients
  fit_df <- data.table("coef" = names(lme4::fixef(fm)),
		       "fit_y" = lme4::fixef(fm)) %>%
    mutate(Genotype = gsub("Genotype","",sapply(strsplit(coef,"\\:"),"[",1))) %>%
    mutate(BioFraction=gsub("BioFraction","",sapply(strsplit(coef,"\\:"),"[",2)))

  # combine module data and fitted values
  df <- left_join(df, fit_df,by=c("Genotype","BioFraction"))

  # again, insure factor order is correct
  df$Genotype <- factor(df$Genotype,levels=c("Control","Mutant"))
  df$BioFraction <- factor(df$BioFraction,
			   levels=c("F4","F5","F6","F7","F8","F9","F10"))
  
  # get marginal r2 for annot plot title
  r2 <- r.squaredGLMM.merMod(fm)[,"R2m"]
  r2_anno <- paste("(",paste(paste(c("R2.Fixef = "),
			     round(r2,3)),collapse=" | "),")")

  # generate the plot
  plot <- ggplot(df)
  plot <- plot + aes(x = BioFraction)
  plot <- plot + aes(y = scale_Intensity)
  plot <- plot + aes(group = interaction(Genotype,Protein))
  plot <- plot + aes(colour = Genotype)
  plot <- plot + aes(shape = Genotype)
  plot <- plot + aes(fill = Genotype)
  plot <- plot + aes(shade = Genotype)
  plot <- plot + geom_line(alpha=0.25)
  plot <- plot + theme(legend.position = "none")
  plot <- plot + ggtitle(paste0("WASH Complex ", "(n = ",nprots,")\n",r2_anno))
  plot <- plot + ylab("Scaled Protein Intensity")
  plot <- plot + scale_y_continuous(breaks=scales::pretty_breaks(n=5))
  plot <- plot + theme(axis.text.x = element_text(color="black", size=11))
  plot <- plot + theme(axis.text.x = element_text(angle = 0, hjust = 1)) 
  plot <- plot + theme(axis.text.x = element_text(family = "Arial"))
  plot <- plot + theme(axis.text.y = element_text(color="black", size=11))
  plot <- plot + theme(axis.text.y = element_text(angle = 0, hjust = 1)) 
  plot <- plot + theme(axis.text.y = element_text(family = "Arial"))
  plot <- plot + theme(panel.background = element_blank())
  plot <- plot + theme(axis.line.x=element_line())
  plot <- plot + theme(axis.line.y=element_line())
  plot <- plot + geom_line(aes(y=fit_y, group=interaction("fit",Genotype)),
			   linetype="dashed",alpha=1,size=0.75)
  plot <- plot + scale_colour_manual(values=c(wt_color,mut_color))

  return(plot)
} #EOF


## ---- generate plot

plot <- plotWASHC(swip_tmt)


## ---- save as pdf

myfile <- file.path(figsdir,"WASHC_profile.pdf")
ggsave(myfile, plot, height=5, width=5)
message("saved: ", myfile)
