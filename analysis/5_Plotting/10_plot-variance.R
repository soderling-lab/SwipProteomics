#!/usr/bin/env Rscript

# title: SwipProteomics
# author: twab
# description: examine variance attributable to major experimental covariates

## ---- prepare the env

root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# output directory for figures
figsdir <- file.path(root,"figs","Variance")
if (!dir.exists(figsdir)) { dir.create(figsdir,recursive=TRUE) }

# set plotting theme
ggtheme()
set_font("Arial",font_path=file.path(root,"fonts"))


## ---- load the data

data(swip)
data(swip_tmt)
data(swip_partition)

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(gridExtra)
  library(data.table)
  library(doParallel)
})


## ---- work

# register parallel backend
doParallel::registerDoParallel(parallel::detectCores() - 1)

# FIXME: add subject
swip_tmt <- swip_tmt %>% mutate(Subject = as.numeric(interaction(Mixture, Genotype)))

# loop through all proteins, fit the model (all effects modeled as random effects)
fx <- log2(Intensity) ~ (1|BioFraction) + (1|Genotype) + (1|Mixture) + (1|Mixture:Subject) + (1|Subject)

# all prots
proteins <- unique(swip_tmt$Protein)

# loop
pve_list <- foreach(prot = proteins) %dopar% {
  lmer_control <- lme4::lmerControl(check.conv.singular="ignore")
  fm <- lme4::lmer(fx, swip_tmt %>% subset(Protein == prot), control = lmer_control)
  vp <- getVariance(fm)
  pve <- vp/sum(vp)
  pve_df <- as.data.table(t(pve)) %>% mutate(Protein = prot) %>% select(-Fixed)
  return(pve_df)
} #EOL


# collect results
prot_pve <- do.call(rbind, pve_list)

# tidy
df <- prot_pve %>%
	reshape2::melt(id.var="Protein", variable.name = "Parameter", value.name = "Variance")

# summary
df %>% group_by(Parameter) %>%
	summarize(Median=median(Variance), .groups="drop") %>%
	knitr::kable()

# set the order
xpos <- df %>% filter(Parameter == "Residual") %>%
	arrange(-desc(Variance)) %>% select(Protein) %>% unlist(use.names=FALSE)

df$xpos <- match(df$Protein,xpos)

# generate the plot
plot <- ggplot(df)
plot <- plot + aes(x=xpos)
plot <- plot + aes(y=Variance)
plot <- plot + aes(group=Parameter)
plot <- plot + aes(stat=Variance)
plot <- plot + aes(fill=Parameter)
plot <- plot + geom_area()
plot <- plot + xlab("Protein")
plot <- plot + ylab("Percentage of Variance")
plot <- plot + theme(axis.line.x=element_line())
plot <- plot + theme(axis.line.y=element_line())
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(axis.text.x = element_text(color="black",size=11))
plot <- plot + theme(axis.text.x = element_text(angle=0,hjust=1,family="Arial"))
plot <- plot + theme(axis.text.y = element_text(color="black",size=11))
plot <- plot + theme(axis.text.y = element_text(angle=0,hjust=1,family="Arial"))
plot <- plot + theme(panel.border = element_rect(colour = "black", fill=NA))
plot <- plot + scale_x_continuous(expand=c(0,0))
plot <- plot + scale_y_continuous(expand=c(0,0))


## ---- save plot

myfile <- file.path(root,"figs","Variance","Variance_partition.pdf")
ggsave(myfile, plot,height=4.5,width=4.5)
