# SwipProteomics
This repository contains the raw [TMT](./inst/extdata/TMT.zip) 
and [BioID](./inst/extdata/BioID.zip) proteomics data associated with the
analyses performed by 
[Courtland _et al._, 2021](https://www.biorxiv.org/content/10.1101/2020.08.06.239517v1).

The analysis utilizes functions from
[soderling-lab/tidyProt](https://github.com/soderling-lab/tidyProt) to peform
protein- and module- level comparisions (see example below).  Key data and
results are saved as R objects in `data/`. These can be accessed in R using the
`data()` function. For example, load the Swip spatial proteomics partition with
`data(swip_partition, package="SwipProteomics")`.

The statistical results for protein- and module-level analyses can be found in
the `tables/` directory.

Maps for plasmids used by this study are in the `sequences/` directory.

![wash-module](./elife-image.png)

## Usage Example

```R 
# download the repository as an R package
devtools::install_github("soderling-lab/SwipProteomics")

# install tidyProt for statistical functions
devtools::install_github("soderling-lab/tidyProt")

library(dplyr)

library(tidyProt)
library(SwipProteomics)


# load the normalized TMT data
data(swip_tmt)

# washc4's uniprot ID
data(swip)

# fit a model
fx <- log2(Intensity) ~ 0 + Condition + (1|Mixture)
fm <- lmerTest::lmer(fx,data=swip_tmt %>% subset(Protein==swip))

# create a contrast
LT <- getContrast(fm,"Mutant","Control")

# assess contrast 
res <- lmerTestContrast(fm, LT) %>% mutate(Contrast='Mutant-Control') %>% unique()

knitr::kable(res)

```

|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|       S2|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|--:|--------:|:----------|
|Mutant-Control | -1.401866|      0.3784393| 0.0264791|  -52.94235|      0| 28| 0.007362|TRUE       |

lmerTestContrast returns a data.frame with statistics from the model-based
contrast. The column `isSingular=TRUE` in this case indicates that the variance
attributes to `Mixture` is negligible. 


```R

## fit WASH Complex

library(dplyr)
library(SwipProteomics)

data(washc_prots)

# module-level model includes ranef Protein
fx1 <- log2(rel_Intensity) ~ 0 + Condition + (1|Protein)

# fit the model
fm1 <- lmerTest::lmer(fx1, data = swip_tmt %>% subset(Protein %in% washc_prots))

# assess overall 'Mutant-Control' comparison
res1 <- lmerTestContrast(fm1, LT) %>% mutate(Contrast='Mutant-Control') %>% unique()

knitr::kable(res)

```

|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue|  DF|        S2|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|---:|---------:|:----------|
|Mutant-Control | -1.379633|      0.3843165| 0.0392109|  -35.18497|      0| 151| 0.0645747|FALSE      |
