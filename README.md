# SwipProteomics

This repository contains the raw [TMT](./inst/extdata/TMT.zip) 
and [BioID](./inst/extdata/BioID.zip) proteomics data associated with the
analyses performed by 
[Courtland _et al._, 2021](https://www.biorxiv.org/content/10.1101/2020.08.06.239517v1) [1].

The analysis utilizes functions from
[soderling-lab/tidyProt](https://github.com/soderling-lab/tidyProt) to peform
protein- and module- level comparisions (see example below).  Key data and
results are saved as R objects in `data/`. These can be accessed in R using the
`data()` function. For example, load the Swip spatial proteomics partition with
`data(swip_partition, package="SwipProteomics")`.


## Results

The statistical results for protein- and module-level analyses can be found in
the `tables/` directory.

Maps for plasmids used by this study are in the `sequences/` directory.

A cytoscape graph of the Swip spatial proteome partitioned into 49 modules using
[network enhancement](https://github.com/soderling-lab/neten) [2] and the 
[leiden algorithm](https://github.com/soderling-lab/leiden) [3] is available in
[inst/extdata](./inst/extdata/SwipSpatialProteome.cys).

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


## References

__[1]__ Genetic Disruption of WASHC4 Drives Endo-lysosomal Dysfunction and
Cognitive-Movement Impairments in Mice and Humans.  
Courtland J.L., Bradshaw T.W.A., Waitt G., Soderblom E., Ho T., Rajab A.,
Vancini R., Kim I.H., Soderling S.H. (2021). _eLife_.
[preprint](https://www.biorxiv.org/content/10.1101/2020.08.06.239517v1)

__[2]__ From Louvain to Leiden: guaranteeing well-connected communities.   
Traag, V.A., Waltman. L., Van Eck, N.-J. (2018). _Scientific reports_, 9(1), 5233.
[10.1038/s41598-019-41695-z](https://www.nature.com/articles/s41598-019-41695-z)

__[3]__ Network Enhancement as a general method to denoise weighted biological networks.  
Wang B., Pourshafeie A., Zitnik M., Zhu J., Bustamante C.D., Batzoglou S., Leskovec J. (2018).
_Nature Communications_, 9, 3108. [10.1038/s41467-018-05469-x](https://www.nature.com/articles/s41467-018-05469-x)
