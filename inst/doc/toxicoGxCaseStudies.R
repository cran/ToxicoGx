## ----setup, include = FALSE, cache = FALSE, message = FALSE-------------------

library("knitr")

#opts_knit$set(root.dir=normalizePath('../'))

### Chunk options: see http://yihui.name/knitr/options/ ###

## Text results
opts_chunk$set(echo = TRUE, warning = TRUE, message = FALSE, include = TRUE)

## Cache
opts_chunk$set(cache = 3, cache.path = "output/cache/")

## Plots
opts_chunk$set(fig.path = "output/figures/")


## ---- eval = FALSE, message = FALSE, results = 'hide'-------------------------
#  install.packages("ToxicoGx")

## ---- message = FALSE, fig.width = 8, fig.height = 3--------------------------
library(PharmacoGx)
library(ToxicoGx)
library(ggplot2)

# Load the tset 
data(TGGATESsmall)
ToxicoGx::drugGeneResponseCurve(TGGATESsmall, 
                      duration = c("2", "8", "24"), 
                      cell_line = "Hepatocyte", mDataTypes = "rna", 
                      features = "ENSG00000140465_at",
                      dose = c("Control", "Low", "Middle", "High"),
                      drug = "Carbon tetrachloride",
                      ggplot_args = list(labs(title="Effect of Carbon tetra chloride on CYP1A1")),
                      summarize_replicates = FALSE
                      )


## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics('CS1_published.png')

## ---- results = 'asis'--------------------------------------------------------
library(xtable)
#ata("TGGATESsmall")
# To compute the effect of drug concentration on the molecular profile of the cell
drug.perturbation <- ToxicoGx::drugPerturbationSig(tSet = TGGATESsmall,
                                         mDataType = "rna",
                                         cell_lines = "Hepatocyte",
                                         duration = "24",
                                         dose = c("Control", "Low"),
                                         drugs = c("Omeprazole", "Isoniazid"),
                                         returnValues=c("estimate","tstat", "pvalue", "fdr"),
                                         verbose = FALSE)
data(HCC_sig)
res <- apply(drug.perturbation[,,c("tstat", "fdr")],
             2, function(x, HCC){
               return(PharmacoGx::connectivityScore(x = x,
                                        y = HCC[,2,drop = FALSE],
                                        method = "fgsea", nperm = 100))
             },
             HCC = HCC_sig[1:199,])
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <-  cbind(res,"FDR" = p.adjust(res[,2], method = "fdr"))
res <- res[order(res[,3]),]
xtable::xtable(res,
       caption = 'Connectivity Score results for HCC and TG-GATEs PHH gene
       signature')

