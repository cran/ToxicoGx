## ---- eval = FALSE, message = FALSE, results = 'hide'--------------------
#  devtools::install_github("bhklab/ToxicoGx", ref = "CRAN_submission")

## ---- message = FALSE, fig.width = 8, fig.height = 5---------------------
library(PharmacoGx)
library(ToxicoGx)
# Load the tset 
data(TGGATESsmall)
ToxicoGx::drugGeneResponseCurve(tSets = TGGATESsmall, 
                      duration = c("2", "8", "24"), 
                      cellline = "Hepatocyte", mDataTypes = "rna", 
                      features = "ENSG00000140465_at", 
                      dose = c("Control", "Low", "Middle", "High"),
                      drug = "carbon tetrachloride", 
                      plot.type = "Actual",
                      title = "Effect of Carbon tetra chloride on CYP1A1.",
                      cex = 0.5, cex.main = 1, legend.loc = "topright", 
                      verbose = T)

## ---- results = 'asis'---------------------------------------------------
library(xtable)
data(TGGATESsmall)
# To compute the effect of drug concentration on the molecular profile of the cell
drug.perturbation <- ToxicoGx::drugPerturbationSig(TGGATESsmall,
                                         mDataType = "rna",
                                         cells = "Hepatocyte",
                                         duration = "24",
                                         dose = c("Control", "Low"),
                                         drugs = c("omeprazole", "isoniazid"),
                                         verbose = FALSE)
data(HCC_sig)
res <- apply(drug.perturbation[,,c("tstat", "fdr")],
             2, function(x, HCC){
               return(PharmacoGx::connectivityScore(x = x,
                                        y = HCC[, 2, drop = FALSE],
                                        method = "fgsea", nperm = 100))
             }, HCC = HCC_sig[1:199,])
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <-  cbind(res,"FDR" = p.adjust(res[,2], method = "fdr"))
res <- res[order(res[,3]),]
xtable::xtable(res,
       caption = 'Connectivity Score results for HCC and TG-GATEs PHH gene
       signature.')

