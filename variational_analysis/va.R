library(scater)
library(Rcpp)

source("../cavi/phenotime_cavi.R")
sourceCpp("../cavi/cavi.cpp")

load("../data/sc_coad_tumour_only.Rdata")

scale_vec <- function(x) (x - mean(x)) / sd(x)

# Select genes ------------------------------------------------------------

var_exprs <- matrixStats::rowVars(exprs(sc_tumour))

to_use <- var_exprs > 2

sc <- sc_tumour[to_use, ]


# Construct covariates ----------------------------------------------------

x_metastasis <- 1 * (sc$m_stage != "m0")
x_metastasis[is.na(x_metastasis)] <- mean(x_metastasis, na.rm = TRUE)
x_metastasis <- scale_vec(x_metastasis)

x_msi <- scale_vec( 1 * (sc$msi_status == "msi-h") )

x <- cbind(x_msi, x_metastasis)

y <- scale(t(exprs(sc)))


# Call phenotime ----------------------------------------------------------

pcavi <- phenotime_cavi(y, x, elbo_tol = 0.01)

# Save results
retained_fnames <- featureNames(sc)
save(pcavi, retained_fnames, file = "va_results.Rdata")