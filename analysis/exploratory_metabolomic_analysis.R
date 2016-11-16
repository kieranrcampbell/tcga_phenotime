library(readr)
library(dplyr)
library(scater)
library(biomaRt)
library(gplots)
library(viridis)
library(matrixStats)
library(GO.db)

# markers <- TP53, PTEN, PIK3CA, KRAS, APC, VHL, NF1, NOTCH1

# source("http://bioconductor.org/biocLite.R")
# biocLite(c("RTCGA", "RTCGA.rnaseq"))

library(RTCGA.rnaseq)
library(RTCGA.clinical)
library(RTCGA)
browseVignettes("RTCGA.rnaseq")

data("OV.rnaseq")
data("OV.clinical")


# Metabolic heatmap -------------------------------------------------------

gex <- OV.rnaseq[,-1]
gex <- log2(as.matrix(gex) + 1)
patient_barcode <- OV.rnaseq[,1]

split <- strsplit(colnames(gex), "|", fixed = TRUE)
hgnc_symbols <- sapply(split, `[[`, 1)
entrez_gene_id <- sapply(split, `[[`, 2)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

go <- "GO:0006629" # lipid metabolic process
go <- "GO:0008152" # metabolic process
go <- "GO:0016049" # cell growth
go <- "GO:0044237" # cellular metabolic process

xx <- as.list(GOBPCHILDREN)
ch <- xx[[go]]

lm_genes <- getBM("entrezgene", filters = "go_id",
                   values = ch, mart = ensembl)$entrezgene

is_lipid_gene <- entrez_gene_id %in% lm_genes

gex_lipid <- gex[, is_lipid_gene]

mean_lipid_expr <- colMeans(gex_lipid)
var_lipid_expr <- colVars(gex_lipid)
qplot(mean_lipid_expr, var_lipid_expr) + xlab("Mean") + ylab("Variance")

# lipid_genes_to_use <- var_lipid_expr > 2 & mean_lipid_expr < 7.5
lipid_genes_to_use <- var_lipid_expr > 2 & mean_lipid_expr < 10

cor_mat <- cor(gex_lipid[, lipid_genes_to_use])
diag(cor_mat) <- NA

pst <- prcomp(gex_lipid[, lipid_genes_to_use])$x[,2]

heatmap.2(gex_lipid[, lipid_genes_to_use], trace = "none", col = "viridis")
heatmap.2(cor_mat, trace = "none", col = "viridis")
heatmap.2(gex_lipid[order(pst), lipid_genes_to_use], trace = "none", col = "viridis",
          Rowv = FALSE, Colv = FALSE)


# genes <- c("CIDEA|1149", "PLIN1|5346", "LEP|3952", "CD36|948")
# 
# heatmap.2(gex_lipid[order(pst), colnames(gex_lipid) %in% genes], trace = "none", col = "viridis",
#           Rowv = FALSE, Colv = FALSE)

days_to_death <- as.numeric(dclinical$days_to_death)
# days_to_death[is.na(days_to_death)] <- 5000

qplot(pst, days_to_death) + geom_rug()
cor(pst, days_to_death, use = "na")


# ----

rnaseq_barcodes <- OV.rnaseq$bcr_patient_barcode
clinical_barcodes <- as.character(OV.clinical$patient.bcr_patient_barcode)

barcodes_rnaseq_to_clinical <- function(bcode) {
  bcode <- as.character(bcode)
  bc <- strsplit(bcode, "-")[[1]]
  bc <- paste(bc[1:3], collapse = "-")
  bc <- tolower(bc)
  return(bc)
}

## Convert RNA-seq barcodes to clinical (short) format
rnaseq_bc_short_format <- sapply(rnaseq_barcodes, barcodes_rnaseq_to_clinical)
OV.rnaseq$barcode <- rnaseq_bc_short_format

## Remove duplicates from RNA-seq data
OV.rnaseq <- OV.rnaseq[-which(duplicated(OV.rnaseq$barcode)), ]

## Rearrange clinical to RNA-seq order
rna_in_clinical <- match(OV.rnaseq$barcode, clinical_barcodes)
stopifnot(!any(is.na(rna_in_clinical)))
OV.clinical <- OV.clinical[rna_in_clinical, ]

stopifnot(all(OV.rnaseq$barcode == as.character(OV.clinical$patient.bcr_patient_barcode)))

dclinical <- OV.clinical %>%
  dplyr::select(days_to_death = as.numeric(patient.days_to_death),
                barcode = patient.bcr_patient_barcode,
                age = as.numeric(patient.age_at_initial_pathologic_diagnosis))
rownames(dclinical) <- dclinical$barcode

gex <- OV.rnaseq 
rownames(gex) <- OV.rnaseq$barcode
gex <- dplyr::select(gex, -bcr_patient_barcode, -barcode)

## Construct sceset
sce <- newSCESet(fpkmData = as.matrix(t(gex)),
                 phenoData = new("AnnotatedDataFrame", dclinical))
sce$dies <- !is.na(sce$days_to_death)
sce <- calculateQCMetrics(sce)

## Filter out low expression genes
sce <- sce[rowSums(fpkm(sce)) > 1, ]

## Parse gene names
fn <- featureNames(sce)
fn2 <- sapply(strsplit(fn, "|", fixed = TRUE), `[`, 1)

sce <- sce[-which(duplicated(fn2)), ]
featureNames(sce) <- fn2[-which(duplicated(fn2))]

sce <- plotPCA(sce, colour_by = "pct_dropout", ncomponents = 3, return_SCESet = TRUE)
plotPCA(sce, colour_by = "dies", ncomponents = 3)


plotQC(sce, type = "find", variable = "dies")
plotQC(sce, type = "expl", variable = "dies")

## Focus on patients who die

sced <- filter(sce, dies == TRUE)
sced <- sced[rowSums(fpkm(sced)) > 1, ]
pData(sced)$days_to_death <- as.numeric(sced$days_to_death)
pData(sced)$log_days_to_death <- log2(sced$days_to_death)

sced <- sced[, sced$days_to_death < 3000] # remove outliers 

correlations <- sapply(seq_len(nrow(sced)), function(i) cor(sced$days_to_death, exprs(sced)[i,]))

plot(sced$days_to_death, exprs(sced)[which.max(correlations), ])
plot(sced$days_to_death, exprs(sced)[which.min(correlations), ])

plotPCA(sced, colour_by = "days_to_death", ncomponents = 3)



markers <- c("TP53", "BRCA1", "CSMD3", "NF1", "CDK12", 
             "FAT3", "GABRA6", "BRCA2", "RB1")

sm <- sce[match(markers, featureNames(sce)), ]
stopifnot(all(featureNames(sm) == markers))

col_colors <- plyr::mapvalues(sce$dies, from = c(T,F), to = c("blue", "red"))
heatmap.2(t(redDim(sce)), trace = "none", ColSideColors = col_colors)

plotQC(sm, type = "find", variable = "dies")

plotPCA(sm, colour_by = "log_days_to_death", ncomponents = 3)
plotQC(sm, type = "find", variable = "log_days_to_death")


s <- data.frame(t(exprs(sm)), dtd = sm$log_days_to_death)


xgrep <- function(w, x) x[grep(w, x)]
cnames <- names(BRCA.clinical)
xgrep("death", cnames)
xgrep("tumor", cnames)

View(BRCA.clinical[, grep("death", cnames)])


#data("COAD.rnaseq")

