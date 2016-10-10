library(scater)
library(tidyverse)
library(magrittr)



library(RTCGA.clinical)
data("OV.clinical")

## Read in Kallisto quantified data
ov <- read_tsv("../data/TCGA_OV_tpm.tsv.gz")
sample_names <- names(ov)[-1]
ov <- data.frame(ov)
names(ov)[1] <- "feature_id"
rownames(ov) <- ov$feature_id; ov$feature_id <- NULL
names(ov) <- sample_names


## Find feature names (ensembl transcript id + hgnc symbol) and gene types
id_split <- strsplit(rownames(ov), "|", fixed = TRUE)
feature_names <- sapply(id_split, function(x) paste0(x[1], "_", x[6]))
gene_type <- sapply(id_split, `[`, 8)
rownames(ov) <- feature_names

## Match CGHubAnalysisID with comparible ID in OV.clinical

id_map <- read_csv("../data/TCGA_ID_MAP.csv") %>% filter(Disease == "OV")

to_tcga_barcode <- function(x) tolower(paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "-"))
id_map %<>% mutate(patient_barcode = sapply(AliquotBarcode, to_tcga_barcode))

## Get clinical data associated with IDs
ov_clinical <- OV.clinical %>% 
  filter(patient.bcr_patient_barcode %in% id_map$patient_barcode) %>% 
  select(patient_barcode = patient.bcr_patient_barcode,
         patient.days_to_death,
         patient.age_at_initial_pathologic_diagnosis, 
         patient.days_to_last_followup, patient.days_to_tumor_progression,
         patient.days_to_tumor_recurrence,
         patient.stage_event.clinical_stage,
         patient.tumor_stage) %>% 
  mutate(censored = is.na(patient.days_to_death))

ov_clinical %<>% inner_join(select(id_map, patient_barcode, CGHubAnalysisID, AliquotBarcode),
                            by = c("patient_barcode"))

## Tidy up classes
for(covariate_index in grep("days", names(ov_clinical)))
  ov_clinical[[covariate_index]] <- as.numeric(ov_clinical[[covariate_index]])

plate <- factor(sapply(strsplit(ov_clinical$AliquotBarcode, "-"), `[`, 6))
centre <- factor(sapply(strsplit(ov_clinical$AliquotBarcode, "-"), `[`, 7))
ov_clinical$plate <- plate
ov_clinical$centre <- centre

## Put clinical data in order corresponding to expression data
ov_clinical <- ov_clinical[match(colnames(ov), ov_clinical$CGHubAnalysisID), ]
stopifnot(all.equal(colnames(ov), ov_clinical$CGHubAnalysisID))

ov_clinical_pd <- data.frame(ov_clinical)
rownames(ov_clinical_pd) <- ov_clinical_pd$CGHubAnalysisID

sce <- newSCESet(tpmData = as.matrix(ov), 
                 phenoData = AnnotatedDataFrame(ov_clinical_pd))

sce <- sce[, sce$centre != 31] # get rid of the crappy centre
sce$plate <- factor(as.character(sce$plate))

save(sce, file = "../data/sce_ovarian_kallisto.Rdata")

#voom_cpm <- limma::voom(counts(sce))
#exprs(sce) <- voom_cpm$E

## Batch effect
plotPCA(sce, colour_by = "plate", ncomponents = 3)
plotQC(sce, type = "find", var = "plate")

# Biological effect
plotPCA(sce, colour_by = "patient.days_to_birth", ncomponents = 3)
plotPCA(sce, colour_by = "censored", ncomponents = 3)

plotQC(sce, type = "find", var = "patient.days_to_birth")
plotQC(sce, type = "find", var = "censored")


## QC Metrics
is_exprs(sce) <- exprs(sce) > 0
sce <- calculateQCMetrics(sce)

# Exploratory analysis ----------------------------------------------------
sc <- sce[rowMeans(exprs(sce)) > 1 & matrixStats::rowVars(exprs(sce)) > 0, ]
sc$random_two_factor <- sample(c(T,F), ncol(sc), replace = TRUE)

plotQC(sc, type = "find", var = "censored")
plotQC(sc, type = "find", var = "random_two_factor")
plotQC(sc, type = "expl", var = c("censored", "random_two_factor"))


# Regress out plate effect ------------------------------------------------

## first remove low varying genes
qplot(matrixStats::rowVars(exprs(sce)))
sc <- sce[matrixStats::rowVars(exprs(sce)) > 1, ]

batch <- sc$plate
design <- model.matrix(~ 1, data = pData(sc))

exprs_combat <- ComBat(exprs(sc), batch, design)

sc2 <- sc
exprs(sc2) <- exprs_combat

plotPCA(sc2, ncomponents = 3, colour_by = "plate")
plotPCA(sc2, ncomponents = 3, colour_by = "censored")
plotPCA(sc2, ncomponents = 3, colour_by = "patient.days_to_birth")
plotPCA(sc2, ncomponents = 3, colour_by = "patient.days_to_death")
plotPCA(sc2, colour_by = "patient.stage_event.clinical_stage", ncomponents = 3)

plotQC(sc2, type = 'find', var = 'plate')
plotQC(sc2, type = 'find', var = 'censored')
plotQC(sc2, type = 'find', var = 'patient.days_to_birth')
plotQC(sc2, type = 'find', var = 'patient.days_to_death')
plotQC(sc2, type = 'find', var = 'patient.stage_event.clinical_stage')


# Mean variance analysis --------------------------------------------------
library(ggrepel)

hgnc_symbols <- sapply(strsplit(featureNames(sc2), "_"), `[`, 2)

dvar <- data_frame(mean = rowMeans(exprs(sc2)),
                   var = matrixStats::rowVars(exprs(sc2)),
                   hgnc_symbols, gene = featureNames(sc2))

ggplot(dvar, aes(x = mean, y = var)) + geom_point(color = "grey") + xlab("Mean") +
  ylab("Variance") +
  geom_text_repel(data = filter(dvar, var > 8.5), aes(label = hgnc_symbols),
                  color = 'blue')

is_exprs(sce) <- exprs(sce) > 0
sce <- calculateQCMetrics(sce)
plotHighestExprs(sce, col_by_variable = "plate")

plotQC(sc2, type = "expl", variables = c("plate", "censored", "patient.age_at_initial_pathologic_diagnosis"))#,
                                        




# Try some sweet phenotime ------------------------------------------------
library(rstan)
library(coda)
library(MCMCglmm)

X <- model.matrix(~ censored + patient.age_at_initial_pathologic_diagnosis, pData(sc2))
X <- t( scale( X[,-1] ) ) # remove intercept

Y <- t( scale( t( exprs(sc2)[matrixStats::rowVars(exprs(sc2)) > 5, ] ) ) )

X <- rbind(X, rnorm(ncol(X))) # REMOVE

N <- ncol(Y)
# stopifnot(N == nrow(X))

G <- nrow(Y); P = nrow(X)

dlist <- list(N = N, G = G, P = P, X = X, Y = Y)

model <- stan_model("~/oxford/phenot/synthetic/phenot.stan",
                    model_name = "phenotime")

fit <- vb(model, dlist, grad_samples = 3)

tmap <- posterior.mode(mcmc(extract(fit, "pst")$pst))
pca <- prcomp(t(Y))

df <- data_frame(tmap, pc1 = pca$x[,1], pc2 = pca$x[,2], age = sc2$patient.age_at_initial_pathologic_diagnosis,
                 censored = sc2$censored,
                 death = sc2$patient.days_to_death)

ggplot(df, aes(x = pc1, y = tmap, color = age)) + geom_point() + viridis::scale_color_viridis()

cowplot::plot_grid( ggplot(df, aes(x = pc1, y = tmap, color = age)) + 
                    geom_point() + viridis::scale_color_viridis(),
                    ggplot(df, aes(x = pc1, y = tmap, color = censored)) + geom_point() + 
                      scale_color_brewer(palette = "Set1") )

ggplot(df, aes(x = pc1, y = tmap, color = censored)) + geom_point() + 
  scale_color_brewer(palette = "Set1")


ggplot(df, aes(x = pc1, y = tmap, color = age)) + 
  geom_point() + viridis::scale_color_viridis(name = "Age") +
  theme_classic() + xlab("Principal component 1") +
  ylab("MAP covariate latent trajectory")

ggsave("~/Dropbox/Oxford/talks/confirmation_of_status/figs/phenotime1.png",
       width = 6, height = 4)

tidy_beta <- function(beta, n) {
  names(beta) <- rownames(Y)
  beta_tidy <- gather(beta, gene, value)
  beta_tidy$beta <- n
  return(beta_tidy)
}

beta_1_tidy <- tidy_beta( as_data_frame( extract(fit, "beta")$beta[,,1] ), "beta_1" )
beta_2_tidy <- tidy_beta( as_data_frame( extract(fit, "beta")$beta[,,2] ), "beta_2" )
beta_3_tidy <- tidy_beta( as_data_frame( extract(fit, "beta")$beta[,,3] ), "beta_3" )


beta <- bind_rows(beta_1_tidy, beta_2_tidy, beta_3_tidy)

beta$beta <- plyr::mapvalues(beta$beta, from = c("beta_1", "beta_2", "beta_3"),
                             to = c("censored", "patient_age", "random"))


filter(beta, gene %in% rownames(Y)[1:20]) %>% 
  ggplot(aes(x = gene, y = value, color = beta)) + geom_boxplot()

beta_sig <- beta %>% group_by(gene, beta) %>% 
  summarise(lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>% 
  mutate(significant = upper < 0 | lower > 0)


## Plot for chris
sig_genes <- beta_sig %>% filter(significant) %>% extract2("gene")

beta_for_plot <- filter(beta, gene %in% sig_genes)
beta_for_plot$gene <- sapply(strsplit(beta_for_plot$gene, "_"), `[`, 2)

ggplot(beta_for_plot, aes(x = beta, y = value, color = beta)) + 
  geom_boxplot() + facet_wrap(~ gene, scales = "free_y") + 
  theme(legend.position = "none") + geom_hline(yintercept = 0, linetype =2)


mean_expr_df <- data_frame(mean = rowMeans(Y), var = matrixStats::rowVars(Y),
                           gene = rownames(Y), tau = tau_map)

beta_sig <- inner_join(beta_sig, mean_expr_df, by = "gene")

ggplot(beta_sig, aes(x = significant, y = mean)) + geom_boxplot()
ggplot(beta_sig, aes(x = significant, y = var)) + geom_boxplot()
ggplot(beta_sig, aes(x = significant, y = tau)) + geom_boxplot()

View(filter(beta_sig, beta == "censored", significant))


# Some plots for chris ----------------------------------------------------

plotPCA(sce, colour_by = "ENST00000464611.1_ACTB")

plotExpression(sc2, x = "plate", feature = "ENST00000464611.1_ACTB")


mean_var_df <- data_frame(gene = featureNames(sce), 
                          mean = rowMeans(exprs(sce)),
                          var = matrixStats::rowVars(exprs(sce)))

ggplot(mean_var_df, aes(x = mean, y = var)) +
  geom_point(color = "grey") +
  geom_point(data = filter(mean_var_df, grepl("TP53", mean_var_df$gene)),
             color = "blue") +
  ggtitle("blue genes all tp53 transcripts")

# ACTB vs most variable gene

mostvar_genes <- arrange(mean_var_df, desc(var)) %>% 
   extract2("gene")


exprs_df <- data.frame(t(exprs(sce)[mostvar_genes[1:10], ]))
names(exprs_df) <- sapply(strsplit(mostvar_genes[1:10], "_"), `[`, 2)
pairs(exprs_df)
