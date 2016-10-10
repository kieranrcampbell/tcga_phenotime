library(scater)
library(tidyverse)
library(magrittr)
library(matrixStats)
library(ggrepel)

load("../data/sce_ovarian_kallisto.Rdata")

# Mean variance analysis --------------------------------------------------

sc <- sce[matrixStats::rowVars(exprs(sce)) > 1, ]; rm(sce)
is_exprs(sc) <- exprs(sc) > 0
sc <- calculateQCMetrics(sc)

hgnc_symbols <- sapply(strsplit(featureNames(sc), "_"), `[`, 2)

dvar <- data_frame(mean_exprs = rowMeans(exprs(sc)),
                   var_exprs = matrixStats::rowVars(exprs(sc)),
                   hgnc_symbols, gene = featureNames(sc))
fData(sc) <- cbind(fData(sc), dvar)

## group factors to avoid low frequency ones
gt_table <- table(gene_type)
new_gt_names <- names(gt_table)
new_gt_names[gt_table < 50] <- "other"
fData(sc)$short_gene_type <- plyr::mapvalues(fData(sc)$gene_type,
                                             from = names(gt_table),
                                             to = new_gt_names)


ggplot(fData(sc), aes(x = mean_exprs, y = var_exprs, color = short_gene_type)) + 
  geom_point() + xlab("Mean") +
  ylab("Variance") +
  geom_text_repel(data = filter(fData(sc), var_exprs > 8.5), 
                  aes(label = hgnc_symbols))

ggplot(fData(sc), aes(x = mean_exprs, y = var_exprs)) + 
  geom_point(alpha = 0.5) + xlab("Mean") +
  ylab("Variance") +
  geom_point(data = filter(fData(sc), gene_type == "nonsense_mediated_decay"), 
                  aes(label = hgnc_symbols, fill = short_gene_type), 
             color = 'black', shape = 21, size = 3)

ggplot(fData(sc), aes(x = mean_exprs, y = var_exprs)) + 
  geom_point(alpha = 0.5) + xlab("Mean") +
  ylab("Variance") +
  geom_point(data = filter(fData(sc), gene_type == "IG_V_gene"), 
             aes(label = hgnc_symbols, fill = short_gene_type), 
             color = 'black', shape = 21, size = 3)



# we find "ENST00000464611.1_ACTB" as highly expressed and most likely
# a technical effect:

plotExpression(sc, feature = "ENST00000464611.1_ACTB", x = "plate")

## Let's pick out other genes correlated with it that might also show tech effects

actb_exprs <- exprs(sc)["ENST00000464611.1_ACTB", ]
cors_with_actb <- as.numeric(cor(t(exprs(sc)), actb_exprs))
null_cors_with_actb <- as.numeric(cor(t(exprs(sc)), sample(actb_exprs)))

data_frame(actb = cors_with_actb, null = null_cors_with_actb) %>% 
  gather(type, correlation) %>% 
  ggplot(aes(x = correlation, fill = type)) + geom_density(alpha = 0.7) 
  
plotExpression(sc, x = "plate", feature = featureNames(sc)[order(cors_with_actb, decreasing = T)][2])
plotExpression(sc, x = "plate", feature = featureNames(sc)[order(cors_with_actb)][1])


fData(sc)$correlation_with_actb <- as.numeric(cors_with_actb)

ggplot(fData(sc), aes(x = correlation_with_actb, y = var_exprs,
                      color = short_gene_type)) +
  geom_point() + stat_density_2d(data  = filter(fData(sc), short_gene_type %in% c("nonsense_mediated_decay", "IG_V_gene")))

## Let's do dumb NMD and IGV scores

term <- "nonsense_mediated_decay"

get_mean_exprs <- function(term) {
  pathway_genes <- filter(fData(sc), gene_type %in% term) %>% extract2("gene")
  mean_exprs <- colMeans(exprs(sc)[pathway_genes, ])
}

sc$nmd_score <- get_mean_exprs("nonsense_mediated_decay")
sc$igv_score <- get_mean_exprs("IG_V_gene")

plotPhenoData(sc, aes_string(x = "nmd_score", y = "igv_score"))

plotPCA(sc, colour_by = "pct_dropout", ncomponents = 4)
plotPCA(sc, colour_by = "igv_score", ncomponents = 4)
plotPCA(sc, colour_by = "nmd_score", ncomponents = 4)





