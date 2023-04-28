# ==============================================================================
# Author(s) : Heini M Natri, heini.natri@gmail.com
# Date: 12/18/2020
# Description: mashr analysis for the Indonesian and European eQTLs
# ==============================================================================

# ======================================
# Environment parameters
# ======================================

# Working directory
setwd("/scratch/hnatri/Indonesian/mashr/")

# ======================================
# Load libraries
# ======================================

library(dplyr)
library(plyr)
library(dmetar)
library(tidyverse)
library(mashr)
library(data.table)
library(dotwhisker)

# ======================================
# Importing eQTL summary statistics
# ======================================

indo_nom <- read.table("/scratch/hnatri/Indonesian/eQTL_methylQTL_result/Indonesian_eQTL_lifted_ALL_CHRS_nominal1.txt")
colnames(indo_nom) <- c("gene", "chr", "start", "end", "strand", "length", "total_tested_vars", "var_id", "var_chr", "var_start", "var_end", "pval", "beta", "is_top_snp")
indo_maf <- "/home/hnatri/Indonesia/MAF_varID_removemissing.tsv"
indo_varinfo <- "/home/hnatri/Indonesia/info_cols.tsv"
indo_nom$gene_rsID <- paste0(indo_nom$gene, "_", indo_nom$var_id)
indo_nom <- indo_nom[!(is.na(indo_nom$var_id) | indo_nom$var_id=="."),]

# Calculating SE
se <- se.from.p(indo_nom$beta, indo_nom$pval, 115, effect.size.type = 'difference',
                calculate.g = FALSE)
indo_nom <- cbind(indo_nom, se)
colnames(indo_nom) <- gsub("StandardError", "se", colnames(indo_nom))

# MAF info
varInfo <- read.table(indo_maf, header=F)
colnames(varInfo) <- c("chr", "pos", "maf", "variant_ID")

infoCols <- read.table(indo_varinfo, header=T)
infoCols$variant_ID <- paste(infoCols$CHROM, infoCols$POS, sep="_")

varInfoCols <- merge(varInfo, infoCols, by="variant_ID")

# From ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv
#eur_eqtls <- c("GTEx_ge_blood", "Lepik_2017_ge_blood", "TwinsUK_ge_blood.all")
gtex_nom <- read.table(gzfile("/scratch/hnatri/Indonesian/eQTL_coloc/GTEx_ge_blood.all.tsv.gz"), header=T)
lepik_nom <- read.table(gzfile("/scratch/hnatri/Indonesian/eQTL_coloc/Lepik_2017_ge_blood.all.tsv.gz"), header=T)
twinsuk_nom <- read.table(gzfile("/scratch/hnatri/Indonesian/eQTL_coloc/TwinsUK_ge_blood.all.tsv.gz"), header=T)

gtex_nom$gene_rsID <- paste0(gtex_nom$molecular_trait_id, "_", gtex_nom$rsid)
gtex_nom <- gtex_nom[!(is.na(gtex_nom$rsid) | gtex_nom$rsid=="."),]
lepik_nom$gene_rsID <- paste0(lepik_nom$molecular_trait_id, "_", lepik_nom$rsid)
lepik_nom <- lepik_nom[!(is.na(lepik_nom$rsid) | lepik_nom$rsid=="."),]
twinsuk_nom$gene_rsID <- paste0(twinsuk_nom$molecular_trait_id, "_", twinsuk_nom$rsid)
twinsuk_nom <- twinsuk_nom[!(is.na(twinsuk_nom$rsid) | twinsuk_nom$rsid=="."),]

# Inputs for mashr
beta_se_p_input_list <- list("indo"=indo_nom, "gtex"=gtex_nom, "lepik"=lepik_nom, "twinsuk"=twinsuk_nom)
beta_input_list <- lapply(beta_se_p_input_list, function(x) x[,(names(x) %in% c("gene_rsID", "beta"))])
beta_input_df <- purrr::reduce(beta_input_list, full_join, by = "gene_rsID")

rownames(beta_input_df) <- beta_input_df$gene_rsID
beta_input_df <- as.matrix(beta_input_df[,!colnames(beta_input_df) %in% c("gene_rsID")])
colnames(beta_input_df) <- names(beta_se_p_input_list)

se_input_list <- lapply(beta_se_p_input_list, function(x) x[,(names(x) %in% c("gene_rsID", "se"))])
sapply(se_input_list, head, n = 6)
se_input_df <- purrr::reduce(se_input_list, full_join, by = "gene_rsID")

rownames(se_input_df) <- se_input_df$gene_rsID
se_input_df <- as.matrix(se_input_df[,!colnames(se_input_df) %in% c("gene_rsID")])
colnames(se_input_df) <- names(beta_se_p_input_list)

# Checking for NAs and NaNs
#beta_input_df[rowSums(is.na(beta_input_df)) > 0, ]
#se_input_df[rowSums(is.na(se_input_df)) > 0, ]
#beta_input_df[rowSums(is.nan(beta_input_df)) > 0, ]
#se_input_df[rowSums(is.nan(se_input_df)) > 0, ]

# Removing rows with NAs
beta_matrix <- na.omit(beta_input_df)
dim(beta_matrix)
se_matrix <- na.omit(se_input_df)
se_matrix <- abs(se_matrix)

# Writing to files
write.table(beta_matrix, "/scratch/hnatri/Indonesian/mashr/beta_matrix.tsv", sep="\t", quote=F)
write.table(se_matrix, "/scratch/hnatri/Indonesian/mashr/se_matrix.tsv", sep="\t", quote=F)
                        
# Removing rows with 0 values (zero in both beta and SE for the same individials causes problems)
beta_matrix_filter_rows <- as.data.frame(beta_matrix) %>% filter_all(any_vars(. %in% c(0)))
beta_matrix_filter <- beta_matrix[!(rownames(beta_matrix) %in% c(rownames(beta_matrix_filter_rows))),]
# Selecting a random subset for testing, half of all tests?
# 100000
beta_matrix_filter <- beta_matrix_filter[sample(nrow(beta_matrix_filter), 11666162), ]
se_matrix_filter <- se_matrix[rownames(beta_matrix_filter),]
# Sanity check
identical(rownames(beta_matrix_filter), rownames(se_matrix_filter))

# Genes shared between all cell populations (included in testing)
genes_included <- unique(sapply(strsplit(rownames(beta_matrix_filter),"_"), `[`, 1))

# ======================================
# Coordinate liftover
# ======================================

# Lifting the eQTL catalog coordinated from GRCh38 to GRCh37
#coordinates <- GRanges(seqnames = gtf_gene$V1, strand=gtf_gene$V7, ranges=IRanges(start=gtf_gene$V4, end=gtf_gene$V5), info=gtf_gene$V9)

# Lifting the start and end coordinates from GRCh38 to hg19
#path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
#ch = import.chain(path)
#ch

#seqlevelsStyle(coordinates) = "UCSC"  # necessary
#coordinates19 = liftOver(coordinates, ch)
#class(coordinates19)
#
#ranges(coordinates19)
#mcols(coordinates19)
#
#coordinates19 <- unlist(coordinates19)
#
#coordinates19_df <- data.frame(seqnames=seqnames(coordinates19),
#                               starts=start(coordinates19),
#                               ends=end(coordinates19),
#                               info=mcols(coordinates19)$info)
#
## Splitting the info column to get the gene id
#coordinates19_df <- separate(data = coordinates19_df, col = info, into = c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype"), sep = "\\;")
#coordinates19_df <- separate(data = coordinates19_df, col = gene_id, into = c("dummy", "gene_id_ensemble"), sep = " ")
#coordinates19_df <- coordinates19_df[,c("seqnames", "starts", "ends", "gene_id_ensemble")]
#colnames(coordinates19_df) <- gsub("gene_id_ensemble", "Geneid", colnames(coordinates19_df))

# ======================================
# Running mashr
# ======================================

mashr_data = mash_set_data(beta_matrix_filter, se_matrix_filter)

# TODO: run with only the canonical covariance matrix
U.c = cov_canonical(mashr_data)
print(names(U.c))

# Running mashr
# The Crucial Rule: fitting the mash model must be performed with either all the
# tests you performed, or – if that causes computational challenges – a large 
# random subset of tests.
# Start the clock!
ptm <- proc.time()
m = mash(mashr_data, U.c, algorithm = "R")
# Stop the clock
proc.time() - ptm
print(get_loglik(m), digits = 10)

# Save an object to a file
saveRDS(m, file = "/scratch/hnatri/Indonesian/mashr/mashr_canoncovs.rds")
                        
# Restore the object
m <- readRDS(file = "/scratch/hnatri/Indonesian/mashr/mashr_canoncovs.rds")

# Extracting posterior summaries
dim(get_lfsr(m))
head(get_lfsr(m))
head(get_pm(m))
head(get_psd(m))

# Significant results
head(get_significant_results(m))
print(length(get_significant_results(m)))
write.table(get_significant_results(m), "/scratch/hnatri/Indonesian/mashr/significant_canoncov.tsv", sep="\t", quote=F)

# Significant for a condition
print(head(get_significant_results(m, conditions=1)))

# Pairwise sharing
print(get_pairwise_sharing(m))
write.table(get_pairwise_sharing(m), "/scratch/hnatri/Indonesian/mashr/pairwise_canoncov.tsv", sep="\t", quote=F)

# Sharing the same sign?
print(get_pairwise_sharing(m, factor=0))
write.table(get_pairwise_sharing(m, factor=0), "/scratch/hnatri/Indonesian/mashr/pairwise_shared_sign_canoncov.tsv", sep="\t", quote=F)

# Estimated mixture proportions
print(get_estimated_pi(m))
barplot(get_estimated_pi(m), las = 2)

# Metaplot on the most significant result
mash_plot_meta(m, get_significant_results(m)[1])

w=5
h=5
pdf(sprintf("/scratch/hnatri/Indonesian/mashr/mashr_metaplot_canoncov.pdf"), width=w, height=h)
mash_plot_meta(m, get_significant_results(m)[1])
dev.off()

# Setting up the data driven covariance matrix

# Step 1: finding strong signals
m.1by1 = mash_1by1(mashr_data)
# 0.05
strong = get_significant_results(m.1by1, 0.05)
length(strong)

# Step 2. obtaining the initial data-driven covariance matrices
# TODO: more PCs?
# https://stephenslab.github.io/mashr/articles/intro_mash_dd.html
U.pca = cov_pca(mashr_data, 2, subset=strong)
print(names(U.pca))

# Step 3: apply extreme deconvolution
U.ed = cov_ed(mashr_data, U.pca, subset=strong)

# Running mashr
m.ed = mash(mashr_data, U.ed, algorithm = "R")
print(get_loglik(m.ed), digits = 10)

# Running mashr by combining canonical and data driven covariaces
m = mash(mashr_data, c(U.c, U.ed), algorithm = "R")
print(get_loglik(m), digits = 10)

# Save an object to a file
saveRDS(m, file = "/scratch/hnatri/Indonesian/mashr/mashr_canoncovs_datacovs.rds")
# Restore the object
m <- readRDS(file = "/scratch/hnatri/Indonesian/mashr/mashr_canoncovs_datacovs.rds")

# Extracting posterior summaries
head(get_lfsr(m))
head(get_pm(m))
head(get_psd(m))

# Significant results
head(get_significant_results(m))
print(length(get_significant_results(m)))
write.table(get_significant_results(m), "/scratch/hnatri/Indonesian/mashr/significant_canoncov_datacov.tsv", sep="\t", quote=F)

# Significant for a condition
print(head(get_significant_results(m, conditions=1)))
length(get_significant_results(m, conditions=1))
indo_specific <- get_significant_results(m, conditions=1)
writeLines(names(indo_specific), "/scratch/hnatri/Indonesian/mashr/indo_specific.tsv")

# Pairwise sharing
print(get_pairwise_sharing(m))
write.table(get_pairwise_sharing(m), "/scratch/hnatri/Indonesian/mashr/pairwise_canoncov_datacov.tsv", sep="\t", quote=F)

# Sharing the same sign?
print(get_pairwise_sharing(m, factor=0))
write.table(get_pairwise_sharing(m, factor=0), "/scratch/hnatri/Indonesian/mashr/pairwise_shared_sign_canoncov_datacov.tsv", sep="\t", quote=F)

# Estimated mixture proportions
print(get_estimated_pi(m))
barplot(get_estimated_pi(m), las = 2)

# Metaplot on the most significant result
mash_plot_meta(m, get_significant_results(m)[1])

w=5
h=5
pdf(sprintf("/scratch/hnatri/Indonesian/mashr/mashr_metaplot_canoncov_datacov.pdf"), width=w, height=h)
mash_plot_meta(m, get_significant_results(m)[1])
dev.off()

