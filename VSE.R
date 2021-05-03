# ==============================================================================
# Author(s) : Heini M Natri, heini.natri@gmail.com
# Date: May 2020
# Description: Variant set enrichment analysis
# ==============================================================================

# ======================================
# Load libraries
# ======================================

# Dependencies require >R/3.6.0
library(VSE)
library(GenomicRanges)

# ======================================
# Environment parameters
# ======================================

# Working directory
setwd("/scratch/hnatri/Indonesian/VSE/")

# ======================================
# Load data
# ======================================

#ld<-loadLd(file.path(system.file("extdata", "ld_BCa_raggr.csv", package="VSE")), type="raggr")
#avs<-makeAVS(ld)
#avs_size <- avsSize(avs)
#head(avs_size)

# Linked SNPs (R2>0.8). Must contain at least five columns: chr, start, end, LD_snp_id, tag_snp_id.
# Pairwise LD values were obtained with PLINK:
# sort eQTL_fdrp001_rsIDs.tsv | uniq > eQTL_fdrp001_rsIDs_unique.tsv;
# plink --vcf ALL_CHRS_WGS_SNParray_phased_imputed_n115_DR2_095_maxmissing03_maf005_maxmaf095_mac2_renamed_APformatheader_GTGPonly_fillGPfield_GP090_fixed_forQTLtools_onlyrsID.vcf --r2 --ld-snp-list eQTL_fdrp001_rsIDs_unique.tsv --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8
qtlvariant_ld <- read.table("eVariants_plink.ld", header=TRUE)
#qtlvariant_ld <- head(qtlvariant_ld, length=100000)

# Note: colnames must include idTag and idLd
colnames(qtlvariant_ld) <- c("tag_snp_chr", "tag_snp_start", "idTag", "LD_snp_chr", "LD_snp_start", "idLd", "R2")
#test <- unique(qtlvariant_ld$idTag)[1200:2400]
#qtlvariant_ld <- qtlvariant_ld[qtlvariant_ld$idTag %in% test ,]
qtlvariant_ld <- qtlvariant_ld[,colnames(qtlvariant_ld) %in% c("LD_snp_chr", "LD_snp_start", "idTag", "idLd")]
#qtlvariant_ld$LD_snp_start <- as.numeric(qtlvariant_ld$LD_snp_start)
#qtlvariant_ld$LD_snp_chr <- as.numeric(qtlvariant_ld$LD_snp_chr)
#qtlvariant_ld$LD_snp_chr <- paste("chr", qtlvariant_ld$LD_snp_chr, sep="")
qtlvariant_ld_GRanges <- makeGRangesFromDataFrame(qtlvariant_ld,
                                                  seqnames.field="LD_snp_chr",
                                                  start.field="LD_snp_start",
                                                  end.field="LD_snp_start",
                                                  ignore.strand=TRUE,
                                                  keep.extra.columns=TRUE)

# Concatenating lists within the GRanges object
#qtlvariant_ld_GRanges <- c(qtlvariant_ld_GRanges[[1]],
#                         qtlvariant_ld_GRanges[[2]],
#                         qtlvariant_ld_GRanges[[3]],
#                         qtlvariant_ld_GRanges[[4]],
#                         qtlvariant_ld_GRanges[[5]],
#                         qtlvariant_ld_GRanges[[6]])
#
# ======================================
# Creating the associated variant set (AVS) 
# ======================================

# Creating the AVS
qtlvariant_AVS <- makeAVS(qtlvariant_ld_GRanges)

qtlvariant_AVS <- readRDS(file = "qtlvariant_AVS.rds")

length(qtlvariant_ld_GRanges)
dim(qtlvariant_AVS)
length(unique(qtlvariant_ld$idTag))

# Saving AVS GRanges
#save(qtlvariant_AVS, file="eVariant_AVS.Rda")
#qtlvariant_AVS <- load("eVariant_AVS.Rda")

# Check the size of each LD block
qtlvariant_AVS_size <- avsSize(qtlvariant_AVS)
head(qtlvariant_AVS_size)

# ======================================
# Constructing Matched Random Variant Sets (MRVSs)
# ============================= =========

eVariant_mrvs.100b <- makeMRVS(qtlvariant_AVS, bgSize=5, mc.cores = 16)
head(eVariant_mrvs.100b)

# Save MRVS for future use, GRanges object
save(eVariant_mrvs.100b, file="eVariant_mrvs.100b.Rda")
load("eVariant_mrvs.100b")

# ======================================
# Loading genomic features
# ======================================

# Loading sample sheet
samples <- read.table("roadmap_encode_samples.tsv", header = TRUE)
head(samples)
samples <- samples[samples$Tissue == "eroadmap", ]

# ======================================
# Enrichment analysis
# ======================================

# Plotting the intersects
qtlvariant_intersect <- intersectMatrix(qtlvariant_AVS, regions = samples,
                                        col = c("white", "grey10"), 
                                        scale = "none", margins = c(5, 5), 
                                        cexRow = 1, cexCol = 0.5, 
                                        Rowv = NA, Colv = NA)

qtlvariant_VSE <- variantSetEnrichment(qtlvariant_AVS, qtlvariant_mrvs.100, samples)

# variantSetEnrichment output is a list of three matrices, in which, the second
# matrix is the normalized null distribution. It is essential that the null
# distribution is a normal distribution. Normality of the null can be confirmed
# using a QQ-plot.

par.original <- par(no.readonly = TRUE)
par(mfrow = c(ceiling(length(samples$Peaks)/3), 3), mai = c(0.1, 0.1, 0.1, 0.1))
VSEqq(qtlvariant_VSE)

w=8
h=8
filename <- "eVariant_nulldist_bg100_encode_epiroadmap.pdf"
pdf(sprintf(filename), width=w, height=h)
print(VSEqq(qtlvariant_VSE))
dev.off()

qtlvariant_VSE_res <- VSESummary(qtlvariant_VSE)
qtlvariant_VSE_res

# Writing results to a file
write.csv(qtlvariant_VSE_res, "eVariant_VSE_res_bg100_encode_epiroadmap.csv")

# Plotting
VSEplot(qtlvariant_VSE, las = 2, pch = 20, cex = 1, cex.main = 0.6, padj = 0.05)

w=12
h=8
filename <- "eVariant_enrichment_bg100_encode_epiroadmap.pdf"
pdf(sprintf(filename), width=w, height=h)
print(VSEplot(qtlvariant_VSE, las = 2, pch = 20, cex = 1, cex.main = 0.6, padj = 0.05))
dev.off()
