library(MatchIt)
library(plyr)
library(ggpubr)
library(grid)
library(gridExtra)

setwd("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Plotting/")

eQTLtoolsResPath <- "/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indonesian_eQTL_lifted_ALL_CHRS_perm10k_FDR001_significant.txt"
methylQTLtoolsResPath <- "/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indonesian_methylQTL_ALL_CHUNKS_perm10k_FDR001_significant.txt"

eQTLannotationsCounts <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indo_eQTL_site_annotation_counts.tsv", header=T)
methylQTLannotationsCounts <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indo_methylQTL_site_annotation_counts.tsv", header=T)

eQTLannotations <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indo_eQTL_site_annotation.tsv", header=T)
methylQTLannotations <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indo_methylQTL_site_annotation.tsv", header=T)
eQTLannotations$var <- paste(eQTLannotations$seqnames, eQTLannotations$start, sep="_")
eQTLannotations$var <- gsub("chr", "", eQTLannotations$var)
methylQTLannotations$var <- paste(methylQTLannotations$seqnames, methylQTLannotations$start, sep="_")
methylQTLannotations$var <- gsub("chr", "", methylQTLannotations$var)

eQTLtoolsRes <- read.table(eQTLtoolsResPath)
colnames(eQTLtoolsRes) <- c("molecular_trait_id", "chromosome", "target_start",
                            "target_end", "strand", "nvars_tested", "var_distance",
                            "var_rsID", "var_chromosome", "var_start", "var_end",
                            "dfs", "dummy1", "1st_param_betadist", "2nd_param_betadist", "nominalp",
                            "slope", "adjp", "adjp_beta", "fdrp", "dunno")
eQTLtoolsRes$var <- paste(eQTLtoolsRes$var_chromosome, eQTLtoolsRes$var_start, sep="_")
#eQTLtoolsRes$var <- paste("chr", eQTLtoolsRes$var, sep="")
#eQTLtoolsResSig <- eQTLtoolsRes[eQTLtoolsRes$fdrp<0.01 , ]
eQTLtoolsResSig <- eQTLtoolsRes
dim(eQTLtoolsResSig)
hist(eQTLtoolsRes$slope)
hist(eQTLtoolsRes$var_distance)
mean(abs(eQTLtoolsRes$slope))

methylQTLtoolsRes <- read.table(methylQTLtoolsResPath)
colnames(methylQTLtoolsRes) <- c("molecular_trait_id", "chromosome", "target_start",
                                 "target_end", "strand", "nvars_tested", "var_distance",
                                 "var_rsID", "var_chromosome", "var_start", "var_end",
                                 "dfs", "dummy1", "1st_param_betadist", "2nd_param_betadist", "nominalp",
                                 "slope", "adjp", "adjp_beta", "fdrp", "dunno")
methylQTLtoolsRes$var <- paste(methylQTLtoolsRes$var_chromosome, methylQTLtoolsRes$var_start, sep="_")
methylQTLtoolsRes$var <- paste("chr", methylQTLtoolsRes$var, sep="")
methylQTLtoolsResSig <- methylQTLtoolsRes[methylQTLtoolsRes$fdrp<0.01 , ]
hist(methylQTLtoolsResSig$slope)

mean(abs(methylQTLtoolsResSig$slope))

# Fractions of QTLs on promoters/enhancers
NeQTLvar <- length(unique(eQTLannotations$var))
NmethylQTLvar <- length(unique(methylQTLannotations$var))
NeQTLpromoter <- length(unique(eQTLannotations[eQTLannotations$annot.type=="hg19_genes_promoters",]$var))
NmethylQTLpromoter <- length(unique(methylQTLannotations[methylQTLannotations$annot.type=="hg19_genes_promoters",]$var))
NeQTLenhancer <- length(unique(eQTLannotations[eQTLannotations$annot.type=="hg19_enhancers_fantom",]$var))
NmethylQTLenhancer <- length(unique(methylQTLannotations[methylQTLannotations$annot.type=="hg19_enhancers_fantom",]$var))

NeQTLpromoter/NeQTLvar
NmethylQTLpromoter/NmethylQTLvar

NeQTLenhancer/NeQTLvar
NmethylQTLenhancer/NmethylQTLvar

# Effect sizes of promoter/enhancer QTLs
eQTLtoolsResPromoter <- eQTLtoolsRes[(eQTLtoolsRes$var %in% eQTLannotations[eQTLannotations$annot.type=="hg19_genes_promoters",]$var) , ]
methylQTLtoolsResPromoter <- methylQTLtoolsRes[(methylQTLtoolsRes$var %in% methylQTLannotations[methylQTLannotations$annot.type=="hg19_genes_promoters",]$var) , ]

eQTLtoolsResEnhancer <- eQTLtoolsRes[(eQTLtoolsRes$var %in% eQTLannotations[eQTLannotations$annot.type=="hg19_enhancers_fantom",]$var) , ]
methylQTLtoolsResEnhancer <- methylQTLtoolsRes[(methylQTLtoolsRes$var %in% methylQTLannotations[methylQTLannotations$annot.type=="hg19_enhancers_fantom",]$var) , ]

eQTLtoolsResOther <- eQTLtoolsRes[!(eQTLtoolsRes$var %in% eQTLannotations[eQTLannotations$annot.type %in% c("hg19_genes_promoters", "hg19_enhancers_fantom"),]$var) , ]
methylQTLtoolsResOther <- methylQTLtoolsRes[!(methylQTLtoolsRes$var %in% methylQTLannotations[methylQTLannotations$annot.type %in% c("hg19_genes_promoters", "hg19_enhancers_fantom"),]$var) , ]

mean(eQTLtoolsResPromoter$slope)
mean(eQTLtoolsResEnhancer$slope)
mean(eQTLtoolsResOther$slope)

mean(methylQTLtoolsResPromoter$slope)
mean(methylQTLtoolsResEnhancer$slope)
mean(methylQTLtoolsResOther$slope)

# Distance from TSS histogram
eQTLtoolsResSig$distances_kb <- eQTLtoolsResSig$var_distance/1000

distance_eqtl <- ggplot(eQTLtoolsResSig, aes(x=distances_kb)) + geom_histogram(fill="coral1", binwidth=100) +
  theme_classic() +
  theme(axis.text.x = element_text(color="#000000", size=10),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_text(color="#000000", size=12),
        axis.title.y=element_text(color="#000000", size=12),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linetype = 1, lineend = "square")) +
  xlab("Distance from TSS") +
  ylab("Frequency")

distance_eqtl

methylQTLtoolsResSig$distances_kb <- methylQTLtoolsResSig$var_distance/1000

distance_methylqtl <- ggplot(methylQTLtoolsResSig, aes(x=distances_kb)) + geom_histogram(fill="dodgerblue2", binwidth=100) +
  theme_classic() +
  theme(axis.text.x = element_text(color="#000000", size=10),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_text(color="#000000", size=12),
        axis.title.y=element_text(color="#000000", size=12),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linetype = 1, lineend = "square")) +
  xlab("Distance from TSS") +
  ylab("Frequency")

distance_methylqtl

# Effect size histogram

effects_eqtl <- ggplot(eQTLtoolsResSig, aes(x=slope)) + geom_histogram(fill="coral1", binwidth=0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(color="#000000", size=10),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_text(color="#000000", size=12),
        axis.title.y=element_text(color="#000000", size=12),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linetype = 1, lineend = "square")) +
  xlab("Effect size") +
  ylab("Frequency")

effects_eqtl

effects_methylqtl <- ggplot(methylQTLtoolsResSig, aes(x=slope)) + geom_histogram(fill="dodgerblue2", binwidth=0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(color="#000000", size=10),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_text(color="#000000", size=12),
        axis.title.y=element_text(color="#000000", size=12),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linetype = 1, lineend = "square")) +
  xlab("Effect size") +
  ylab("Frequency")

effects_methylqtl

# Barplots of annotations
#Turn your 'annot.type' column into a character vector
eQTLannotationsCounts$annot.type <- as.character(eQTLannotationsCounts$annot.type)
#Then turn it back into a factor with the levels in the correct order
eQTLannotationsCounts$annot.type <- factor(eQTLannotationsCounts$annot.type, levels=unique(eQTLannotationsCounts$annot.type))
eQTLannotationsCounts$annot.type <- gsub("hg19_", "", eQTLannotationsCounts$annot.type)

# Grouped
eQTL_bars <- ggplot(eQTLannotationsCounts, aes(y=n, x=reorder(annot.type, -n))) +
  geom_bar(position="dodge", stat="identity", color="coral1", fill="coral1") + 
  #scale_fill_manual(values = c("azure3", "deepskyblue4", "firebrick")) +
  #scale_y_continuous(breaks=c(0,300,600,900,1200), limits = c(0,1200)) +
  theme_classic() +
  theme(axis.text.x = element_text(color="#000000", size=10, angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linetype = 1, lineend = "square"))


eQTL_bars

methylQTLannotationsCounts$annot.type <- as.character(methylQTLannotationsCounts$annot.type)
methylQTLannotationsCounts$annot.type <- factor(methylQTLannotationsCounts$annot.type, levels=unique(methylQTLannotationsCounts$annot.type))
methylQTLannotationsCounts$annot.type <- gsub("hg19_", "", methylQTLannotationsCounts$annot.type)

# Grouped
methylQTL_bars <- ggplot(methylQTLannotationsCounts, aes(y=n, x=reorder(annot.type, -n))) +
  geom_bar(position="dodge", stat="identity", color="dodgerblue2", fill="dodgerblue2") + 
  #scale_fill_manual(values = c("azure3", "deepskyblue4", "firebrick")) +
  #scale_y_continuous(breaks=c(0,300,600,900,1200), limits = c(0,1200)) +
  theme_classic() +
  theme(axis.text.x = element_text(color="#000000", size=10, angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linetype = 1, lineend = "square"))

methylQTL_bars

# Saving
ggsave("QTL_effects_distances_annotations.pdf", 
       grid.arrange(distance_methylqtl, distance_eqtl,
                    effects_methylqtl, effects_eqtl,
                    methylQTL_bars, eQTL_bars,
                    ncol=2, widths=c(4, 4)),
       width = 8, height = 8)

methylQTLannotationsCounts$group <- "methylQTL"
eQTLannotationsCounts$group <- "eQTL"

eQTL_methylQTLannotationsCounts <- rbind(eQTLannotationsCounts,
                                         methylQTLannotationsCounts)

eQTL_methylQTLannotationsCounts$group <- factor(eQTL_methylQTLannotationsCounts$group, levels=c("methylQTL", "eQTL"))

# Grouped
eQTL_methylQTL_bars <- ggplot(eQTL_methylQTLannotationsCounts, aes(y=n, x=reorder(annot.type, -n), color=group, fill=group)) +
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values = c("dodgerblue2", "coral1")) +
  scale_color_manual(values = c("dodgerblue2", "coral1")) +
  #scale_y_continuous(breaks=c(0,300,600,900,1200), limits = c(0,1200)) +
  theme_classic() +
  theme(axis.text.x = element_text(color="#000000", size=10, angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.ticks.length=unit(.25, "cm"),
        axis.ticks = element_line(linetype = 1, lineend = "square")) +
  ylim(0,1200)

eQTL_methylQTL_bars

ggsave("eQTL_methylQTL_annot_bars_bottom.pdf", 
       eQTL_methylQTL_bars,
       width = 8, height = 3)


# Plotting the effect sizes of the LAI deriven QTLs

pap_eqtl <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/eQTL_QTLfdrp010_genotype_varexplained_EASPAP.tsv")
deni_eqtl <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/eQTL_QTLfdrp010_genotype_varexplained_DENI.tsv")
nean_eqtl <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/eQTL_QTLfdrp010_genotype_varexplained_NEAN.tsv")
pap_methylqtl <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/methylQTL_QTLfdrp010_genotype_varexplained_EASPAP.tsv")
deni_methylqtl <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/methylQTL_QTLfdrp010_genotype_varexplained_DENI.tsv")
nean_methylqtl <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/methylQTL_QTLfdrp010_genotype_varexplained_NEAN.tsv")

pap_eqtl$var <- rownames(pap_eqtl)
deni_eqtl$var <- rownames(deni_eqtl)
nean_eqtl$var <- rownames(nean_eqtl)
pap_methylqtl$var <- rownames(pap_methylqtl)
deni_methylqtl$var <- rownames(deni_methylqtl)
nean_methylqtl$var <- rownames(nean_methylqtl)

# Getting QTL effect sizes of ancestry driven QTLs
pap_eqtl_sig <- pap_eqtl[pap_eqtl$fdrp<0.01,]
eQTLtoolsRes_pap <- eQTLtoolsRes[eQTLtoolsRes$var %in% rownames(pap_eqtl_sig),]
hist(abs(eQTLtoolsRes_pap$slope))

deni_eqtl_sig <- deni_eqtl[deni_eqtl$fdrp<0.01,]
eQTLtoolsRes_deni <- eQTLtoolsRes[eQTLtoolsRes$var %in% rownames(deni_eqtl_sig),]
hist(abs(eQTLtoolsRes_deni$slope))

nean_eqtl_sig <- nean_eqtl[nean_eqtl$fdrp<0.01,]
eQTLtoolsRes_nean <- eQTLtoolsRes[eQTLtoolsRes$var %in% rownames(nean_eqtl_sig),]
hist(abs(eQTLtoolsRes_nean$slope))

pap_methylqtl_sig <- pap_methylqtl[pap_methylqtl$fdrp<0.01,]
methylQTLtoolsRes_pap <- methylQTLtoolsRes[methylQTLtoolsRes$var %in% rownames(pap_methylqtl_sig),]
hist(abs(methylQTLtoolsRes_pap$slope))

deni_methylqtl_sig <- deni_methylqtl[deni_methylqtl$fdrp<0.01,]
methylQTLtoolsRes_deni <- methylQTLtoolsRes[methylQTLtoolsRes$var %in% rownames(deni_methylqtl_sig),]
hist(abs(methylQTLtoolsRes_deni$slope))

nean_methylqtl_sig <- nean_methylqtl[nean_methylqtl$fdrp<0.01,]
methylQTLtoolsRes_nean <- methylQTLtoolsRes[methylQTLtoolsRes$var %in% rownames(nean_methylqtl_sig),]
hist(abs(methylQTLtoolsRes_nean$slope))

# Allele frequency matching
# Adding MAF of each variant
maf <- read.table("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Colocalization/MAF_varID_removemissing.tsv")
colnames(maf) <- c("chr", "pos", "maf", "var")
maf$var <- paste("chr", maf$var, sep="")
eQTLtoolsRes_maf <- merge(eQTLtoolsRes, maf, by="var")
eQTLtoolsRes_pap_maf <- merge(eQTLtoolsRes_pap, maf, by="var")
eQTLtoolsRes_deni_maf <- merge(eQTLtoolsRes_deni, maf, by="var")
eQTLtoolsRes_nean_maf <- merge(eQTLtoolsRes_nean, maf, by="var")
methylQTLtoolsRes_maf <- merge(methylQTLtoolsRes, maf, by="var")
methylQTLtoolsRes_pap_maf <- merge(methylQTLtoolsRes_pap, maf, by="var")
methylQTLtoolsRes_deni_maf <- merge(methylQTLtoolsRes_deni, maf, by="var")
methylQTLtoolsRes_nean_maf <- merge(methylQTLtoolsRes_nean, maf, by="var")

eQTLtoolsRes_maf_notarc <- eQTLtoolsRes_maf[!(eQTLtoolsRes_maf$var %in% c(eQTLtoolsRes_deni$var, eQTLtoolsRes_nean$var)),]
methylQTLtoolsRes_maf_notarc <- methylQTLtoolsRes_maf[!(methylQTLtoolsRes_maf$var %in% c(methylQTLtoolsRes_deni$var, methylQTLtoolsRes_nean$var)),]

mean(abs(eQTLtoolsRes_maf_notarc$slope))
mean(abs(eQTLtoolsRes_pap$slope))
mean(abs(eQTLtoolsRes_deni$slope))
mean(abs(eQTLtoolsRes_nean$slope))

mean(abs(methylQTLtoolsRes_maf_notarc$slope))
mean(abs(methylQTLtoolsRes_pap$slope))
mean(abs(methylQTLtoolsRes_deni$slope))
mean(abs(methylQTLtoolsRes_nean$slope))

t.test(abs(methylQTLtoolsRes_maf_notarc$slope),abs(methylQTLtoolsRes_deni_maf$slope))
t.test(abs(methylQTLtoolsRes_maf_notarc$slope),abs(methylQTLtoolsRes_nean_maf$slope)) 

t.test(abs(eQTLtoolsRes_maf_notarc$slope),abs(eQTLtoolsRes_deni_maf$slope))
t.test(abs(eQTLtoolsRes_maf_notarc$slope),abs(eQTLtoolsRes_nean_maf$slope)) 

mean(methylQTLtoolsRes_maf_notarc$maf)
sd(methylQTLtoolsRes_maf_notarc$maf)
mean(methylQTLtoolsRes_deni_maf$maf)
sd(methylQTLtoolsRes_deni_maf$maf)
mean(methylQTLtoolsRes_nean_maf$maf)
sd(methylQTLtoolsRes_nean_maf$maf)

mean(eQTLtoolsRes_maf_notarc$maf)
sd(eQTLtoolsRes_maf_notarc$maf)
mean(eQTLtoolsRes_deni_maf$maf)
sd(eQTLtoolsRes_deni_maf$maf)
mean(eQTLtoolsRes_nean_maf$maf)
sd(eQTLtoolsRes_nean_maf$maf)


# Subsetting matched MAFs
# eQTL deni
eQTLtoolsRes_deni_maf$pop <- 1
eQTLtoolsRes_maf_notarc$pop <- 0
eQTL_noarc_deni <- rbind.fill(eQTLtoolsRes_deni_maf, eQTLtoolsRes_maf_notarc)
eQTL_noarc_deni <- eQTL_noarc_deni[,c("maf", "var", "pop")]
eQTL_deni_match <- matchit(pop ~ maf, data = eQTL_noarc_deni, method="nearest")
eQTL_noarc_deni_matched <- match.data(eQTL_deni_match)
#eQTL_noarc_deni_matched[eQTL_noarc_deni_matched$pop==1,]$var
eQTL_noarc_deni_matched_deni_slopes <- eQTLtoolsRes_deni_maf[eQTLtoolsRes_deni_maf$var %in% eQTL_noarc_deni_matched[eQTL_noarc_deni_matched$pop==1,]$var , ]
eQTL_noarc_deni_matched_nodeni_slopes <- eQTLtoolsRes_maf_notarc[eQTLtoolsRes_maf_notarc$var %in% eQTL_noarc_deni_matched[eQTL_noarc_deni_matched$pop==0,]$var , ]

# eQTL nean
eQTLtoolsRes_nean_maf$pop <- 1
eQTLtoolsRes_maf_notarc$pop <- 0
eQTL_noarc_nean <- rbind.fill(eQTLtoolsRes_nean_maf, eQTLtoolsRes_maf_notarc)
eQTL_noarc_nean <- eQTL_noarc_nean[,c("maf", "var", "pop")]
eQTL_nean_match <- matchit(pop ~ maf, data = eQTL_noarc_nean, method="nearest")
eQTL_noarc_nean_matched <- match.data(eQTL_nean_match)
#eQTL_noarc_deni_matched[eQTL_noarc_deni_matched$pop==1,]$var
eQTL_noarc_nean_matched_nean_slopes <- eQTLtoolsRes_nean_maf[eQTLtoolsRes_nean_maf$var %in% eQTL_noarc_nean_matched[eQTL_noarc_nean_matched$pop==1,]$var , ]
eQTL_noarc_nean_matched_nonean_slopes <- eQTLtoolsRes_maf_notarc[eQTLtoolsRes_maf_notarc$var %in% eQTL_noarc_nean_matched[eQTL_noarc_nean_matched$pop==0,]$var , ]

t.test(eQTL_noarc_nean_matched_nean_slopes$maf, eQTL_noarc_nean_matched_nonean_slopes$maf)

# methylQTL deni
methylQTLtoolsRes_deni_maf$pop <- 1
methylQTLtoolsRes_maf_notarc$pop <- 0
methylQTL_noarc_deni <- rbind.fill(methylQTLtoolsRes_deni_maf, methylQTLtoolsRes_maf_notarc)
methylQTL_noarc_deni <- methylQTL_noarc_deni[,c("maf", "var", "pop")]
methylQTL_deni_match <- matchit(pop ~ maf, data = methylQTL_noarc_deni, method="nearest")
methylQTL_noarc_deni_matched <- match.data(methylQTL_deni_match)
#methylQTL_noarc_deni_matched[methylQTL_noarc_deni_matched$pop==1,]$var
methylQTL_noarc_deni_matched_deni_slopes <- methylQTLtoolsRes_deni_maf[methylQTLtoolsRes_deni_maf$var %in% methylQTL_noarc_deni_matched[methylQTL_noarc_deni_matched$pop==1,]$var , ]
methylQTL_noarc_deni_matched_nodeni_slopes <- methylQTLtoolsRes_maf_notarc[methylQTLtoolsRes_maf_notarc$var %in% methylQTL_noarc_deni_matched[methylQTL_noarc_deni_matched$pop==0,]$var , ]

t.test(methylQTL_noarc_deni_matched_deni_slopes$maf, methylQTL_noarc_deni_matched_nodeni_slopes$maf)

# methylQTL nean
methylQTLtoolsRes_nean_maf$pop <- 1
methylQTLtoolsRes_maf_notarc$pop <- 0
methylQTL_noarc_nean <- rbind.fill(methylQTLtoolsRes_nean_maf, methylQTLtoolsRes_maf_notarc)
methylQTL_noarc_nean <- methylQTL_noarc_nean[,c("maf", "var", "pop")]
methylQTL_nean_match <- matchit(pop ~ maf, data = methylQTL_noarc_nean, method = "nearest")
methylQTL_noarc_nean_matched <- match.data(methylQTL_nean_match)
#methylQTL_noarc_deni_matched[methylQTL_noarc_deni_matched$pop==1,]$var
methylQTL_noarc_nean_matched_nean_slopes <- methylQTLtoolsRes_nean_maf[methylQTLtoolsRes_nean_maf$var %in% methylQTL_noarc_nean_matched[methylQTL_noarc_nean_matched$pop==1,]$var , ]
methylQTL_noarc_nean_matched_nonean_slopes <- methylQTLtoolsRes_maf_notarc[methylQTLtoolsRes_maf_notarc$var %in% methylQTL_noarc_nean_matched[methylQTL_noarc_nean_matched$pop==0,]$var , ]

t.test(methylQTL_noarc_nean_matched_nean_slopes$maf, methylQTL_noarc_nean_matched_nonean_slopes$maf)
t.test(methylQTL_noarc_nean_matched_nean_slopes$slope, methylQTL_noarc_nean_matched_nonean_slopes$slope)

# A dataframe for plotting
# eQTLs
eQTLtoolsRes_maf_notarc$pop <- "Other"
eQTLtoolsRes_deni_maf$pop <- "Denisovan driven"
eQTLtoolsRes_nean_maf$pop <- "Neanderthal driven"

#eQTL_noarc_deni_matched_deni_slopes$pop <- "deni"
eQTL_noarc_deni_matched_nodeni_slopes$pop <- "Other, MAF matched"
#eQTL_noarc_deni_matched_nean_slopes$pop <- "nean"
eQTL_noarc_nean_matched_nonean_slopes$pop <- "Other, MAF matched"

methylQTLtoolsRes_maf_notarc$pop <- "Other"
methylQTLtoolsRes_deni_maf$pop <- "Denisovan driven"
methylQTLtoolsRes_nean_maf$pop <- "Neanderthal driven"

#methylQTL_noarc_deni_matched_deni_slopes$pop <- "deni"
methylQTL_noarc_deni_matched_nodeni_slopes$pop <- "Other, MAF matched"
#methylQTL_noarc_nean_matched_nean_slopes$pop <- "nean"
methylQTL_noarc_nean_matched_nonean_slopes$pop <- "Other, MAF matched"

#methylQTL_noarc_deni_matched_deni_slopes$abseffect <- abs(methylQTL_noarc_deni_matched_deni_slopes$slope)
#methylQTL_noarc_deni_matched_nodeni_slopes$abseffect <- abs(methylQTL_noarc_deni_matched_nodeni_slopes$slope)
#t.test(methylQTL_noarc_deni_matched_deni_slopes$abseffect,methylQTL_noarc_deni_matched_nodeni_slopes$abseffect)

methylqtl_notmatched_plotdata <- rbind.fill(methylQTLtoolsRes_maf_notarc, methylQTLtoolsRes_deni_maf, methylQTLtoolsRes_nean_maf)
methylqtl_notmatched_plotdata$abseffect <- abs(methylqtl_notmatched_plotdata$slope)
methylqtl_deni_matched_plotdata <- rbind.fill(methylQTL_noarc_deni_matched_nodeni_slopes, methylQTLtoolsRes_deni_maf)
methylqtl_deni_matched_plotdata$abseffect <- abs(methylqtl_deni_matched_plotdata$slope)
methylqtl_nean_matched_plotdata <- rbind.fill(methylQTL_noarc_nean_matched_nonean_slopes, methylQTLtoolsRes_nean_maf)
methylqtl_nean_matched_plotdata$abseffect <- abs(methylqtl_nean_matched_plotdata$slope)
#t.test(abs(noarc_deni_matched_d_slopes$slope),abs(noarc_deni_matched_no_slopes$slope))

t.test(methylQTL_noarc_deni_matched_deni_slopes$maf, methylQTL_noarc_deni_matched_nodeni_slopes$maf)

# Visualizing
# Colors
palette <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2")

# Specify the comparisons
my_comparisons <- list( c("Other", "Denisovan driven"), c("Other", "Neanderthal driven"), c("Denisovan driven", "Neanderthal driven") )

# For the nonmatched
methyl_effect <- ggboxplot(methylqtl_notmatched_plotdata, x = "pop", y = "abseffect", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("Absolute effect size") +
  rremove("legend")
  #stat_compare_means(label.y = 50)     # Add global p-value

methyl_maf <- ggboxplot(methylqtl_notmatched_plotdata, x = "pop", y = "maf", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("MAF") +
  rremove("legend")

# For the Denisovan matched
my_comparisons <- list( c("Other, MAF matched", "Denisovan driven") )
matched_deni_methyl_effect <- ggboxplot(methylqtl_deni_matched_plotdata, x = "pop", y = "abseffect", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("Absolute effect size") +
  rremove("legend")
#stat_compare_means(label.y = 50)     # Add global p-value

matched_deni_methyl_maf <- ggboxplot(methylqtl_deni_matched_plotdata, x = "pop", y = "maf", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("MAF") +
  rremove("legend")

# For the Neanderthal matched
my_comparisons <- list( c("Other, MAF matched", "Neanderthal driven") )

matched_nean_methyl_effect <- ggboxplot(methylqtl_nean_matched_plotdata, x = "pop", y = "abseffect", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("Absolute effect size") +
  rremove("legend")
#stat_compare_means(label.y = 50)     # Add global p-value

matched_nean_methyl_maf <- ggboxplot(methylqtl_nean_matched_plotdata, x = "pop", y = "maf", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("MAF") +
  rremove("legend")



ggsave("methylQTL_effects_MAFs.pdf", 
       grid.arrange(methyl_maf, methyl_effect,
                    matched_deni_methyl_maf, matched_deni_methyl_effect,
                    matched_nean_methyl_maf, matched_nean_methyl_effect,
                    ncol=2, widths=c(5, 5)),
                    width = 10, height = 15)


# Plotting eQTLs
# Specify the comparisons
eqtl_notmatched_plotdata <- rbind.fill(eQTLtoolsRes_maf_notarc, eQTLtoolsRes_deni_maf, eQTLtoolsRes_nean_maf)
eqtl_notmatched_plotdata$abseffect <- abs(eqtl_notmatched_plotdata$slope)
eqtl_deni_matched_plotdata <- rbind.fill(eQTL_noarc_deni_matched_nodeni_slopes, eQTLtoolsRes_deni_maf)
eqtl_deni_matched_plotdata$abseffect <- abs(eqtl_deni_matched_plotdata$slope)
eqtl_nean_matched_plotdata <- rbind.fill(eQTL_noarc_nean_matched_nonean_slopes, eQTLtoolsRes_nean_maf)
eqtl_nean_matched_plotdata$abseffect <- abs(eqtl_nean_matched_plotdata$slope)

my_comparisons <- list( c("Other", "Denisovan driven"), c("Other", "Neanderthal driven"), c("Denisovan driven", "Neanderthal driven") )

# For the nonmatched
e_effect <- ggboxplot(eqtl_notmatched_plotdata, x = "pop", y = "abseffect", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("Absolute effect size") +
  rremove("legend")
#stat_compare_means(label.y = 50)     # Add global p-value

e_maf <- ggboxplot(eqtl_notmatched_plotdata, x = "pop", y = "maf", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("MAF") +
  rremove("legend")

# For the Denisovan matched
my_comparisons <- list( c("Other, MAF matched", "Denisovan driven") )
matched_deni_e_effect <- ggboxplot(eqtl_deni_matched_plotdata, x = "pop", y = "abseffect", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("Absolute effect size") +
  rremove("legend")
#stat_compare_means(label.y = 50)     # Add global p-value

matched_deni_e_maf <- ggboxplot(eqtl_deni_matched_plotdata, x = "pop", y = "maf", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("MAF") +
  rremove("legend")

# For the Neanderthal matched
my_comparisons <- list( c("Other, MAF matched", "Neanderthal driven") )

matched_nean_e_effect <- ggboxplot(eqtl_nean_matched_plotdata, x = "pop", y = "abseffect", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("Absolute effect size") +
  rremove("legend")
#stat_compare_means(label.y = 50)     # Add global p-value

matched_nean_e_maf <- ggboxplot(eqtl_nean_matched_plotdata, x = "pop", y = "maf", color = "pop", palette = palette, add="jitter") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
  #title("") +
  xlab("Type of QTL") +
  ylab("MAF") +
  rremove("legend")


ggsave("eQTL_effects_MAFs.pdf", 
       grid.arrange(e_maf, e_effect,
                    matched_deni_e_maf, matched_deni_e_effect,
                    matched_nean_e_maf, matched_nean_e_effect,
                    ncol=2, widths=c(5, 5)),
       width = 10, height = 15)






