# ==============================================================================
# Author(s) : Heini M Natri, heini.natri@gmail.com
# Date: June 2020
# Description: Calculating the variance in QTL genotype explained by local
# ancestry
# ==============================================================================

# ======================================
# Load libraries
# ======================================

library(GenomicRanges)
library(wesanderson)
library(ggplot2)
library(grid)
library(gridExtra)

# ======================================
# Environment parameters
# ======================================

# Working directory
setwd("/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/")

# Colors
# Mappi
MPIcol <- wes_palette("Zissou1", 20, type = "continuous")[20]
# Mentawai
MTWcol <- wes_palette("Zissou1", 20, type = "continuous")[1]
# Sumba
SMBcol <- wes_palette("Zissou1", 20, type = "continuous")[11]

# ======================================
# Data
# ======================================

# exp or methyl
dataType <- "exp"

QTLtoolsResPath <- NULL
QTLALTallelesPath <- NULL # Genotypes in QTLs, FDR-p<0.1, subset_qtl_vars.sh
molPhenoDataPath <- NULL
outputPrefix <- NULL
LAIQTLresEASPAPpath <- NULL
LAIQTLresDENIpath <- NULL
LAIQTLresNEANpath <- NULL

# TODO: use nominal pass?
if (dataType=="exp"){
  QTLtoolsResPath <- "/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indonesian_eQTL_lifted_ALL_CHRS_perm10k_FDR010_significant.txt"
  QTLALTallelesPath <- "eQTLs_FDRp010.GT.FORMAT"
  molPhenoDataPath <- "/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/RNAseq/lifted_hg19_n123_fpkm_filtered_normalized_sample_gene_for_QTLtools_geneID.tsv"
  outputPrefix <- "eQTL"
  LAIQTLresEASPAPpath <- "eQTL_QTLfdrp001_genotype_varexplained_EASPAP.tsv"
  LAIQTLresDENIpath <- "eQTL_QTLfdrp001_genotype_varexplained_DENI.tsv"
  LAIQTLresNEANpath <- "eQTL_QTLfdrp001_genotype_varexplained_NEAN.tsv"
} else {
  QTLtoolsResPath <- "/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/result_29peer_5gt/Indonesian_methylQTL_ALL_CHUNKS_perm10k_FDR001_significant.txt"
  QTLALTallelesPath <- "methylQTLs_FDRp010.GT.FORMAT"
  molPhenoDataPath <- "/Users/hnatri/Dropbox (Personal)/Indonesian_eQTL/Ancestry/methylQTL_mvals_plotvartargets_wheader.bed"
  outputPrefix <- "methylQTL"
  LAIQTLresEASPAPpath <- "methylQTL_QTLfdrp001_genotype_varexplained_EASPAP.tsv"
  LAIQTLresDENIpath <- "methylQTL_QTLfdrp001_genotype_varexplained_DENI.tsv"
  LAIQTLresNEANpath <- "methylQTL_QTLfdrp001_genotype_varexplained_NEAN.tsv"
}


# ======================================
# Helper functions
# ======================================

lmp <- function (lm) {
  if (class(lm) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(lm)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

get_legend <- function(plot){
  tmp <- ggplot_gtable(ggplot_build(p1))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# ======================================
# Modern ancestry
# ======================================

# Selecting QTLs to test
QTLtoolsRes <- read.table(QTLtoolsResPath)
colnames(QTLtoolsRes) <- c("molecular_trait_id", "chromosome", "target_start",
                       "target_end", "strand", "length", "var_distance",
                       "var_rsID", "var_chromosome", "var_start", "var_end",
                       "dfs", "dummy1", "1st_param_betadist", "2nd_param_betadist", "nominalp",
                       "slope", "adjp", "adjp_beta", "fdrp", "fdrp1")
QTLtoolsRes$var <- paste(QTLtoolsRes$var_chromosome, QTLtoolsRes$var_start, sep="_")
QTLtoolsRes$var <- paste("chr", QTLtoolsRes$var, sep="")
QTLtoolsResSig <- QTLtoolsRes[QTLtoolsRes$fdrp<0.1 , ]

# Significant (FDR-p<0.01) eQTL genotypes extracted from the VCF
# Dataframe with N of ALT alleles in each QTL genotype for each individual
QTLALTalleles <- read.table(QTLALTallelesPath, header=T)
QTLALTalleles$CHROM <- paste("chr", QTLALTalleles$CHROM, sep="")
rownames(QTLALTalleles) <- paste(QTLALTalleles$CHROM, QTLALTalleles$POS, sep = "_")

QTLpositions <- data.frame(matrix(NA, nrow = nrow(QTLALTalleles), ncol = 3))
colnames(QTLpositions) <- c("chr", "start", "stop")
chr_pos <- strsplit(rownames(QTLALTalleles), "_")
chr <- unlist(chr_pos)[2*(1:length(rownames(QTLALTalleles)))-1]
pos <- unlist(chr_pos)[2*(1:length(rownames(QTLALTalleles)))  ]
QTLpositions$chr <- chr
QTLpositions$start <- as.numeric(pos)
QTLpositions$stop <- as.numeric(pos)
# Converting to GRanges
QTLpositions_GRanges <- makeGRangesFromDataFrame(QTLpositions,
                                                   seqinfo=NULL,
                                                   seqnames.field="chr",
                                                   start.field="start",
                                                   end.field=, "stop",
                                                   ignore.strand=TRUE,
                                                   keep.extra.columns=FALSE)

# Formatting
QTLALTalleles[] <- lapply(QTLALTalleles, as.character)
QTLALTalleles[QTLALTalleles=="0|0"] <- 0
QTLALTalleles[QTLALTalleles=="1|0"] <- 1
QTLALTalleles[QTLALTalleles=="0|1"] <- 1
QTLALTalleles[QTLALTalleles=="1|1"] <- 2
QTLALTalleles[QTLALTalleles=="."] <- NA
QTLALTalleles <- QTLALTalleles[ , !(colnames(QTLALTalleles) %in% c("CHROM", "POS"))]

# Samples with LAI
individuals <- readLines("LAI_individuals.txt")

# Dataframes with N of PAP an EAS alleles in each QTL position for each individual
N_PAP <- data.frame(matrix(NA, nrow = nrow(QTLALTalleles), ncol = length(individuals)))
colnames(N_PAP) <- individuals
rownames(N_PAP) <- rownames(QTLALTalleles)
N_EAS <- data.frame(matrix(NA, nrow = nrow(QTLALTalleles), ncol = length(individuals)))
colnames(N_EAS) <- individuals
rownames(N_EAS) <- rownames(QTLALTalleles)

find_overlaps_modern <- function(ind, pop) {
  hap1 <- read.table(paste("heini_wgs_cp/", ind, ".", pop, ".1", sep=""), col.names=c("chr", "start", "stop"))
  # Converting to GRanges
  hap1_GRanges <- makeGRangesFromDataFrame(hap1,
                                           seqinfo=NULL,
                                           seqnames.field="chr",
                                           start.field="start",
                                           end.field=, "stop",
                                           ignore.strand=TRUE,
                                           keep.extra.columns=FALSE)
  
  # Overlapping the significant QTL positions with the ancestry informative regions
  Sub <- subsetByOverlaps(QTLpositions_GRanges,hap1_GRanges)
  QTLs_overlapping_anc_alleles_1 <- data.frame(seqnames=seqnames(Sub),
                                                starts=start(Sub),
                                                ends=end(Sub))
  vector_of_overlapping_anc_alleles_1 <- paste(QTLs_overlapping_anc_alleles_1$seqnames,
                                               QTLs_overlapping_anc_alleles_1$start,
                                               sep = "_")
  
  
  hap2 <- read.table(paste("heini_wgs_cp/", ind, ".", pop, ".2", sep=""), col.names=c("chr", "start", "stop"))
  # Converting to GRanges
  hap2_GRanges <- makeGRangesFromDataFrame(hap2,
                                           seqinfo=NULL,
                                           seqnames.field="chr",
                                           start.field="start",
                                           end.field=, "stop",
                                           ignore.strand=TRUE,
                                           keep.extra.columns=FALSE)
  
  # Overlapping the significant QTL positions with the EAS, PAP, DENI, and NEAN regions
  # Findign overlaps
  Sub <- subsetByOverlaps(QTLpositions_GRanges,hap2_GRanges)
  QTLs_overlapping_anc_alleles_2 <- data.frame(seqnames=seqnames(Sub),
                                                starts=start(Sub),
                                                ends=end(Sub))
  vector_of_overlapping_anc_alleles_2 <- paste(QTLs_overlapping_anc_alleles_2$seqnames,
                                               QTLs_overlapping_anc_alleles_2$start,
                                               sep = "_")
  return(list(vector_of_overlapping_anc_alleles_1=vector_of_overlapping_anc_alleles_1,
              vector_of_overlapping_anc_alleles_2=vector_of_overlapping_anc_alleles_2))

}

# Iterating through all individuals, counting N of alleles
for (ind in individuals) {
  # Finding vectors of overlapping EAS and PAP alleles
  vector_of_overlapping_EAS_alleles_1 <- find_overlaps_modern(ind, "eas")[[1]]
  vector_of_overlapping_EAS_alleles_2 <- find_overlaps_modern(ind, "eas")[[2]]
  vector_of_overlapping_PAP_alleles_1 <- find_overlaps_modern(ind, "pap")[[1]]
  vector_of_overlapping_PAP_alleles_2 <- find_overlaps_modern(ind, "pap")[[2]]

  #intersect(vector_of_overlapping_anc_alleles_2,vector_of_overlapping_anc_alleles_1)
  
  # Counting EAS alleles in each significant QTL position
  for (var in rownames(QTLALTalleles)){
    N_EAS_alleles <- 0
    if (var %in% vector_of_overlapping_EAS_alleles_1){
      N_EAS_alleles <- N_EAS_alleles+1
    }
    if (var %in% vector_of_overlapping_EAS_alleles_2){
      N_EAS_alleles <- N_EAS_alleles+1
    }
    # Setting the N of alleles to the result dataframe
    N_EAS[var, ind] <- N_EAS_alleles
  }
  
  # Counting PAP alleles in each significant QTL position
  for (var in rownames(QTLALTalleles)){
    N_PAP_alleles <- 0
    if (var %in% vector_of_overlapping_PAP_alleles_1){
      N_PAP_alleles <- N_PAP_alleles+1
    }
    if (var %in% vector_of_overlapping_PAP_alleles_2){
      N_PAP_alleles <- N_PAP_alleles+1
    }
    # Setting the N of alleles to the result dataframe
    N_PAP[var, ind] <- N_PAP_alleles
  }
}

# Fixing the IDs of some individuals
head(N_PAP)
colnames(N_PAP) <- gsub("MTW0", "MTW-0", colnames(N_PAP))
colnames(N_EAS) <- gsub("MTW0", "MTW-0", colnames(N_EAS))

# Saving to a file
#write.table(N_PAP, "N_PAP_methylQTL_FDRp001_sites.tsv")
#write.table(N_EAS, "N_EAS_methylQTL_FDRp001_sites.tsv")
#N_PAP <- read.table("N_PAP_methylQTL_FDRp010_sites.tsv")
#N_EAS <- read.table("N_EAS_methylQTL_FDRp010_sites.tsv")
colnames(N_PAP) <- gsub("\\.", "-", colnames(N_PAP))
colnames(N_EAS) <- gsub("\\.", "-", colnames(N_EAS))

# Subsetting WGS samples from the QTL ALT-allele dataframes
head(QTLALTalleles)
colnames(QTLALTalleles) <- gsub("\\.", "-", colnames(QTLALTalleles))
QTLALTalleles <- QTLALTalleles[ , colnames(QTLALTalleles) %in% colnames(N_PAP)]
# Sorting columns to match
N_PAP <- N_PAP[ , colnames(N_PAP) %in% colnames(QTLALTalleles)]
N_PAP <- N_PAP[ , match(colnames(QTLALTalleles), colnames(N_PAP))]
N_EAS <- N_EAS[ , colnames(N_EAS) %in% colnames(QTLALTalleles)]
N_EAS <- N_EAS[ , match(colnames(QTLALTalleles), colnames(N_EAS))]

# Result dataframe
LAIQTLres <- data.frame(matrix(NA, nrow = length(unique(QTLtoolsResSig$var)), ncol = 4))
rownames(LAIQTLres) <- unique(QTLtoolsResSig$var)
colnames(LAIQTLres) <- c("R2", "adjR2", "p", "target")

# Running linear regression for each variant
for (var in rownames(LAIQTLres)) {
  V <- as.numeric(QTLALTalleles[var,])
  PAP <- as.numeric(N_PAP[var,])
  EAS <- as.numeric(N_EAS[var,])
  if (length(unique(PAP))==1 | length(unique(EAS))==1){
    next
  }
  lm <- lm(V ~ PAP + EAS)
  # Get the R^2 and p-values
  R2 <- summary(lm)$r.squared
  adjR2 <- summary(lm)$adj.r.squared
  
  p <- lmp(lm)
  
  # Finding target
  target <- QTLtoolsResSig[QTLtoolsResSig$var == var , ]$molecular_trait_id
  if (length(target)>1){
    target <- paste(target, collapse = ", ")
  }
  
  # Saving results
  LAIQTLres[var, "R2"] <- R2
  LAIQTLres[var, "adjR2"] <- adjR2
  LAIQTLres[var, "p"] <- p
  LAIQTLres[var, "target"] <- as.character(target)
}

# Removing NAs
#LAIQTLres <- LAIQTLres[complete.cases(LAIQTLres) , ]
head(LAIQTLres)
LAIQTLres <- LAIQTLres[complete.cases(LAIQTLres) , ]
LAIQTLres$fdrp <- p.adjust(LAIQTLres$p, method = "fdr", n=length(LAIQTLres$p))
#LAIQTLres_sig <- LAIQTLres[LAIQTLres$fdrp<0.01 , ]
LAIQTLres_sig <- LAIQTLres[abs(LAIQTLres$R2)>0.7 , ]
dim(LAIQTLres_sig)

# Writing results to a file
write.table(LAIQTLres, LAIQTLresEASPAPpath, sep="\t")
LAIQTLres <- read.table(LAIQTLresEASPAPpath, sep="\t")
LAIQTLres_sig <- LAIQTLres[abs(LAIQTLres$R2)>0.7 , ]

# ======================================
# Archaic ancestry
# ======================================

find_overlaps_archaic <- function(ind, pop) {
  hap1 <- read.table(paste("hconf_61Af_binary_trim_hmm_skov_cparchaicalt_20042020/", ind, "_", pop, "_H1.ssAf61.gapmask.highconf_HMMSkov001_CPAlt999.bed", sep=""), col.names=c("chr", "start", "stop"))
  # Converting to GRanges
  hap1_GRanges <- makeGRangesFromDataFrame(hap1,
                                           seqinfo=NULL,
                                           seqnames.field="chr",
                                           start.field="start",
                                           end.field=, "stop",
                                           ignore.strand=TRUE,
                                           keep.extra.columns=FALSE)
  
  # Overlapping the significant QTL positions with the ancestry informative regions
  Sub <- subsetByOverlaps(QTLpositions_GRanges,hap1_GRanges)
  QTLs_overlapping_anc_alleles_1 <- data.frame(seqnames=seqnames(Sub),
                                               starts=start(Sub),
                                               ends=end(Sub))
  vector_of_overlapping_anc_alleles_1 <- paste(QTLs_overlapping_anc_alleles_1$seqnames,
                                               QTLs_overlapping_anc_alleles_1$starts,
                                               sep = "_")
  
  
  hap2 <- read.table(paste("hconf_61Af_binary_trim_hmm_skov_cparchaicalt_20042020/", ind, "_", pop, "_H2.ssAf61.gapmask.highconf_HMMSkov001_CPAlt999.bed", sep=""), col.names=c("chr", "start", "stop"))
  # Converting to GRanges
  hap2_GRanges <- makeGRangesFromDataFrame(hap2,
                                           seqinfo=NULL,
                                           seqnames.field="chr",
                                           start.field="start",
                                           end.field=, "stop",
                                           ignore.strand=TRUE,
                                           keep.extra.columns=FALSE)
  
  # Overlapping the significant QTL positions with the ancestry informartive regions
  # Findign overlaps
  Sub <- subsetByOverlaps(QTLpositions_GRanges,hap2_GRanges)
  QTLs_overlapping_anc_alleles_2 <- data.frame(seqnames=seqnames(Sub),
                                               starts=start(Sub),
                                               ends=end(Sub))
  vector_of_overlapping_anc_alleles_2 <- paste(QTLs_overlapping_anc_alleles_2$seqnames,
                                               QTLs_overlapping_anc_alleles_2$starts,
                                               sep = "_")
  return(list(vector_of_overlapping_anc_alleles_1=vector_of_overlapping_anc_alleles_1,
              vector_of_overlapping_anc_alleles_2=vector_of_overlapping_anc_alleles_2))
  
}

# Iterating through all individuals, counting N of introgressed Denisovan and 
# Neanderthal alleles
get_arc_alleles <- function(arc) {
  # An empty df for N of archaic introgressed alleles in each position and individual
  N_ARC <- data.frame(matrix(NA, nrow = nrow(QTLALTalleles), ncol = length(individuals)))
  rownames(N_ARC) <- rownames(QTLALTalleles)
  colnames(N_ARC) <- individuals
  for (ind in individuals) {
    # Finding vectors of overlapping EAS and PAP alleles
    vector_of_overlapping_arc_alleles_1 <- find_overlaps_archaic(ind, arc)[[1]]
    vector_of_overlapping_arc_alleles_2 <- find_overlaps_archaic(ind, arc)[[2]]
    
    # Counting archaic alleles in each significant QTL position
    for (var in rownames(QTLALTalleles)){
      N_alleles <- 0
      if (var %in% vector_of_overlapping_arc_alleles_1){
        N_alleles <- N_alleles+1
      }
      if (var %in% vector_of_overlapping_arc_alleles_2){
        N_alleles <- N_alleles+1
      }
      # Setting the N of alleles to the result dataframe
      N_ARC[var, ind] <- N_alleles
    }
  }
  return (N_ARC)
}


# Running linear regression for each variant
get_arc <- function(arc, N_ARC) {
  # Result dataframe
  LAIQTLres_arc <- data.frame(matrix(NA, nrow = length(unique(QTLtoolsResSig$var)), ncol = 5))
  rownames(LAIQTLres_arc) <- unique(QTLtoolsResSig$var)
  colnames(LAIQTLres_arc) <- c("R2", "adjR2", "p", "target", "note")

  for (var in rownames(LAIQTLres_arc)) {
    V <- as.numeric(QTLALTalleles[var,])
    ARC <- as.numeric(N_ARC[var,])
    #NOT_ARC <- as.numeric(chartr("012","210", ARC))
    if (length(unique(ARC))==1){
      LAIQTLres_arc[var, "note"] <- "monomorphic"
      next
    } 
    
    lm <- lm(V ~ ARC)
    # Get the R^2 and p-values
    R2 <- summary(lm)$r.squared
    adjR2 <- summary(lm)$adj.r.squared
    
    p <- lmp(lm)
    
    # Finding target
    target <- QTLtoolsResSig[QTLtoolsResSig$var == var , ]$molecular_trait_id
    if (length(target)>1){
      target <- paste(target, collapse = ", ")
    }
    
    # Saving results
    LAIQTLres_arc[var, "R2"] <- R2
    LAIQTLres_arc[var, "adjR2"] <- adjR2
    LAIQTLres_arc[var, "p"] <- p
    LAIQTLres_arc[var, "target"] <- as.character(target)
    LAIQTLres_arc[var, "note"] <- "polymorphic"
  }
  return (LAIQTLres_arc)
}

# Getting N of archaic introgressed alleles for each variant
N_DENI <- get_arc_alleles("D")
# Fixing the IDs of some individuals
colnames(N_DENI) <- gsub("MTW0", "MTW-0", colnames(N_DENI))
# Sorting columns to match the QTL genotypes
N_DENI <- N_DENI[ , colnames(N_DENI) %in% colnames(QTLALTalleles)]
N_DENI <- N_DENI[ , match(colnames(QTLALTalleles), colnames(N_DENI))]

# Getting results for Denisovan ancestry
LAIQTLresDENI <- get_arc("D", N_DENI)
#LAIQTLresDENI <- LAIQTLres_arc
head(LAIQTLresDENI)
# Removing NAs
LAIQTLresDENI <- LAIQTLresDENI[complete.cases(LAIQTLresDENI) , ]
head(LAIQTLresDENI)
dim(LAIQTLresDENI)
LAIQTLresDENI$fdrp <- p.adjust(LAIQTLresDENI$p, method = "fdr", n=length(LAIQTLresDENI$p))
#LAIQTLresDENIsig <- LAIQTLresDENI[LAIQTLresDENI$fdrp<0.01 , ]
LAIQTLresDENIsig <- LAIQTLresDENI[abs(LAIQTLresDENI$R2)>0.7 , ]
dim(LAIQTLresDENIsig)

# chr1_33231380
#N_DENI["chr1_33231380",]

# Same for Neanderthal ancestry

# Getting N of archaic introgressed alleles for each variant
N_NEAN <- get_arc_alleles("N")
# Fixing the IDs of some individuals
colnames(N_NEAN) <- gsub("MTW0", "MTW-0", colnames(N_NEAN))
# Sorting columns to match the QTL genotypes
N_NEAN <- N_NEAN[ , colnames(N_NEAN) %in% colnames(QTLALTalleles)]
N_NEAN <- N_NEAN[ , match(colnames(QTLALTalleles), colnames(N_NEAN))]

#write.table(N_NEAN, "N_NEAN_eQTL_FDRp001_sites.tsv")

LAIQTLresNEAN <- get_arc("N", N_NEAN)

# Removing NAs
LAIQTLresNEAN <- LAIQTLresNEAN[complete.cases(LAIQTLresNEAN) , ]
LAIQTLresNEAN$fdrp <- p.adjust(LAIQTLresNEAN$p, method = "fdr", n=length(LAIQTLresNEAN$p))
#LAIQTLresNEANsig <- LAIQTLresNEAN[LAIQTLresNEAN$fdrp<0.01 , ]
LAIQTLresNEANsig <- LAIQTLresNEAN[abs(LAIQTLresNEAN$R2)>0.7 , ]
dim(LAIQTLresNEANsig)


# Writing N of archaic alleles to a file
#write.table(N_DENI, "N_DENI_methylQTL_FDRp010_sites.tsv")
#write.table(N_NEAN, "N_NEAN_methylQTL_FDRp010_sites.tsv")

#N_DENI <- read.table("N_DENI_methylQTL_FDRp010_sites.tsv")
#colnames(N_DENI) <- gsub("\\.", "-", colnames(N_DENI))
#N_NEAN <- read.table("N_NEAN_methylQTL_FDRp010_sites.tsv")
#colnames(N_NEAN) <- gsub("\\.", "-", colnames(N_NEAN))


# Writing results to a file
write.table(LAIQTLresDENI, LAIQTLresDENIpath, sep="\t")
write.table(LAIQTLresNEAN, LAIQTLresNEANpath, sep="\t")

LAIQTLresDENI <- read.table(LAIQTLresDENIpath, sep="\t")
LAIQTLresNEAN <- read.table(LAIQTLresNEANpath, sep="\t")

# ======================================
# Plotting
# ======================================

# Individuals to include
inds <- colnames(QTLALTalleles)

# QTL-LAI association results
LAIQTLres <- read.table(LAIQTLresNEANpath)
LAIQTLres$var <- rownames(LAIQTLres)

# Merging with QTL results
LAIQTLres_QTLtoolsRes <- merge(LAIQTLres, QTLtoolsResSig, by="var")
colnames(LAIQTLres_QTLtoolsRes) <- gsub(".y", "_qtl", colnames(LAIQTLres_QTLtoolsRes))
colnames(LAIQTLres_QTLtoolsRes) <- gsub(".x", "_laiqtl", colnames(LAIQTLres_QTLtoolsRes))

# Selecting variants to plot: a large QTL effect size and significant correlation
# with LAI
LAIQTLres_sig <- LAIQTLres_QTLtoolsRes[abs(LAIQTLres_QTLtoolsRes$R2)>0.7 , ]
dim(LAIQTLres_sig)
sorted <- LAIQTLres_sig[order(-LAIQTLres_sig$R2),]
head(sorted)
# chr14_23411054

# QTL, but no LAI correlation
sigQTL_noLAI <- LAIQTLres_QTLtoolsRes[abs(LAIQTLres_QTLtoolsRes$R2)<0.001 , ]
dim(sigQTL_noLAI)
sorted <- sigQTL_noLAI[order(-sigQTL_noLAI$R2),]
tail(sorted)
# chr7_73505188 

sorted[sorted$var=="chr3_118968730",]
# cg05408284
sorted[sorted$target=="cg05408284",]

# Plotting top 10 eQTLs
#plotvars <- rownames(head(sorted, n=10))
plotvars <- c("chr3_118968730")

# Getting targets of variants to plot
# Subsetting these from the M-val data on command line to molPhenoData
# QTLtoolsResSig[QTLtoolsResSig$var %in% plotvars , ]$molecular_trait_id
# Or all tested positions
#QTLtoolsResSig_bed <- QTLtoolsResSig[, colnames(QTLtoolsResSig) %in% c("chromosome", "target_start", "target_end", "molecular_trait_id")]
#QTLtoolsResSig_bed <- QTLtoolsResSig_bed[,match(c("chromosome", "target_start", "target_end", "molecular_trait_id"), colnames(QTLtoolsResSig_bed))]
#write.table(QTLtoolsResSig_bed, "methylQTL_mvals_plotvartargets.bed", quote=FALSE, row.names=FALSE, sep="\t")

molPhenoData <- read.table(molPhenoDataPath, header = TRUE)
rownames(molPhenoData) <- molPhenoData$id
molPhenoData <- molPhenoData[ , !(colnames(molPhenoData) %in% c("chr", "start", "end", "id"))]
colnames(molPhenoData) <- gsub("\\.", "-", colnames(molPhenoData))
molPhenoData <- molPhenoData[ , (colnames(molPhenoData) %in% inds)]
molPhenoData <- molPhenoData[ , match(inds, colnames(molPhenoData))]

for (var in plotvars){
  # Plotting association between QTL genotype and ancestry
  V <- as.character(QTLALTalleles[var,])
  PAP <- as.character(N_PAP[var,])
  EAS <- as.character(N_EAS[var,])
  plotData <- data.frame(V, PAP)
  plotData$island <- sapply(strsplit(inds,"-"), `[`, 1)
  plotData$island <- gsub("MPI", "Korowai", plotData$island)
  plotData$island <- gsub("MTW", "Mentawai", plotData$island)
  plotData$island <- gsub("SMB", "Sumba", plotData$island)
  plotData$island <- factor(plotData$island, levels=c("Mentawai", "Sumba", "Korowai"))
  
  # Dotplot with jitter
  p1 <- ggplot(plotData, aes(x=V, y=PAP, color=island)) +
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="N of QTL minor alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_discrete(name="N of inferred Papuan alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.title = element_blank(),
          legend.text=element_text(size=10))
  
  # Plotting genotype and target molecular trait
  # Finding target gene
  QTLtoolsResSig$var <- paste(QTLtoolsResSig$var_chromosome, QTLtoolsResSig$var_start, sep="_")
  QTLtoolsResSig$var <- paste("chr", QTLtoolsResSig$var, sep="")
  # If there are multiple targets, this only takes the first one
  target <- QTLtoolsResSig[QTLtoolsResSig$var == var , ]$molecular_trait_id[1]
  targetMolPhenoData <- molPhenoData[rownames(molPhenoData)==target , ]
  plotData$molPhenoData <- as.numeric(targetMolPhenoData)
  
  # Boxplot with jitter
  p2 <- ggplot(plotData, aes(x=V, y=molPhenoData, color=island)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="N of QTL minor alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_continuous(name=paste("Normalized", target, sep=" ")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.position = "none")
  
  p3 <- ggplot(plotData, aes(x=PAP, y=molPhenoData, color=island)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="N of inferred Papuan alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_continuous(name=paste("Normalized", target, sep=" ")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.position = "none")
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position="none")
  
  varpos <- gsub("_", ":", var)
  rsid <- as.character(QTLtoolsResSig[QTLtoolsResSig$var == var , ]$var_rsID)[1]
  
  #grid.arrange(p1, p2, p3, legend, ncol=4, widths=c(4, 4, 4, 1.5))

  ggsave(paste(outputPrefix, "_QTLfdrp001_EASPAP_plotPAP_", var, "_", target, ".pdf", sep=""), 
         grid.arrange(p1, p2, p3, legend, ncol=4, widths=c(4, 4, 4, 1.5),
                      top = textGrob(paste(varpos, rsid, sep=", "), x = 0, hjust = 0,
                                     gp = gpar(fontsize=12, font=1))),
         width = 12.5, height = 4)
}

# For archaic ancestry
# TODO: run for both N and D

# QTL-LAI association results
LAIQTLresDENI <- read.table(LAIQTLresDENIpath)
LAIQTLresDENIsig <- LAIQTLresDENI[LAIQTLresDENI$fdrp<0.01 , ]
LAIQTLresDENIsig <- LAIQTLresDENIsig[abs(LAIQTLresDENIsig$R2)>0.7 , ]
dim(LAIQTLresDENIsig)
DENIsorted <- LAIQTLresDENIsig[order(-LAIQTLresDENIsig$R2),]
head(DENIsorted, n=10)

LAIQTLresNEAN <- read.table(LAIQTLresNEANpath)
LAIQTLresNEANsig <- LAIQTLresNEAN[LAIQTLresNEAN$fdrp<0.01 , ]
LAIQTLresNEANsig <- LAIQTLresNEANsig[abs(LAIQTLresNEANsig$R2)>0.7 , ]
dim(LAIQTLresNEANsig)
NEANsorted <- LAIQTLresNEANsig[order(-LAIQTLresNEANsig$R2),]
head(NEANsorted, n=10)

# Plotting top 10 QTLs
plotvars <- rownames(head(DENIsorted, n=10))
plotvars <- c("chr19_9697958")
# Checking if R2=1 SNPs in Deni and Nean are overlapping
#LAIQTLresNEANsig_weird <- LAIQTLresNEANsig[LAIQTLresNEANsig$R2==1 , ]
#LAIQTLresDENIsig_weird <- LAIQTLresDENIsig[LAIQTLresDENIsig$R2==1 , ]
#intersect(rownames(LAIQTLresNEANsig_weird), rownames(LAIQTLresDENIsig_weird))

#plotvars <- c("chr5_118701550", "chr1_168239293")
#LAIQTLresDENI[rownames(LAIQTLresDENI)=="chr1_168239293" , ]

#plotvars <- rownames(LAIQTLresDENIsig_weird)

# Getting targets of variants to plot
# Subsetting these from the M-val data (with BEDtools) to molPhenoData
# LAIQTLresDENIsig[LAIQTLresDENIsig$var %in% plotvars , ]$molecular_trait_id
# LAIQTLresNEANsig[LAIQTLresNEANsig$var %in% plotvars , ]$molecular_trait_id

N_DENI <- get_arc_alleles("D")
N_NEAN <- get_arc_alleles("N")
#N_DENI <- read.table("N_DENI_methylQTL_FDRp010_sites.tsv")
colnames(N_DENI) <- gsub("\\.", "-", colnames(N_DENI))
#N_NEAN <- read.table("N_NEAN_methylQTL_FDRp010_sites.tsv")
colnames(N_NEAN) <- gsub("\\.", "-", colnames(N_NEAN))

head(N_DENI)
N_DENI <- N_DENI[ , colnames(N_DENI) %in% colnames(N_PAP)]
N_NEAN <- N_NEAN[ , colnames(N_NEAN) %in% colnames(N_PAP)]

# Subsetting WGS samples from the QTL ALT-allele dataframes
head(QTLALTalleles)
colnames(QTLALTalleles) <- gsub("\\.", "-", colnames(QTLALTalleles))
QTLALTalleles <- QTLALTalleles[ , colnames(QTLALTalleles) %in% colnames(N_PAP)]

for (var in plotvars){
  # Plotting association between QTL genotype and ancestry
  V <- as.character(QTLALTalleles[var,])
  # Getting N of archaic introgressed alleles for each variant
  N_ARC <- N_DENI
  # Fixing the IDs of some individuals
  #colnames(N_ARC) <- gsub("MTW0", "MTW-0", colnames(N_ARC))
  # Sorting columns to match
  #N_ARC <- N_ARC[ , colnames(N_ARC) %in% colnames(QTLALTalleles)]
  #N_ARC <- N_ARC[ , match(colnames(QTLALTalleles), colnames(N_ARC))]
  
  ARC <- as.character(N_ARC[var,])
  #NOT_ARC <- as.numeric(chartr("012","210", ARC))
  plotData <- data.frame(V, ARC)
  plotData$island <- sapply(strsplit(inds,"-"), `[`, 1)
  plotData$island <- gsub("MPI", "Korowai", plotData$island)
  plotData$island <- gsub("MTW", "Mentawai", plotData$island)
  plotData$island <- gsub("SMB", "Sumba", plotData$island)
  plotData$island <- factor(plotData$island, levels=c("Mentawai", "Sumba", "Korowai"))
  
  # Dotplot with jitter
  p1 <- ggplot(plotData, aes(x=V, y=ARC, color=island)) +
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="N of QTL minor alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_discrete(name="N of inferred Denisovan alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.title = element_blank(),
          legend.text=element_text(size=10))
  
  # Plotting genotype and target molecular trait
  # Finding target gene
  QTLtoolsResSig$var <- paste(QTLtoolsResSig$var_chromosome, QTLtoolsResSig$var_start, sep="_")
  QTLtoolsResSig$var <- paste("chr", QTLtoolsResSig$var, sep="")
  # If there are multiple targets, this only takes the first one
  target <- QTLtoolsResSig[QTLtoolsResSig$var == var , ]$molecular_trait_id[1]
  targetMolPhenoData <- molPhenoData[rownames(molPhenoData)==target , ]
  plotData$molPhenoData <- as.numeric(targetMolPhenoData)
  
  # Boxplot with jitter
  p2 <- ggplot(plotData, aes(x=V, y=molPhenoData, color=island)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="N of QTL minor alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_continuous(name=paste("Normalized", target, sep=" ")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.position = "none")
  
  p3 <- ggplot(plotData, aes(x=ARC, y=molPhenoData, color=island)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="N of inferred Denisovan alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_continuous(name=paste("Normalized", target, sep=" ")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.position = "none")
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position="none")
  
  varpos <- gsub("_", ":", var)
  rsid <- as.character(QTLtoolsResSig[QTLtoolsResSig$var == var , ]$var_rsID)[1]
  
  ggsave(paste(outputPrefix, "_QTLfdrp001_DENI_", var, "_", target, ".pdf", sep=""), 
         grid.arrange(p1, p2, p3, legend, ncol=4, widths=c(5, 5, 5, 1.5),
                      top = textGrob(paste(varpos, rsid, sep=", "), x = 0, hjust = 0,
                                     gp = gpar(fontsize=12, font=1))),
         width = 15.5, height = 5)
  
}

# For Neanderthal
#plotvars <- rownames(head(NEANsorted, n=10))
plotvars <- c("chr5_134241411")

for (var in plotvars){
  # Plotting association between QTL genotype and ancestry
  V <- as.character(QTLALTalleles[var,])
  # Getting N of archaic introgressed alleles for each variant
  N_ARC <- N_NEAN
  # Fixing the IDs of some individuals
  #colnames(N_ARC) <- gsub("MTW0", "MTW-0", colnames(N_ARC))
  # Sorting columns to match
  #N_ARC <- N_ARC[ , colnames(N_ARC) %in% colnames(QTLALTalleles)]
  #N_ARC <- N_ARC[ , match(colnames(QTLALTalleles), colnames(N_ARC))]
  
  ARC <- as.character(N_ARC[var,])
  #NOT_ARC <- as.numeric(chartr("012","210", ARC))
  plotData <- data.frame(V, ARC)
  plotData$island <- sapply(strsplit(inds,"-"), `[`, 1)
  plotData$island <- gsub("MPI", "Korowai", plotData$island)
  plotData$island <- gsub("MTW", "Mentawai", plotData$island)
  plotData$island <- gsub("SMB", "Sumba", plotData$island)
  plotData$island <- factor(plotData$island, levels=c("Mentawai", "Sumba", "Korowai"))
  
  # Dotplot with jitter
  p1 <- ggplot(plotData, aes(x=V, y=ARC, color=island)) +
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="# of QTL minor alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_discrete(name="# of inferred Neanderthal alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.title = element_blank(),
          legend.text=element_text(size=10))
  
  # Plotting genotype and target molecular trait
  # Finding target gene
  QTLtoolsResSig$var <- paste(QTLtoolsResSig$var_chromosome, QTLtoolsResSig$var_start, sep="_")
  QTLtoolsResSig$var <- paste("chr", QTLtoolsResSig$var, sep="")
  # If there are multiple targets, this only takes the first one
  target <- QTLtoolsResSig[QTLtoolsResSig$var == var , ]$molecular_trait_id[1]
  targetMolPhenoData <- molPhenoData[rownames(molPhenoData)==target , ]
  plotData$molPhenoData <- as.numeric(targetMolPhenoData)
  
  # Boxplot with jitter
  p2 <- ggplot(plotData, aes(x=V, y=molPhenoData, color=island)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="# of QTL minor alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_continuous(name=paste("Normalized", target, sep=" ")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.position = "none")
  
  p3 <- ggplot(plotData, aes(x=ARC, y=molPhenoData, color=island)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(shape=19, size=3, width = 0.25, height = 0.25) +
    scale_color_manual(values=c(MTWcol, SMBcol, MPIcol)) + 
    theme_bw() +
    scale_x_discrete(name="# of inferred Neanderthal alleles", breaks=c(0,1,2), labels=c("0","1","2")) +
    scale_y_continuous(name=paste("Normalized", target, sep=" ")) +
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12),
          legend.position = "none")
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position="none")
  
  varpos <- gsub("_", ":", var)
  rsid <- as.character(QTLtoolsResSig[QTLtoolsResSig$var == var , ]$var_rsID)[1]
  
  ggsave(paste(outputPrefix, "_QTLfdrp001_NEAN_", var, "_", target, ".pdf", sep=""), 
         grid.arrange(p1, p2, p3, legend, ncol=4, widths=c(4, 4, 4, 1.5),
                      top = textGrob(paste(varpos, rsid, sep=", "), x = 0, hjust = 0,
                                     gp = gpar(fontsize=12, font=1))),
         width = 12.5, height = 4)
  
}


# Plotting overall variance explained by local ancestry

# QTL-LAI association results
LAIeQTLresEASPAP <- read.table("eQTL_QTLfdrp001_genotype_varexplained_EASPAP.tsv")
LAIeQTLresEASPAP$anc <- "EASPAP"
LAIeQTLresEASPAP$note <- NA
LAIeQTLresDENI <- read.table("eQTL_QTLfdrp001_genotype_varexplained_DENI.tsv")
LAIeQTLresDENI$anc <- "DENI"
LAIeQTLresNEAN <- read.table("eQTL_QTLfdrp001_genotype_varexplained_NEAN.tsv")
LAIeQTLresNEAN$anc <- "NEAN"

LAIeQTLresALL <- rbind(LAIeQTLresEASPAP, LAIeQTLresDENI, LAIeQTLresNEAN)

eqtl_easpap_sig <- LAIeQTLresEASPAP[LAIeQTLresEASPAP$R2>0.7,]
dim(eqtl_easpap_sig)

eqtl_deni_sig <- LAIeQTLresDENI[LAIeQTLresDENI$R2>0.7,]
dim(eqtl_deni_sig)

eqtl_nean_sig <- LAIeQTLresNEAN[LAIeQTLresNEAN$R2>0.7,]
dim(eqtl_nean_sig)

unique_tested <- unique(c(rownames(LAIeQTLresEASPAP), rownames(LAIeQTLresDENI), rownames(LAIeQTLresNEAN)))
length(unique_tested)
unique_sig <- unique(c(rownames(eqtl_easpap_sig), rownames(eqtl_deni_sig), rownames(eqtl_nean_sig)))
length(unique_sig)/length(unique_tested)

LAImethylQTLresEASPAP <- read.table("methylQTL_QTLfdrp001_genotype_varexplained_EASPAP.tsv")
LAImethylQTLresEASPAP$anc <- "EASPAP"
LAImethylQTLresEASPAP$note <- NA
LAImethylQTLresDENI <- read.table("methylQTL_QTLfdrp001_genotype_varexplained_DENI.tsv")
LAImethylQTLresDENI$anc <- "DENI"
LAImethylQTLresNEAN <- read.table("methylQTL_QTLfdrp001_genotype_varexplained_NEAN.tsv")
LAImethylQTLresNEAN$anc <- "NEAN"

LAImethylQTLresALL <- rbind(LAImethylQTLresEASPAP, LAImethylQTLresDENI, LAImethylQTLresNEAN)
LAImethylQTLresALL$anc <- as.factor(LAImethylQTLresALL$anc)
LAImethylQTLresALL$R2 <- as.numeric(LAImethylQTLresALL$R2)

methylqtl_easpap_sig <- LAImethylQTLresEASPAP[LAImethylQTLresEASPAP$R2>0.7,]
dim(methylqtl_easpap_sig)

methylqtl_deni_sig <- LAImethylQTLresDENI[LAImethylQTLresDENI$R2>0.7,]
dim(methylqtl_deni_sig)

methylqtl_nean_sig <- LAImethylQTLresNEAN[LAImethylQTLresNEAN$R2>0.7,]
dim(methylqtl_nean_sig)

unique_tested <- unique(c(rownames(LAImethylQTLresEASPAP), rownames(LAImethylQTLresDENI), rownames(LAImethylQTLresNEAN)))
length(unique_tested)
unique_sig <- unique(c(rownames(methylqtl_easpap_sig), rownames(methylqtl_deni_sig), rownames(methylqtl_nean_sig)))
length(unique_sig)/length(unique_tested)

p_all <- ggplot(LAIeQTLresALL, aes(x=R2, fill=anc)) +
  geom_density(alpha=0.3, position = 'identity') +
  scale_fill_manual(values=c("green", pinkish, "purple")) + 
  theme_bw() +
  geom_vline(xintercept=0.7, linetype="dashed") +
  scale_x_continuous(name="Variance in eQTL genotype explained ancestry") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "right")

p_all

ggsave("eQTL_LAI_QTLfdrp001_density.pdf", 
       p_all, width = 8, height = 5)

p_all <- ggplot(LAImethylQTLresALL, aes(x=R2, fill=anc)) +
  geom_density(alpha=0.3, position = 'identity') +
  scale_fill_manual(values=c("green", pinkish, "purple")) + 
  theme_bw() +
  geom_vline(xintercept=0.7, linetype="dashed") +
  scale_x_continuous(name="Variance in methylQTL genotype explained ancestry") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "right")

p_all

ggsave("methylQTL_LAI_QTLfdrp001_density.pdf", 
       p_all, width = 8, height = 5)

# eQTLs

p1 <- ggplot(LAImethylQTLresALL, aes(x=R2, fill=anc)) +
  geom_density(alpha=0.3, position = 'identity') +
  theme_bw() +
  scale_x_continuous(name="Variance in eQTL genotype explained by modern LA") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "none")

p2 <- ggplot(LAIeQTLresDENI, aes(x=R2)) +
  geom_histogram() + 
  theme_bw() +
  scale_x_continuous(name="Variance in eQTL genotype explained by Denisovan LA") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "none")

p3 <- ggplot(LAIeQTLresNEAN, aes(x=R2)) +
  geom_histogram() + 
  theme_bw() +
  scale_x_continuous(name="Variance in eQTL genotype explained by Neanderthal LA") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "none")

ggsave("eQTL_QTLfdrp001_histograms.pdf", 
       grid.arrange(p1, p2, p3, nrow=3, heights=c(3, 3, 3)),
       width = 5, height = 9)

# methylQTLs

p1 <- ggplot(LAImethylQTLresEASPAP, aes(x=R2)) +
  geom_histogram() + 
  theme_bw() +
  scale_x_continuous(name="Variance in methylQTL genotype explained by modern LA") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "none")

p2 <- ggplot(LAImethylQTLresDENI, aes(x=R2)) +
  geom_histogram() + 
  theme_bw() +
  scale_x_continuous(name="Variance in methylQTL genotype explained by Denisovan LA") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "none")

p3 <- ggplot(LAImethylQTLresNEAN, aes(x=R2)) +
  geom_histogram() + 
  theme_bw() +
  scale_x_continuous(name="Variance in methylQTL genotype explained by Neanderthal LA") +
  scale_y_continuous(name="Frequency") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12),
        legend.position = "none")

ggsave("methylQTL_QTLfdrp001_histograms.pdf", 
       grid.arrange(p1, p2, p3, nrow=3, heights=c(3, 3, 3)),
       width = 5, height = 9)


