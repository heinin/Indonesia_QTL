library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(tidyr)
library(edgeR)
library(tibble)

setwd("/Users/hnatri/Dropbox (ASU)/Indonesian_eQTL/RNAseq/")

rnaseq_data <- read.table("counts_allGENCODE.tsv", header=T)
gtf <- read.table("Homo_sapiens.GRCh38.90.gtf", header=F, skip=5, sep="\t")
gtf_gene <- gtf[gtf$V3=="gene",]

coordinates <- GRanges(seqnames = gtf_gene$V1, strand=gtf_gene$V7, ranges=IRanges(start=gtf_gene$V4, end=gtf_gene$V5), info=gtf_gene$V9)

# Lifting the start and end coordinates from GRCh38 to hg19
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch

seqlevelsStyle(coordinates) = "UCSC"  # necessary
coordinates19 = liftOver(coordinates, ch)
class(coordinates19)

ranges(coordinates19)
mcols(coordinates19)

coordinates19 <- unlist(coordinates19)

coordinates19_df <- data.frame(seqnames=seqnames(coordinates19),
                                     starts=start(coordinates19),
                                     ends=end(coordinates19),
                                     info=mcols(coordinates19)$info)

# Splitting the info column to get the gene id
coordinates19_df <- separate(data = coordinates19_df, col = info, into = c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype"), sep = "\\;")
coordinates19_df <- separate(data = coordinates19_df, col = gene_id, into = c("dummy", "gene_id_ensemble"), sep = " ")
coordinates19_df <- coordinates19_df[,c("seqnames", "starts", "ends", "gene_id_ensemble")]
colnames(coordinates19_df) <- gsub("gene_id_ensemble", "Geneid", colnames(coordinates19_df))

# Merging with the expression data
rnaseq_data_info <- merge(coordinates19_df, rnaseq_data, by="Geneid")
head(rnaseq_data_info)

# Adding QTLtools columns
colnames(rnaseq_data_info) <- gsub("Geneid", "pid", colnames(rnaseq_data_info))
colnames(rnaseq_data_info) <- gsub("ends", "end", colnames(rnaseq_data_info))
colnames(rnaseq_data_info) <- gsub("starts", "start", colnames(rnaseq_data_info))
colnames(rnaseq_data_info) <- gsub("seqnames", "Chr", colnames(rnaseq_data_info))
rnaseq_data_info$gid <- rnaseq_data_info$pid
rnaseq_data_info$strand <- "+"

infocols <- c("Chr", "start", "end", "pid", "gid", "strand")
samplecols <- setdiff(colnames(rnaseq_data_info), infocols)
rnaseq_data_info <- rnaseq_data_info[,c(infocols, samplecols)]

# Sorting
#rnaseq_data_info_sorted <- rnaseq_data_info[with(rnaseq_data_info, order(Chr, start)),]
#write.table(rnaseq_data_info_sorted, "lifted19_counts_allGENCODE_QTLtools_sorted.bed", quote=F, row.names = F, col.names = T)

# Adding some genes to the old data
old_rnaseq <- read.table("indoRNA_full_data.tsv", header=T)
old_gene_info <- read.table("gene_info.tsv", header=T)

old_gene_info_coordinates <- GRanges(seqnames = old_gene_info$Chr, ranges=IRanges(start=old_gene_info$Start, end=old_gene_info$End), info=old_gene_info$Geneid)

# Lifting the start and end coordinates from GRCh38 to hg19

seqlevelsStyle(old_gene_info_coordinates) = "UCSC"  # necessary
old_gene_info_coordinates19 = liftOver(old_gene_info_coordinates, ch)
class(old_gene_info_coordinates19)

ranges(old_gene_info_coordinates19)
mcols(old_gene_info_coordinates19)

old_gene_info_coordinates19 <- unlist(old_gene_info_coordinates19)

old_gene_info_coordinates19_df <- data.frame(seqnames=seqnames(old_gene_info_coordinates19),
                               starts=start(old_gene_info_coordinates19),
                               ends=end(old_gene_info_coordinates19),
                               info=mcols(old_gene_info_coordinates19)$info)
colnames(old_gene_info_coordinates19_df) <- gsub("info", "Geneid", colnames(old_gene_info_coordinates19_df))

# Merging with the expression data
old_rnaseq_info <- merge(old_gene_info_coordinates19_df, old_rnaseq, by="Geneid")
head(old_rnaseq_info)

# Adding QTLtools columns
colnames(old_rnaseq_info) <- gsub("Geneid", "pid", colnames(old_rnaseq_info))
colnames(old_rnaseq_info) <- gsub("ends", "end", colnames(old_rnaseq_info))
colnames(old_rnaseq_info) <- gsub("starts", "start", colnames(old_rnaseq_info))
colnames(old_rnaseq_info) <- gsub("seqnames", "Chr", colnames(old_rnaseq_info))
old_rnaseq_info$gid <- old_rnaseq_info$pid
old_rnaseq_info$strand <- "+"

infocols <- c("Chr", "start", "end", "pid", "gid", "strand")
samplecols <- setdiff(colnames(old_rnaseq_info), infocols)
old_rnaseq_info <- old_rnaseq_info[,c(infocols, samplecols)]

# Adding rows
# Genes of interest
genes_to_add <- c("ENSG00000158856", "ENSG00000206341", "ENSG00000204632")
#rnaseq_data_info[rnaseq_data_info$gid=="ENSG00000204632",]
#old_rnaseq_info[old_rnaseq_info$gid=="ENSG00000204632",]

rnaseq_data_info_selected <- rnaseq_data_info[rnaseq_data_info$gid %in% genes_to_add,]


setdiff(colnames(rnaseq_data_info_selected), colnames(old_rnaseq_info))
setdiff(colnames(old_rnaseq_info), colnames(rnaseq_data_info_selected))

#colnames(old_rnaseq_info) <- gsub("SMB.HPM016", "SMB.HPM.016", colnames(old_rnaseq_info))

old_rnaseq_info_with_extra <- rbind(old_rnaseq_info, rnaseq_data_info_selected)

# To FPKM and filtering
gene_info <- old_rnaseq_info_with_extra[,1:4]
gene_info$length <- gene_info$end-gene_info$start
counts_only <- old_rnaseq_info_with_extra[,7:length(colnames(old_rnaseq_info_with_extra))]

# Creating the DGEList object
dge <- DGEList(counts=counts_only, genes=gene_info)
colnames(dge) <- colnames(counts_only)
dge$samples$samplenames <- colnames(counts_only)

#gene_indices <- as.vector(rownames(dge$genes))

# Calculating RPKM values
genes.length <- as.vector(dge$genes$length)
rpkm <- rpkm(dge, genes.length, log=FALSE)
colnames(rpkm) <- colnames(counts_only)
rownames(rpkm) <- gene_info$pid
dim(rpkm)

# Calculating mean expression values
rpkm_mean = apply(rpkm,1,mean,na.rm=TRUE)

# New dataframe with gene names mean expression values
rpkm_mean_exp <- data.frame(matrix(nrow=nrow(rpkm)))
rpkm_mean_exp$pid <- dge$genes$pid
rpkm_mean_exp$chr <- dge$genes$chr
rpkm_mean_exp$mean_exp <- rpkm_mean

# Keeping genes that have mean FPKM of at least 0.5 and at least 6 reads in at least 50 samples
#keep <- rpkm_mean_exp$mean_exp > 0.5
#dge <- dge[keep,,keep.lib.size=FALSE]
keep <- rpkm_mean_exp$mean_exp > 0.1
dge <- dge[keep,,keep.lib.size=FALSE]
keep <- rowSums(dge$counts > 6) >= 10
dge <- dge[keep,,keep.lib.size=FALSE]
# Removing NA counts
keep <- rowSums(is.na(dge$counts) == TRUE) == 0
dge <- dge[keep,,keep.lib.size=FALSE]
dge <- calcNormFactors(dge)
keep <- rowSums(is.nan(dge$counts) == TRUE) == 0
dge <- dge[keep,,keep.lib.size=FALSE]
dge <- calcNormFactors(dge)
keep <- rowSums(is.infinite(dge$counts) == TRUE) == 0
dge <- dge[keep,,keep.lib.size=FALSE]
dge <- calcNormFactors(dge)

# To FPKM
filtered_genes <- dge$genes
filtered_genes.length <- as.vector(filtered_genes$length)
filtered_rpkm <- rpkm(dge, filtered_genes.length, log=FALSE)

filtered_gene_names <- dge$genes$pid
rownames(filtered_rpkm) <- filtered_gene_names
dim(filtered_rpkm)

#####################
### Normalizing
#####################

# Fitting samples to normal quantiles
data_qqnorm <- (apply(filtered_rpkm, 2, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
colnames(data_qqnorm) <- colnames(counts_only)

# Fitting genes to normal quantiles
data_qqnorm <- t(apply(data_qqnorm, 1, function(xx){qqnorm(rank(xx, ties.method = "random"), plot = F)$x}))
colnames(data_qqnorm) <- colnames(counts_only)
rownames(data_qqnorm) <- rownames(filtered_rpkm)

data_qqnorm_1 <- data.matrix(data_qqnorm)

# Adding QTLtools bed columns
chr_col = as.vector(dge$genes$Chr)
start_col = as.vector(dge$genes$start)
end_col = as.vector(dge$genes$end)
Geneid_col = as.vector(dge$genes$pid)

first_sample <- as.vector(colnames(counts_only))[1]

data_qqnorm_2 <- as_tibble(data_qqnorm_1)
data_qqnorm_2 <- add_column(data_qqnorm_2, chr=chr_col, .before = first_sample)
data_qqnorm_2 <- add_column(data_qqnorm_2, start=start_col, .before = first_sample)
data_qqnorm_2 <- add_column(data_qqnorm_2, end=end_col, .before = first_sample)
data_qqnorm_2 <- add_column(data_qqnorm_2, pid=Geneid_col, .before = first_sample)
data_qqnorm_2 <- add_column(data_qqnorm_2, gid=Geneid_col, .before = first_sample)
data_qqnorm_2 <- add_column(data_qqnorm_2, strand="+", .before = first_sample)

data_qqnorm_3 <- as.data.frame(data_qqnorm_2)

# Writing to file
write.table(data_qqnorm_3, file="liftedhg19_withExtraGenes_n123_fpkm_filtered_normalized_sample_gene_for_QTLtools_geneID.csv", row.names=FALSE, quote = FALSE, sep="\t")

data_qqnorm_3 <- read.table("liftedhg19_withExtraGenes_n123_fpkm_filtered_normalized_sample_gene_for_QTLtools_geneID.csv", sep="\t", header=T)

matching_names <- c("MPI-025", "MPI-048", "MPI-051", "MPI-061", "MPI-065", "MPI-242", "MPI-296", "MPI-318", "MPI-333", "MPI-335", "MPI-336", "MPI-369", "MPI-376", "MPI-379", "MPI-381", "MPI-383", "MPI-389", "MTW-MDB-001", "MTW-MDB-002", "MTW-MDB-009", "MTW-MDB-011", "MTW-MDB-014", "MTW-MDB-015", "MTW-MDB-016", "MTW-MDB-019", "MTW-MDB-020", "MTW-MDB-024", "MTW-MDB-025", "MTW-MDB-026", "MTW-MDB-028", "MTW-MDB-032", "MTW-MDB-034", "MTW-MDB-035", "MTW-TLL-007", "MTW-TLL-008", "MTW-TLL-016", "MTW-TLL-017", "MTW-TLL-019", "MTW-TLL-020", "MTW-TLL-022", "MTW-TLL-023", "MTW-TLL-024", "MTW-TLL-025", "MTW-TLL-026", "MTW-TLL-029", "MTW-TLL-034", "MTW-TLL-038", "MTW-TLL-039", "MTW-TLL-043", "SMB-ANK-004", "SMB-ANK-005", "SMB-ANK-006", "SMB-ANK-007", "SMB-ANK-009", "SMB-ANK-011", "SMB-ANK-013", "SMB-ANK-015", "SMB-ANK-019", "SMB-ANK-024", "SMB-ANK-026", "SMB-ANK-027", "SMB-ANK-028", "SMB-ANK-029", "SMB-ANK-030", "SMB-WNG-002", "SMB-WNG-003", "SMB-WNG-004", "SMB-WNG-006", "SMB-WNG-007", "SMB-WNG-009", "SMB-WNG-012", "SMB-WNG-014", "SMB-WNG-017", "SMB-WNG-018", "SMB-WNG-020", "SMB-WNG-021", "SMB-WNG-022", "SMB-WNG-023", "SMB-WNG-028", "MPI-334", "MPI-345", "MPI-378", "MTW-MDB-003", "MTW-TLL-002", "MTW-TLL-006", "MTW-TLL-010", "MTW-TLL-011", "MTW-TLL-012", "MTW-TLL-013_B2", "MTW-TLL-021", "MTW-TLL-027", "MTW-TLL-030", "MTW-TLL-032", "MTW-TLL-035", "MTW-TLL-037", "MTW-TLL-040", "MTW-TLL-041", "MTW-TLL-042", "SMB-ANK-003", "SMB-ANK-016", "SMB-ANK-027_B2", "SMB-WNG-001", "MPI-381_B3", "MTW-TLL-013_B3", "SMB-ANK-016_B3", "SMB-ANK-027_B3", "SMB-BKB-008", "SMB-HPM-006", "SMB-HPM016", "SMB-HPM-018", "SMB-HPM-021", "SMB-HPM-027", "SMB-PDT-001", "SMB-PDT-007", "SMB-PDT-036", "SMB-PTB-028", "SMB-RIN-003", "SMB-RIN-009", "SMB-RIN-014", "SMB-RIN-016", "SMB-RIN-019", "SMB-WHB-024", "SMB-WNG-021_B3")

setdiff(colnames(data_qqnorm_3), matching_names)
setdiff(matching_names, colnames(data_qqnorm_3))
colnames(data_qqnorm_3) <- gsub("SMB-HPM016", "SMB-HPM-016", colnames(data_qqnorm_3))

write.table(data_qqnorm_3, file="liftedhg19_withExtraGenes_n123_fpkm_filtered_normalized_sample_gene_for_QTLtools_geneID_colnames.bed", row.names=FALSE, quote = FALSE, sep="\t")

# ENSG00000158856
data_qqnorm_3[data_qqnorm_3$pid=="ENSG00000158856",]


