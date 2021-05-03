# Annotating genomic regions using annotatr

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(annotatr)
library(ggplot2)

setwd("/Users/hnatri/Dropbox (ASU)/Annotatr")

# Opening data
snp_regions <- read.table("/Users/hnatri/Dropbox (ASU)/Indonesian_eQTL/result_29peer_5gt/eQTL_fdrp001_positions.bed")
colnames(snp_regions) <- c("chr", "start", "end")
snp_regions$chr <- paste("chr", snp_regions$chr, sep="")

snp_granges <- makeGRangesFromDataFrame(snp_regions,
                                         seqinfo=NULL,
                                         seqnames.field="chr",
                                         start.field="start",
                                         end.field="end",
                                         ignore.strand=TRUE,
                                         keep.extra.columns=TRUE)

# GRanges object. seqnames must start with "chr"!
#snp_regions_granges = read_regions(con = snp_regions, genome = 'hg19', format = 'bed')

# # Chromatin states determined by chromHMM (Ernst and Kellis (2012)) in hg19 are available for
# # nine cell lines (Gm12878, H1hesc, Hepg2, Hmec, Hsmm, Huvec, K562, Nhek, and Nhlf)
# # Create a named vector for the AnnotationHub accession codes with desired names
# # Gm12878 is an immortalized lymphoblastoid cell line
# h3k4me3_codes = c('Gm12878' = 'AH23256')
# # Fetch ah_codes from AnnotationHub and create annotations annotatr understands
# build_ah_annots(genome = 'hg19', ah_codes = h3k4me3_codes, annotation_class = 'H3K4me3')
# # The annotations as they appear in annotatr_cache
# ah_names = c('hg19_H3K4me3_Gm12878')

# Select annotations for intersection with regions
# Note inclusion of custom annotation, and use of shortcuts
#builtin_annotations()
annots <- c('hg19_genes_1to5kb',
             'hg19_genes_promoters',
             'hg19_genes_3UTRs',
             'hg19_genes_5UTRs',
             "hg19_basicgenes",
             'hg19_genes_exons',
             "hg19_genes_firstexons",
             'hg19_genes_intronexonboundaries',
             'hg19_genes_introns',
             'hg19_genes_intergenic',
             "hg19_genes_cds",
             "hg19_cpg_islands",
             "hg19_cpg_shores",
             "hg19_cpg_shelves",
             "hg19_cpg_inter",
             'hg19_enhancers_fantom',
             "hg19_cpgs")


# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)

# Intersect the regions we read in with the annotations
annotated = annotate_regions(
  regions = snp_granges,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# A GRanges object is returned
print(annotated)

# Coerce to a data.frame
df_annotated = data.frame(annotated)

# See the GRanges column of dm_annotaed expanded
print(head(df_annotated))
unique(df_annotated$annot.type)

# Save to file
output_name <- "Indo_eQTL_site_annotation.tsv"
write.table(df_annotated, file = output_name, sep="\t")

# Find the number of regions per annotation type
annsum = summarize_annotations(
  annotated_regions = annotated,
  quiet = TRUE)
print(annsum)

write.table(as.data.frame(annsum), file = "Indo_eQTL_site_annotation_counts.tsv", sep="\t")


# Plot
result_table <- read.table("Indo_eQTL_site_annotation_counts.tsv", sep = "\t", header=TRUE)

#Turn your 'annot.type' column into a character vector
result_table$annot.type <- as.character(result_table$annot.type)
#Then turn it back into a factor with the levels in the correct order
result_table$annot.type <- factor(result_table$annot.type, levels=unique(result_table$annot.type))

# Grouped
bars <- ggplot(result_table, aes(y=n, x=reorder(annot.type, -n))) +
  geom_bar(position="dodge", stat="identity") + 
  theme_bw() +
  #scale_fill_manual(values = c("azure3", "deepskyblue4", "firebrick")) +
  #scale_y_continuous(breaks=c(0,300,600,900,1200), limits = c(0,1200)) +
  theme(axis.text.x = element_text(color="#000000", size=10, angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(color="#000000", size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        plot.margin = unit(c(0.5,0.5,0.5,3), "cm"))
  
bars
