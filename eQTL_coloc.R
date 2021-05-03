# ==============================================================================
# Author(s) : Heini M Natri, heini.natri@gmail.com
# Date: 12/12/2020
# Description: eQTL colocalization analysis
# ==============================================================================

# ======================================
# Load libraries
# ======================================

source("/home/hnatri/utilities.R")
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(httr)
library(jsonlite)
library(tidyverse)
library(coloc)
library(biomaRt)
library(wiggleplotr)
library(GenomicRanges)
library(liftOver)
library(rtracklayer)
#library(locuscomparer)
library(data.table)
library(viridis)
library(UpSetR)

# ======================================
# Environment parameters
# ======================================

# Working directory
setwd("/scratch/hnatri/Indonesian/eQTL_coloc/")

# ======================================
# Read in data
# ======================================

indo_perm <- read.table("/scratch/hnatri/Indonesian/eQTL_methylQTL_result/Indonesian_eQTL_lifted_ALL_CHRS_perm10k.txt")
indo_perm_sig010 <- read.table("/scratch/hnatri/Indonesian/eQTL_methylQTL_result/Indonesian_eQTL_lifted_ALL_CHRS_perm10k_FDR010_significant.txt")
indo_nom <- read.table("/scratch/hnatri/Indonesian/eQTL_methylQTL_result/Indonesian_eQTL_lifted_ALL_CHRS_nominal1.txt")
colnames(indo_nom) <- c("molecular_trait_id", "chr", "start", "end", "strand", "length", "var_distance", "rsid", "var_chr", "var_start", "var_end", "pvalue", "slope", "topsnp")

# Indonesian eGenes for testing
indoeQTL_sig_genes <- indo_perm_sig010$V1

# European eGenes
gtex_perm <- read.table(gzfile(paste0("/scratch/hnatri/Indonesian/eQTL_coloc/GTEx_ge_blood.permuted.tsv.gz")), header=T)
#gtex_perm <- gtex_perm[order(gtex_perm$p_beta, decreasing=F),]
#gtex_sig_genes <- head(unique(gtex_perm$molecular_trait_id), n=3300)
gtex_perm$fdrp <- p.adjust(gtex_perm$p_beta, method="fdr")
gtex_perm_sig_genes <- gtex_perm[which(gtex_perm$fdrp<0.1),]$molecular_trait_object_id

lepik_perm <- read.table(gzfile(paste0("/scratch/hnatri/Indonesian/eQTL_coloc/Lepik_2017_ge_blood.permuted.tsv.gz")), header=T)
#lepik_nom <- lepik_nom[order(lepik_nom$pvalue, decreasing=F),]
#lepik_sig_genes <- head(unique(lepik_nom$molecular_trait_id), n=3300)
lepik_perm$fdrp <- p.adjust(lepik_perm$p_beta, method="fdr")
lepik_perm_sig_genes <- lepik_perm[which(lepik_perm$fdrp<0.1),]$molecular_trait_object_id

twinsuk_perm <- read.table(gzfile(paste0("/scratch/hnatri/Indonesian/eQTL_coloc/TwinsUK_ge_blood.permuted.tsv.gz")), header=T)
#twinsuk_nom <- twinsuk_nom[order(twinsuk_nom$pvalue, decreasing=F),]
#twinsuk_sig_genes <- head(unique(twinsuk_nom$molecular_trait_id), n=3300)
twinsuk_perm$fdrp <- p.adjust(twinsuk_perm$p_beta, method="fdr")
twinsuk_perm_sig_genes <- twinsuk_perm[which(twinsuk_perm$fdrp<0.1),]$molecular_trait_object_id

length(intersect(gtex_perm_sig_genes, intersect(lepik_perm_sig_genes, twinsuk_perm_sig_genes)))

# Venn diagram of the overlap
venn.plot <- VennDiagram:::venn.diagram(
  x = list("Indonesia" = indoeQTL_sig_genes, "GTEx" = gtex_perm_sig_genes, "Lepik" = lepik_perm_sig_genes, "TwinsUK" = twinsuk_perm_sig_genes),
  filename = "top_eGenes_venn.tiff",
  scaled = TRUE,
  col = "transparent",
  fill = c("orchid2", "palegreen2", "pink2", "skyblue2"),
  main.pos = c(0.5, 0.5, 0.5, 0.5),
  cex = 1.5,
  cat.cex = 1.5,
  main.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-15,15,-75,75),
  cat.dist = c(0.05,0.05,0.05,0.05),
  cat.fontfamily = "sans",
  main = "",
  fontfamily = "sans",
  na = "remove",
  inverted = FALSE)

venn.plot

# Upset plot of the overlap
listInput <- list("Indonesia" = indoeQTL_sig_genes, "GTEx" = gtex_perm_sig_genes, "Lepik" = lepik_perm_sig_genes, "TwinsUK" = twinsuk_perm_sig_genes)

upset(fromList(listInput), order.by = "freq", set_size.show=FALSE)
upsetplot <- upset(fromList(listInput), order.by = "freq", set_size.show=FALSE)

# Saving as a pdf
w=5
h=3
pdf(sprintf("/scratch/hnatri/Indonesian/eQTL_coloc/top_eGenes_upset.pdf"), width=w, height=h)
upsetplot
dev.off()

# MAF info for Indonesia
varInfo <- read.table("/scratch/hnatri/Indonesian/MAF_varID_removemissing.tsv", header=F)
colnames(varInfo) <- c("chr", "pos", "maf", "variant_ID")

infoCols <- read.table("/scratch/hnatri/Indonesian/info_cols.tsv", header=T)
infoCols$variant_ID <- paste(infoCols$CHROM, infoCols$POS, sep="_")

varInfoCols <- merge(varInfo, infoCols, by="variant_ID")
head(varInfoCols)
dim(varInfoCols)

# From ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv
#eur_eqtls <- c("GTEx_ge_blood", "Lepik_2017_ge_blood", "TwinsUK_ge_blood")
# Sample sizes
dataset_ns <- list("indonesia"=115, "GTEx_ge_blood"=670, "Lepik_2017_ge_blood"=491, "TwinsUK_ge_blood"=433)
# All pairwise comparisons
pairwise_comparisons_list <- list("indonesia_gtex"=c("indonesia", "GTEx_ge_blood"),
                                  "indonesia_lepik"=c("indonesia", "Lepik_2017_ge_blood"),
                                  "indonesia_twinsuk"=c("indonesia", "TwinsUK_ge_blood"),
                                  "gtex_lepik"=c("GTEx_ge_blood", "Lepik_2017_ge_blood"),
                                  "gtex_twinsuk"=c("GTEx_ge_blood", "TwinsUK_ge_blood"),
                                  "lepik_twinsuk"=c("Lepik_2017_ge_blood", "TwinsUK_ge_blood"))

# ======================================
# coloc functions
# ======================================

# The sensitivity analysis from coloc v4
prior.adjust <- function(summ,newp12,p1=1e-4,p2=1e-4,p12=1e-6) {
  if(is.list(summ) && "summary" %in% names(summ))
    summ <- summ$summary
  if(!identical(names(summ), c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")))
    stop("not a coloc summary vector")
  ## back calculate likelihoods
  f <- function(p12)
    prior.snp2hyp(summ["nsnps"],p12=p12,p1=p1,p2=p2)
  pr1 <- f(newp12)
  pr0 <- matrix(f(p12),nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE)
  ## make everything conformable
  ## if(is.matrix(summ) && nrow(summ)==1) summ <- as.vector(summ)
  ## if(nrow(pr1)==1) pr1 <- as.vector(pr1)
  ## if(nrow(pr1)>1) pr1 <- t(pr1)
  newpp <- matrix(summ[-1],nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE) * pr1/pr0 # prop to, not equal to
  newpp/rowSums(newpp)
}


prior.snp2hyp <- function(nsnp,p12=1e-6,p1=1e-4,p2=1e-4) {
  if(any(p12<p1*p2) || any(p12 > p1) || any(p12 > p2))
    return(NULL)
  tmp <- cbind(nsnp * p1,
               nsnp * p2,
               nsnp * (nsnp-1) * p1 * p2,
               nsnp * p12)
  tmp <- cbind(1-rowSums(tmp),tmp)
  ## if(nrow(tmp)==1) {
  ##     tmp <- c(tmp)
  ##     names(tmp) <- paste0("H",0:4)
  ## } else 
  colnames(tmp) <- paste0("H",0:4)
  tmp
}

manh.plot <- function(df,wh,
                      position=if("position" %in% names(df)) {
                        df$position
                      } else {
                        1:nrow(df)
                      }) {
  znm <- if(wh==1) { "z.df1" } else {"z.df2" }
  ## print(znm)
  ## print(head(df))
  p <- pnorm(abs(df[[znm]]),lower.tail=FALSE)*2
  ## mycol <- ifelse(A$snp %in% nCV, "red","black")
  Pal <- colorRampPalette(c('white','blue'))
  
  ##This adds a column of color values
  ## based on the y values
  Col <- Pal(100)[ceiling(100*df$SNP.PP.H4)]
  plot(position,-log10(p),col="gray20",
       bg = Col, # Fill colour
       pch = 21, # Shape: circles that can filed
       frame.plot = FALSE, # Remove the frame 
       xlab=if("position" %in% names(df)) {
         "Chromosome position"
       } else {
         "SNP number"
       },
       ylab="-log10(p)",
       xaxt='n')
  ## main=paste("Trait",wh))
  axis(side=1,labels=FALSE) 
}


sensitivity <- function(obj, rule="", gene, npoints=100, doplot=TRUE,
                        plot.manhattans=TRUE, preserve.par=FALSE, row=1) {
  stopifnot("coloc_abf" %in% class(obj))
  stopifnot("priors" %in% names(obj))
  stopifnot("summary" %in% names(obj))
  if(rule=="")
    stop("please supply a rule to define colocalisation, eg 'H4 > thr' where thr is some probability of H4 that you accept as colocalisation")
  rule.init <- rule
  rule <- gsub("(H.)","PP.\\1.abf",rule,perl=TRUE)
  
  ## multiple signals?
  if(is.data.table(obj$summary)) {
    if(!(row %in% 1:nrow(obj$summary)))
      stop("row must be between 1 and ",nrow(obj$summary))
    pp <- unlist(c(obj$summary[row,grep("PP|nsnp",names(obj$summary)),with=FALSE]))
    if(paste0("SNP.PP.H4.row",row) %in% names(obj$results)) {
      obj$results[["SNP.PP.H4"]]  <- obj$results[[paste0("SNP.PP.H4.row",row)]]
      obj$results[["z.df1"]]  <- obj$results[[paste0("z.df1.row",row)]]
      obj$results[["z.df2"]]  <- obj$results[[paste0("z.df2.row",row)]]
    } else {
      pp <- unlist(c(obj$summary[row,grep("PP|nsnp",names(obj$summary)),with=FALSE]))
    }
  } else {
    pp <- obj$summary
  }
  
  p12 <- obj$priors["p12"]
  p1 <- obj$priors["p1"]
  p2 <- obj$priors["p2"]
  check <- function(pp) { with(as.list(pp),eval(parse(text=rule))) }
  pass.init <- check(pp)
  
  passrule <- "FAIL"
  if(pass.init){passrule<-"PASS"}
  
  message("Results ",if(check(pp)) { "pass" } else { "fail" }, " decision rule ",rule.init)
  
  testp12 <- 10^seq(log10(p1*p2),log10(min(p1,p1)),length.out=npoints)
  testH <- prior.snp2hyp(pp["nsnps"],p12=testp12,p1=p1,p2=p2)
  testpp <- as.data.frame(prior.adjust(summ=pp,newp12=testp12,p1=p1,p2=p2,p12=p12))
  colnames(testpp) <- gsub("(H.)","PP.\\1.abf",colnames(testpp),perl=TRUE)
  pass <- check(testpp)
  w <- which(pass)
  testpp_pass <- testpp[w,]
  
  
  #width=10
  #height=8
  #filename <- paste("./Plots/", gene, "_colocWithLevik_sensitivity.pdf", sep="")
  #pdf(sprintf(filename), width=width, height=height)
  
  #print(
  #  if(doplot) {
  #    H <- as.character(0:4)
  #    palette(c("#ffffffff",viridis(5,alpha=1)[-1]))
  #    op <- par('mfcol', 'mar', 'mfrow','mar','mgp','las','tck')
  #    on.exit(par(op))
  #    if(!preserve.par) {
  #      if(plot.manhattans)
  #        layout(mat = matrix(1:4,2,2),
  #               heights = c(1, 1), # Heights of the two rows
  #               widths = c(2, 3)) # Widths of the two columns
  #      else
  #        par(mfcol=c(1,2))
  #    }
  #    par(mar = c(3, 3, 2, 1) # Dist' from plot to side of page
  #        ,mgp = c(2, 0.4, 0) # Dist' plot to label
  #        ,las = 1 # Rotate y-axis text
  #        ,tck = -.01 # Reduce tick length
  #    )
  #    if(plot.manhattans) {
  #      manh.plot(obj$results,1)
  #      manh.plot(obj$results,2)
  #    }
  #    m <- list(testH,as.matrix(testpp))
  #    ti <- list("Prior probabilities", "Posterior probabilities")
  #    for(i in 1:2) {
  #      ym <- if(i==1) { max(m[[i]][,-1]) } else { max(m[[i]]) }
  #      matplot(testp12,m[[i]],log="x",xlab="p12",ylab="Prob",
  #              type="b",
  #              bg = 1:5, # Fill colour
  #              pch = 21, # Shape: circles that can filed
  #              col="gray20",
  #              frame.plot = FALSE, # Remove the frame 
  #              panel.first = abline(h = seq(0, 1, 0.2), col = "grey80"),
  #              ylim=c(0,ym))
  #      title(main=ti[[i]],adj=0)
  #      title(sub=paste("shaded region:",rule.init),adj=0)
  #      ## title(main=paste("Acceptance rule (shaded region):",rule.init))
  #      ## legend("topleft",pch=rep(21,5),pt.bg=1:5,legend=paste0("H",0:4))
  #      if(i==1)
  #        legend("left",inset=c(0.1,0),bg="white",pch=rep(21,5),pt.bg=1:5,pt.cex=2,legend=paste0("H",0:4))
  #      abline(v=p12,lty="dashed",col="gray")
  #      text(p12,0.5,"results",srt=90,col="gray40")
  #      if(any(pass))
  #        rect(xleft=testp12[min(w)],ybottom=0,
  #             xright=testp12[max(w)],ytop=1,
  #             col=rgb(0,1,0,alpha=0.1), border="green") 
  #      ## add text showing rule
  #      ## mtext(paste("shaded region:",rule.init),side=3,adj=1)
  #    }
  #  }
  #)
  #dev.off()
  
  if (nrow(testpp_pass)==0){
    p12lower <- NA
    p12upper <- NA
  } else {
    p12lower <- testp12[min(w)]
    p12upper <- testp12[max(w)]
  }
  
  #invisible(cbind(testpp,p12=testp12,pass=pass))
  return(c(passrule, p12lower, p12upper))
}

# ======================================
# Perform colocalisation analysis
# ======================================

gene_coloc <- function(dataset_pair){
  message(dataset_pair)
  # eQTL summary statistics
  dataset1 <- pairwise_comparisons_list[[dataset_pair]][1]
  dataset2 <- pairwise_comparisons_list[[dataset_pair]][2]
  if (dataset1=="indonesia"){
    nom1 <- indo_nom
  } else {
    nom1 <- read.table(gzfile(paste0("/scratch/hnatri/Indonesian/eQTL_coloc/", dataset1, ".all.tsv.gz")), header=T)
  }
  nom2 <- read.table(gzfile(paste0("/scratch/hnatri/Indonesian/eQTL_coloc/", dataset2, ".all.tsv.gz")), header=T)
  
  # Running the colocalization analysis for each gene
  # An empty dataframe for results
  results <- data.frame(matrix(NA, ncol=10, nrow=length(indoeQTL_sig_genes)))
  rownames(results) <- indoeQTL_sig_genes
  colnames(results) <- c("PP0", "PP1", "PP2", "PP3", "PP4", "NsnpsColoc", "sensitivityRule", "sensitivityp12lower", "sensitivityp12upper", "note")
  
  # Counting how many genes have been processed:
  counter <- 0
  for (gene in rownames(results)[1:length(rownames(results))]) {
    counter <- counter+1
    message(paste("Processing gene", counter, "out of", length(indoeQTL_sig_genes), sep=" "))
    
    # Getting the nominal stats for the target gene
    nom1_egene <- nom1[nom1$molecular_trait_id==gene,]
    
    nom2_egene <- nom2[nom2$molecular_trait_id==gene,]
    nom2_egene$rsid <- unlist(nom2_egene$rsid)
    
    # If the df is empty, skipping this gene:
    if (nrow(nom1_egene)==0 | nrow(nom2_egene)==0) {
      results[gene, "note"] <- "Gene not in db"
      message("Gene not in db")
      next
    }
    
    # Adding the Indonesian MAF
    if (dataset1=="indonesia"){
      nom1_egene$variant_ID <- paste(nom1_egene$chr, nom1_egene$var_start, sep = "_")
      testVarInfo <- varInfo[varInfo$variant_ID %in% nom1_egene$variant_ID , ]
      nom1_egene <- nom1_egene[nom1_egene$variant_ID %in% testVarInfo$variant_ID , ]
      testVarInfo <- testVarInfo[match(nom1_egene$variant_ID, testVarInfo$variant_ID) , ]
      
      nom1_egene$maf <- testVarInfo$maf
      nom1_egene$rsid <- as.character(nom1_egene$rsid)
    } else {
      nom1_egene$rsid <- unlist(nom1_egene$rsid)
    }
    
    # Merging input data for colocalization
    input <- merge(nom1_egene, nom2_egene, by="rsid", all = F,
                   suffixes = c("_dataset1","_dataset2"))
    
    # If the merged df is empty (no shared rsIDs), skipping this gene:
    if (nrow(input)==0) {
      results[gene, "note"] <- "No shared variants"
      message("No shared variants")
      next
    }
    
    # Perform colocalisation analysis
    result <- coloc.abf(
      dataset1 = list(pvalues = as.numeric(input$pvalue_dataset1), type = "quant", N = dataset_ns[[dataset1]][1], MAF = as.numeric(input$maf_dataset1)),
      dataset2 = list(pvalues = as.numeric(input$pvalue_dataset2), type = "quant", N = dataset_ns[[dataset2]][1], MAF = as.numeric(input$maf_dataset2)))
    
    # Sensitivity analysis
    rule <- "H4 > 0.8 & H4/H3 >5"
    sensitivityres <- sensitivity(obj=result, rule=rule, gene=gene)
    
    # Adding sensitivity test result to the result df
    results[gene, "sensitivityRule"] <- rule
    results[gene, "sensitivityPassrule"] <- sensitivityres[1]
    results[gene, "sensitivityp12lower"] <- sensitivityres[2]
    results[gene, "sensitivityp12upper"] <- sensitivityres[3]
    
    # Placing results to the dataframe
    results[gene, "PP0"] <- result$summary[[2]]
    results[gene, "PP1"] <- result$summary[[3]]
    results[gene, "PP2"] <- result$summary[[4]]
    results[gene, "PP3"] <- result$summary[[5]]
    results[gene, "PP4"] <- result$summary[[6]]
    results[gene, "NsnpsColoc"] <- result$summary[[1]]
    
  }
  
  # Calculating PP4/PP3
  results$PP4_PP3 <- results$PP4/results$PP3
  
  results$note <- as.character(results$note)
  results$note[is.na(results$note)] <- "Tested"
  
  # Saving to a file
  write.table(results, paste0("/scratch/hnatri/Indonesian/eQTL_coloc/indonesia_permsig010_genes_", dataset1, "_", dataset2, "_res.tsv"), sep="\t")
  
  results
}

# For every pairwise comparison, run colocalization analysis for every gene
result_list <- lapply(names(pairwise_comparisons_list), gene_coloc)
names(result_list) <- names(pairwise_comparisons_list)

# Save the result object to a file
saveRDS(result_list, file = "/scratch/hnatri/Indonesian/eQTL_coloc/Indonesia_permsig010_genes_pairwise_comparisons_result_list")
quit(save="no")
# Restore the object
result_list <- readRDS(file = "/scratch/hnatri/Indonesian/eQTL_coloc/Indonesia_permsig010_genes_pairwise_comparisons_result_list")

# ======================================
# Finding significant colocalizations
# ======================================

#head(result_list[["GTEx_ge_blood"]])
#head(result_list[["Lepik_2017_ge_blood"]])
#head(result_list[["TwinsUK_ge_blood"]])

# Finding significant colocalizations for each pairwise comparisons
# TODO: generalize
coloc_sig <- function(dataset_pair){
  message(dataset_pair)
  
  results <- result_list[[dataset_pair]]
  
  results$note <- as.character(results$note)
  results$note[is.na(results$note)] <- "Tested"
  
  # How many genes were tested
  message("Not in db ", nrow(results[which(results$note=="Gene not in db") , ]))
  message("No shared vars ", nrow(results[which(results$note=="No shared variants") , ]))
  message("Tested ", nrow(results[which(results$note=="Tested") , ]))
  
  results_tested <- results[results$note=="Tested" , ]
  dim(results_tested)
  
  # How many of these fail/pass coloc rule
  message("FAIL ", nrow(results_tested[results_tested$sensitivityPassrule=="FAIL" , ] ))
  message("PASS ", nrow(results_tested[results_tested$sensitivityPassrule=="PASS" , ] ))
  
  results_pass <- results_tested[results_tested$sensitivityPassrule=="PASS" , ]
  dim(results_tested)
  dim(results_pass)
  
  results_pass_relaxed <- results_tested[which(results_tested$PP4>0.5 & results_tested$PP4_PP3>2) , ]
  results_fail_relaxed <- results_tested[!(results_tested$PP4>0.5 | results_tested$PP4_PP3>2) , ]
  message("Pass relaxed ", nrow(results_pass_relaxed))
  message("Fail relaxed ", nrow(results_fail_relaxed))
  
  # How many PASS genes have lower p12 of <1e−06?
  results_pass$sensitivityp12lower <- as.numeric(results_pass$sensitivityp12lower)
  results_pass_robust <- results_pass[results_pass$sensitivityp12lower< 0.000001,]
  message("Pass robust ", nrow(results_pass_robust))
  # Writing results to a file
  #write.table(results_sig, "Indonesian_eQTL_FDRp020_Lepik_p12_075_P4_sig.tsv", sep="\t")
  #write.table(results, "Indonesian_eQTL_FDRp020_Lepik_p12_075_P4.tsv", sep="\t")
  
  # Plotting min p12s
  results_pass$ypos <- "ypos"
  summary(results_pass$sensitivityp12lower)
  
  results_pass_robust
}

result_sig_list <- lapply(names(result_list), coloc_sig)
names(result_sig_list) <- names(result_list)

# Colocalized (pass robust) for all three European datasets?
#GTEx_ge_blood_pass_robust <- result_sig_list[["GTEx_ge_blood"]][result_sig_list[["GTEx_ge_blood"]]$sensitivityp12lower< 0.000001,]
#Lepik_2017_ge_blood_pass_robust <- result_sig_list[["Lepik_2017_ge_blood"]][result_sig_list[["Lepik_2017_ge_blood"]]$sensitivityp12lower< 0.000001,]
#TwinsUK_ge_blood_pass_robust <- result_sig_list[["TwinsUK_ge_blood"]][result_sig_list[["TwinsUK_ge_blood"]]$sensitivityp12lower< 0.000001,]
#
#coloc_robust_pass_intersect <- intersect(intersect(rownames(GTEx_ge_blood_pass_robust), rownames(Lepik_2017_ge_blood_pass_robust)), rownames(TwinsUK_ge_blood_pass_robust))
#length(coloc_robust_pass_intersect)

# Tested in all three European datasets, not colocalized even with a relaxed threshold in any?
coloc_fail_relaxed <- function(dataset_pair){
  message(dataset_pair)
  
  results <- result_list[[dataset_pair]]
  
  results$note <- as.character(results$note)
  results$note[is.na(results$note)] <- "Tested"
  
  results_tested <- results[results$note=="Tested" , ]
  
  results_pass_relaxed <- results_tested[which(results_tested$PP4>0.5 & results_tested$PP4_PP3>2) , ]
  results_fail_relaxed <- results_tested[!(results_tested$PP4>0.5 | results_tested$PP4_PP3>2) , ]
  message("Pass relaxed ", nrow(results_pass_relaxed))
  message("Fail relaxed ", nrow(results_fail_relaxed))
  
  results_fail_relaxed
}

coloc_fail_relaxed_list <- lapply(names(result_list), coloc_fail_relaxed)
names(coloc_fail_relaxed_list) <- names(result_list)

head(coloc_fail_relaxed_list[["GTEx_ge_blood"]])

#coloc_fail_relaxed_intersect <- intersect(intersect(rownames(coloc_fail_relaxed_list[["GTEx_ge_blood"]]), rownames(coloc_fail_relaxed_list[["Lepik_2017_ge_blood"]])), rownames(coloc_fail_relaxed_list[["TwinsUK_ge_blood"]]))
#length(coloc_fail_relaxed_intersect)

# Indonesia specific eQTLs from mashr
indo_specific_mashr <- readLines("/scratch/hnatri/Indonesian/mashr/indo_specific.tsv")
sapply(strsplit(indo_specific_mashr,"_"), `[`, 1)
indo_specific_mashr <- data.frame("gene"=sapply(strsplit(indo_specific_mashr,"_"), `[`, 1), "var"=sapply(strsplit(indo_specific_mashr,"_"), `[`, 2))
indo_specific_mashr$gene_var <- paste0(indo_specific_mashr$var, "_", indo_specific_mashr$gene)
indo_perm_sig010$gene_var <- paste0(indo_perm_sig010$V8, "_", indo_perm_sig010$V1)
length(unique(indo_specific_mashr$gene))
length(intersect(unique(indo_specific_mashr$gene), indoeQTL_sig_genes))
length(intersect(indo_perm_sig010$gene_var, indo_specific_mashr$gene_var))

# Indonesia-specific eQTLs from mashr overlapping genes that show no support for
# colocalization even with a relaxed threshold
intersect(indo_specific_mashr$gene, coloc_fail_relaxed_intersect)
length(intersect(indo_specific_mashr$gene, coloc_fail_relaxed_intersect))

length(unique(coloc_fail_relaxed_intersect))
length(unique(indo_specific_mashr$gene))

# Colocalized (pass robust) for all three European datasets?
gtex_lepik_pass_robust <- eur_pairwise_result_sig_list[["gtex_lepik"]][eur_pairwise_result_sig_list[["gtex_lepik"]]$sensitivityp12lower< 0.000001,]
gtex_twinsuk_pass_robust <- eur_pairwise_result_sig_list[["gtex_twinsuk"]][eur_pairwise_result_sig_list[["gtex_twinsuk"]]$sensitivityp12lower< 0.000001,]
lepik_twinsuk_pass_robust <- eur_pairwise_result_sig_list[["lepik_twinsuk"]][eur_pairwise_result_sig_list[["lepik_twinsuk"]]$sensitivityp12lower< 0.000001,]

# Genes with shared causal variant between all European datasets
eur_pairwise_coloc_robust_pass_intersect <- intersect(intersect(rownames(gtex_lepik_pass_robust), rownames(gtex_twinsuk_pass_robust)), rownames(lepik_twinsuk_pass_robust))
length(eur_pairwise_coloc_robust_pass_intersect)

# Genes that share a causal variant between all European datasets, but have
# a different causal variant in Indonesia
length(intersect(eur_pairwise_coloc_robust_pass_intersect, intersect(indo_specific_mashr$gene, coloc_fail_relaxed_intersect)))

# Genes shared between all European pops and Indonesia
length(intersect(eur_pairwise_coloc_robust_pass_intersect, coloc_robust_pass_intersect))

# GO over-representation analysis
library(clusterProfiler)
#gene.df <- bitr(intersect(indo_specific_mashr$gene, coloc_fail_relaxed_intersect),
#                fromType = "ENSEMBL",
#                toType = c("ENTREZID", "SYMBOL"),
#                OrgDb = org.Hs.eg.db)
#allgenes.df <- bitr(indoeQTL_sig_genes,
#                    fromType = "ENSEMBL",
#                    toType = c("ENTREZID", "SYMBOL"),
#                    OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = intersect(indo_specific_mashr$gene, coloc_fail_relaxed_intersect),
                universe      = indoeQTL_sig_genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

