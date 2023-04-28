# ==============================================================================
# Author(s) : Heini M Natri, heini.natri@gmail.com
# Date : May 2020
# Description: Colocalization analysis for eQTLs and GWAS
# ==============================================================================

# ======================================
# Load libraries
# ======================================

library(dplyr)
library(ggplot2)
#library(readr)
#library(stringr)
#library(httr)
#library(jsonlite)
library(tidyverse)
library(coloc)
library(biomaRt)
library(wiggleplotr)
library(GenomicRanges)
library(biomaRt)
library(rtracklayer)
#library(locuscomparer)
library(data.table)
library(viridis)

# ======================================
# Environment parameters
# ======================================

# Working directory
setwd("/home/hnatri/Colocalization/")

# ======================================
# Read in data
# ======================================

indo_perm_sig010 <- read.table("/scratch/hnatri/Indonesian/eQTL_methylQTL_result/Indonesian_eQTL_lifted_ALL_CHRS_perm10k_FDR010_significant.txt")
indoeQTL_sig_genes <- indo_perm_sig010$V1


# ======================================
# Helper functions
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
    stop("please supply a rule to define colocalisation, eg 'H4 thr' where thr is some probability of H4 that you accept as colocalisation")
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
# Perform colocalisation analysis between eQTLs and GWAS traits
# ======================================

run_gene_trait_coloc <- function(eQTL_nominal, eQTL_N, trait, trait_associations) {
  # An empty dataframe for results
  results <- data.frame(matrix(NA, ncol=10, nrow=length(indoeQTL_sig_genes)))
  rownames(results) <- indoeQTL_sig_genes
  colnames(results) <- c("PP0", "PP1", "PP2", "PP3", "PP4", "NsnpsColoc", "sensitivityRule", "sensitivityp12lower", "sensitivityp12upper", "note")
  results$note <- "ok"
  
  colnames(trait_associations) <- gsub("variant_id", "rsid", colnames(trait_associations))
  colnames(trait_associations) <- gsub("p_value", "pval", colnames(trait_associations))
  
  # Running the colocalization analysis for each gene
  # Counting how many genes have been processed:
  counter <- 0
  for (gene in rownames(results)) {
    counter <- counter+1
    message(paste("Processing gene", counter, "out of", length(rownames(results)), sep=" "))
    # Getting the Indonesian nominal stats for the target gene
    nominal_egene <- eQTL_nominal[eQTL_nominal$molecular_trait_id==gene , ]
    
    if (nrow(nominal_egene)==0) {
      results[gene, "note"] <- "Gene not in the eQTL stats"
      next
    }
    
    trait_associations_sharedpos <- trait_associations[trait_associations$rsid %in% nominal_egene$rsid , ]
    
    # Merging input data for colocalization
    input <- merge(nominal_egene, trait_associations_sharedpos, by="rsid", all = F,
                   suffixes = c("_eqtl","_gwas"))
    
    # If the merged df is empty (no shared rsIDs), skipping this gene:
    if (nrow(input)==0) {
      results[gene, "note"] <- "No shared variants"
      next
    }
    
    # Perform colocalisation analysis
    result <- coloc.abf(
      dataset2 = list(pvalues = as.numeric(input$pval_eqtl), type = "quant", N = eQTL_N, MAF = as.numeric(input$maf)),
      dataset1 = list(pvalues = as.numeric(input$pval_gwas), type = "quant", N = 173480, MAF = as.numeric(input$ma_freq)))
    
    if (is.nan(result[[1]][6])){
      next
    } else {
      # Sensitivity analysis
      rule <- "H4 > 0.8 & H4/H3 > 5"
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
  }
  
  
  # Calculating PP4/PP3
  results$PP4_PP3 <- results$PP4/results$PP3
  
  # Writing results to a file
  message("Saving results")
  result_path <- paste0("/scratch/hnatri/Colocalization/Astle_et_al/eur_res/", eur, "_", trait, "_coloc_res_sensitivity.tsv")
  write.table(results, result_path, sep="\t")

  # Plotting min p12s
  #results_passthreshold <- results[results$sensitivityPassrule=="PASS" , ]
  #results_passthreshold$sensitivityp12lower <- as.numeric(results_passthreshold$sensitivityp12lower)
  #results_passthreshold$ypos <- "ypos"
  
  #summary(results_passthreshold$sensitivityp12lower)
  
  # Proportion with p12 lower bound below 1e−07
  #dim(results_passthreshold[results_passthreshold$sensitivityp12lower < 0.0000010 ,])
  
  #all_results_passthreshold["ENSG00000130592_cg19280572",]
  #summary(log10(all_results_passthreshold$sensitivityp12lower))
  
  #p <- ggplot(results_passthreshold, aes(x=log(sensitivityp12lower))) +
  #  geom_histogram(fill="black", bins = 20) +
  #  theme_bw() +
  #  scale_x_continuous(breaks = log(c(0.000000010, 0.00000010, 0.0000010, 0.000010, 0.00010)) ,
  #                     labels = c("1e−08", "1e−07", "1e−06", "1e−05", "1e−04")) +
  #  xlab("p12 lower bound") +
  #  ylab("Frequency") +
  #  theme(axis.text=element_text(size=10), 
  #        axis.title=element_text(size=12),
  #        legend.position = "none",
  #        axis.text.y=element_blank(),
  #        axis.ticks.y=element_blank())
  #p
  #
  #result_plot_path <- paste("/scratch/hnatri/Colocalization/Astle_et_al/", trait, "_p12_sensitivity_lowerranges_histogram.pdf", sep="")
  
  #ggsave(result_plot_path, p, width = 8, height = 8)
  
}

# Running for all traits
traits <- c("basophil_count_27863252-GCST004618-EFO_0005090",
            "basophil_percentage_of_granulocytes_27863252-GCST004634-EFO_0007995",
            "basophil_percentage_of_leukocytes_27863252-GCST004631-EFO_0007992",
            "eosinophil_count_27863252-GCST004606-EFO_0004842",
            "eosinophil_percentage_of_granulocytes_27863252-GCST004617-EFO_0007996",
            "erythrocyte_count_27863252-GCST004601-EFO_0004305",
            "granulocyte_count_27863252-GCST004614-EFO_0007987",
            "granulocyte_percentage_of_myeloid_white_cells_27863252-GCST004608-EFO_0007997",
            "hematocrit_27863252-GCST004604-EFO_0004348",
            "hemoglobin_measurement_27863252-GCST004615-EFO_0004509",
            "leukocyte_count_27863252-GCST004610-EFO_0004308",
            "lymphocyte_count_27863252-GCST004627-EFO_0004587",
            "lymphocyte_percentage_of_leukocytes_27863252-GCST004632-EFO_0007993",
            "mean_corpuscular_hemoglobin_27863252-GCST004630-EFO_0004527",
            "mean_corpuscular_hemoglobin_concentration_27863252-GCST004605-EFO_0004528",
            "mean_corpuscular_volume_27863252-GCST004602-EFO_0004526",
            "mean_platelet_volume_27863252-GCST004599-EFO_0004584",
            "monocyte_count_27863252-GCST004625-EFO_0005091",
            "monocyte_percentage_of_leukocytes_27863252-GCST004609-EFO_0007989",
            "myeloid_white_cell_count_27863252-GCST004626-EFO_0007988",
            "neutrophil_count_27863252-GCST004629-EFO_0004833",
            "neutrophil_percentage_of_granulocytes_27863252-GCST004623-EFO_0007994",
            "neutrophil_percentage_of_leukocytes_27863252-GCST004633-EFO_0007990",
            "osinophil_percentage_of_leukocytes_27863252-GCST004600-EFO_0007991",
            "platelet_component_distribution_width_27863252-GCST004616-EFO_0007984",
            "platelet_count_27863252-GCST004603-EFO_0004309",
            "platelet_crit_27863252-GCST004607-EFO_0007985",
            "red_blood_cell_distribution_width_27863252-GCST004621-EFO_0005192",
            "reticulocyte_count_27863252-GCST004611-EFO_0007986",
            "reticulocyte_count_27863252-GCST004612-EFO_0007986",
            "reticulocyte_count_27863252-GCST004619-EFO_0007986",
            "reticulocyte_count_27863252-GCST004622-EFO_0007986",
            "reticulocyte_count_27863252-GCST004628-EFO_0007986",
            "sum_of_basophil_and_neutrophil_counts_27863252-GCST004620-EFO_0004833",
            "sum_of_eosinophil_and_basophil_counts_27863252-GCST004624-EFO_0005090",
            "sum_of_neutrophil_and_eosinophil_counts_27863252-GCST004613-EFO_0004833")

# From ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv
eur_eqtls <- c("GTEx_ge_blood", "Lepik_2017_ge_blood", "TwinsUK_ge_blood")

eQTL_Ns <- list("GTEx_ge_blood"=670, "Lepik_2017_ge_blood"=471, "TwinsUK_ge_blood"=1364)

for (eur in eur_eqtls[1]){
  message(eur)
  eQTL_nominal <- read.table(gzfile(paste0("/scratch/hnatri/Indonesian/eQTL_coloc/", eur, ".all.tsv.gz")), header=T)
  colnames(eQTL_nominal) <- gsub("pvalue", "pval", colnames(eQTL_nominal))
  eQTL_nominal <- eQTL_nominal[which(eQTL_nominal$molecular_trait_id %in% indoeQTL_sig_genes),]
  eQTL_N <- eQTL_Ns[[eur]]
  for (trait in traits){
    message(trait)
    # Astle et aL. 2016. The Allelic Landscape of Human Blood Cell Trait Variation 
    # and Links to Common Complex Disease.
    gwas_path <- paste0("/scratch/hnatri/Colocalization/Astle_et_al/", trait, "-Build37.f.tsv.gz")
    trait_associations <- read.table(gzfile(gwas_path), header=TRUE)
    #message("Finished importing the GWAS stats")
    run_gene_trait_coloc(eQTL_nominal, eQTL_N, trait, trait_associations)
  }
}
