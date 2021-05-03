# ==============================================================================
# Author(s) : Heini M Natri, heini.natri@gmail.com
# Date : May 2020
# Description: Colocalization analysis for methylQTLs and GWAS
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
#library(tidyverse)
library(coloc)
#library(biomaRt)
library(wiggleplotr)
library(GenomicRanges)
library(rtracklayer)
#library(locuscomparer)
library(data.table)
#library(viridis)

suppressPackageStartupMessages(library(optparse))

# ======================================
# Environment parameters
# ======================================

# Command line parser options

option.list <- list(
  make_option("--trait", type="character", default="", help="Trait")
  )

opt <- parse_args(OptionParser(option_list=option.list))

# ======================================
# Environment parameters
# ======================================

# Working directory
setwd("/home/hnatri/Colocalization/")

# ======================================
# Read in data
# ======================================

# Significant ecpgs
indoQTL <- read.table("/scratch/hnatri/Indonesian/eQTL_methylQTL_result/Indonesian_methylQTL_ALL_CHUNKS_perm10k_FDR001_significant.txt")
indoQTL_sig_targets <- indoQTL$V1
colnames(indoQTL) <- c("molecular_trait_id", "chromosome", "target_start",
                        "target_end", "strand", "nvars_tested", "var_distance",
                        "var_rsID", "var_chromosome", "var_start", "var_end",
                        "dfs", "dummy1", "1st_param_betadist", "2nd_param_betadist", "nominalp",
                        "slope", "adjp", "adjp_beta", "fdrp", "dunno")

length(indoQTL$molecular_trait_id)

# Nominal stats for all variants in the cis-region of each significant ecpg
indoQTL_nominal <- read.table("/scratch/hnatri/Indonesian/eQTL_methylQTL_result/combined_methylQTL_nominal_stats_permFDRp001_for_colocGWAS.tsv", header=TRUE)
colnames(indoQTL_nominal) <- c("target", "chr", "start", "end", "strand",
                                "length", "var_distance", "rsid", "var_chr",
                                "var_start", "var_end", "pval", "slope", "topsnp")
indoQTL_nominal <- as.data.table(indoQTL_nominal)
indoQTL_nominal$variant_ID <- paste(indoQTL_nominal$chr, indoQTL_nominal$var_start, sep = "_")

# MAF info
varInfo <- read.table("/scratch/hnatri/Indo_data/MAF_varID_removemissing.tsv", header=F)
colnames(varInfo) <- c("chr", "pos", "maf", "variant_ID")
varInfo <- as.data.table(varInfo)

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


sensitivity <- function(obj, rule="", cpg, npoints=100, doplot=TRUE,
                        plot.manhattans=TRUE, preserve.par=FALSE, row=1) {
  #stopifnot("coloc_abf" %in% class(obj))
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
  #filename <- paste("./Plots/", cpg, "_colocWithLevik_sensitivity.pdf", sep="")
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
# Perform colocalisation analysis between Indonesian QTLs and selected GWAS
# traits
# ======================================

#trait <- "mean_corpuscular_hemoglobin_27863252-GCST004630-EFO_0004527"
trait <- opt$trait

# Astle et al. 2016. The Allelic Landscape of Human Blood Cell Trait Variation 
# and Links to Common Complex Disease.
message(trait)
gwas_path <- paste("/scratch/hnatri/Colocalization/Astle_et_al/", trait, "-Build37.f.tsv.gz", sep="")
#trait_associations <- as.data.table(read.table(gwas_path, header=TRUE))
#head(trait_associations)
#dim(trait_associations)
trait_associations <- read.delim(gzfile(gwas_path))
trait_associations <- data.table(trait_associations)

colnames(trait_associations) <- gsub("variant_id", "rsid", colnames(trait_associations))
colnames(trait_associations) <- gsub("p_value", "pval", colnames(trait_associations))
#trait_associations$rsid <- unlist(trait_associations$rsid)
#trait_associations <- data.frame(sapply(trait_associations, function(x) unlist(x) )  )
# Running the colocalization analysis for each cpg
# Counting how many cpgs have been processed:
#counter <- 0
# head(rownames(results))

trait_res <- lapply(indoQTL$molecular_trait_id, function(cpg) {
  # A named vector for results
  #cpg <- "cg21645762"
  results <-rep(NA, 12)
  names(results) <- c("PP0", "PP1", "PP2", "PP3", "PP4", "PP4_PP3", "NsnpsColoc", "sensitivityRule", "sensitivityp12lower", "sensitivityp12upper", "note", "target")
  results["note"] <- "ok"
  results["target"] <- cpg
  
  #counter <- counter+1
  message(paste("Processing CpG", cpg, sep=" "))
  #message(paste("Processing CpG", counter, "out of", length(rownames(results)), sep=" "))
  # Getting the Indonesian nominal stats for the target cpg
  indoQTL_nominal_testCpG <- indoQTL_nominal[target==cpg]
  
  # Adding MAF
  testVarInfo <- varInfo[variant_ID %in% indoQTL_nominal_testCpG$variant_ID]
  indoQTL_nominal_testCpG <- indoQTL_nominal_testCpG[variant_ID %in% testVarInfo$variant_ID]
  testVarInfo <- testVarInfo[match(indoQTL_nominal_testCpG$variant_ID, testVarInfo$variant_ID) , ]
  
  indoQTL_nominal_testCpG$maf <- testVarInfo$maf
  indoQTL_nominal_testCpG$rsid <- as.character(indoQTL_nominal_testCpG$rsid)
  
  trait_associations_sharedpos <- trait_associations[rsid %in% indoQTL_nominal_testCpG$rsid]
  #dim(trait_associations_sharedpos)
  
  # Merging input data for colocalization
  input <- merge(indoQTL_nominal_testCpG, trait_associations_sharedpos, by="rsid", all = F,
                 suffixes = c("_indo","_gwas"))
  
  # If the merged df is empty (no shared rsIDs), skipping this cpg:
  if (nrow(input)==0) {
    results["note"] <- "No shared variants"
    results
  } else {
  
  # TODO: GWAS N?
  
  # Perform colocalisation analysis
  result <- coloc.abf(
    dataset2 = list(pvalues = as.numeric(input$pval_indo), type = "quant", N = 115, MAF = as.numeric(input$maf)),
    dataset1 = list(pvalues = as.numeric(input$pval_gwas), type = "quant", N = 10000, MAF = as.numeric(input$ma_freq)))
  
  #message(class(result))
  #names(result)
  # Sensitivity analysis
  rule <- "H4 > 0.8 & H4/H3 >5"
  sensitivityres <- sensitivity(obj=result, rule=rule, cpg=cpg)
  
  # Adding sensitivity test result to the result df
  results["sensitivityRule"] <- rule
  results["sensitivityPassrule"] <- sensitivityres[1]
  results["sensitivityp12lower"] <- as.numeric(sensitivityres[2])
  results["sensitivityp12upper"] <- as.numeric(sensitivityres[3])
  
  # Placing results to the dataframe
  results["PP0"] <- as.numeric(result$summary[[2]])
  results["PP1"] <- as.numeric(result$summary[[3]])
  results["PP2"] <- as.numeric(result$summary[[4]])
  results["PP3"] <- as.numeric(result$summary[[5]])
  results["PP4"] <- as.numeric(result$summary[[6]])
  results["NsnpsColoc"] <- as.numeric(result$summary[[1]])
  
  # Calculating PP4/PP3
  results["PP4_PP3"] <- as.numeric(results["PP4"])/as.numeric(results["PP3"])
  
  results
  }
})

#results[complete.cases(results),]
#results["ENSG00000079335_cg07146231", ]

names(trait_res) <- indoQTL$molecular_trait_id

# Collapsing to a table
trait_res <- as.data.frame(do.call(rbind, trait_res))
head(trait_res)

# Selecting significant results based on PP4/PP3>5 and PP4>0.8
#all_resultssig <- all_results[all_results$PP4_PP3>5 , ]
#all_resultssig <- all_resultssig[all_resultssig$PP4>0.8 , ]
#dim(all_resultssig)

# Writing results to a file
#"sum_of_neutrophil_and_eosinophil_counts_27863252-GCST004613-EFO_0004833"
result_path <- paste("/scratch/hnatri/Colocalization/Astle_et_al/methylQTL_", trait, "_coloc_res_sensitivity.tsv", sep="")
write.table(trait_res, result_path, sep="\t")
#results <- read.table("/scratch/hnatri/Colocalization/Astle_et_al/sum_of_neutrophil_and_eosinophil_counts_27863252-GCST004613-EFO_0004833_coloc_res_sensitivity.tsv", sep="\t")

# Plotting min p12s
#results_passthreshold <- trait_res[results$sensitivityPassrule=="PASS" , ]
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
#result_plot_path <- paste("/scratch/hnatri/Colocalization/Astle_et_al/methylQTL_", trait, "_p12_sensitivity_lowerranges_histogram.pdf", sep="")
#
#ggsave(result_plot_path, p, width = 8, height = 8)

