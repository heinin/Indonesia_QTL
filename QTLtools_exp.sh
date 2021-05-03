#!/bin/bash
#SBATCH --job-name=QTLtools_cis # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
##SBATCH --mail-user=hnatri@tgen.org # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 96:00:00

module load QTLtools/1.0;

QTLtools cis --vcf /scratch/hnatri/Indo_data/QTL_genotype_data/merged_wgs_snparray/ALL_CHRS_WGS_SNParray_phased_imputed_n115_DR2_095_maxmissing03_maf005_maxmaf095_mac2_renamed_APformatheader_GTGPonly_fillGPfield_GP090_forQTLtools.vcf.gz --bed /scratch/hnatri/Indo_data/lifted_hg19_n123_fpkm_filtered_normalized_sample_gene_for_QTLtools_geneID_sorted_wheader.bed.gz --permute 10000 --region ${CHR} --out Indonesian_eQTL_lifted_chr${CHR}_perm10k.txt --cov /scratch/hnatri/Indo_data/exp_qtl_covariates_29peer_5gtpc.tsv
