# this script combines all trans-eqtls for covid19-related genes across all tissues.

library(ioutil)
library(argparser)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-statdir", help="eqtl stat dir", default="/work-zfs/abattle4/ashis/progres/becca_trans_eqtl/per_tissue_trans_eqtl")
args <- add_argument(args, "-statpfx", help="eqtl stat prefix", default="")
args <- add_argument(args, "-statsfx", help="eqtl stat suffix", default="_crossmap_filtered_trans_eqtls_p_1.txt.gz")
args <- add_argument(args, "-fdr", help="fdr threshold", default=1)
args <- add_argument(args, "-geneannot", help="target gene annotation file (txt)", default="/work-zfs/abattle4/ashis/prog/becca_trans_eqtl/data/becca_target_genes_annot.txt")
args <- add_argument(args, "-snpannot", help="target snp annotation file (txt)", default="/work-zfs/abattle4/ashis/prog/becca_trans_eqtl/data/variant_annot_MAF_05_not_in_repeat.txt.gz")
args <- add_argument(args, "-o", help="output file prefix for combined stats", default="results/becca_combined")

argv = parse_args(args)
eqtl_stat_dir = argv$statdir
eqtl_stat_pfx = argv$statpfx
eqtl_stat_sfx = argv$statsfx
fdr_threshold = argv$fdr
fdr_combined_eqtl_stat_pfx = argv$o
gene_annot_fn = argv$geneannot
snp_annot_fn = argv$snpannot


tissues = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", 
            "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", 
            "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", 
            "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", 
            "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", 
            "Brain_Hippocampus", "Brain_Hypothalamus", 
            "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", 
            "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", 
            "Breast_Mammary_Tissue", 
            "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", 
            "Colon_Sigmoid", "Colon_Transverse", 
            "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", 
            "Esophagus_Muscularis", "Heart_Atrial_Appendage", 
            "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", 
            "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", 
            "Ovary", "Pancreas", "Pituitary", "Prostate", 
            "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", 
            "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", 
            "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")


### read and combine significant eqtls from each tissue
all_sig_trans_eqtls = NULL
all_egene_sig_trans_eqtls = NULL
all_min_egene_trans_eqtls = NULL
for(tissue in tissues){
  fn = sprintf("%s/%s%s%s", eqtl_stat_dir, eqtl_stat_pfx, tissue, eqtl_stat_sfx)
  if(!file.exists(fn)){
    next
  }
  verbose_print(sprintf("combining %s tissue ...", tissue))
  eqtl_stat_df = read.table(file = gzfile(fn), header = T, sep = '\t', quote = "", comment.char = "", stringsAsFactors = F)
  # regular FDR significant
  sig_eqtl_stat_df = eqtl_stat_df[eqtl_stat_df$FDR <= fdr_threshold, , drop=F]
  if(nrow(sig_eqtl_stat_df) > 0 ){
    sig_eqtl_stat_df$tissue = tissue
    all_sig_trans_eqtls = rbind(all_sig_trans_eqtls, sig_eqtl_stat_df)
  }
  # gene-level FDR significant
  egene_sig_eqtl_stat_df = eqtl_stat_df[(eqtl_stat_df$egene_FDR <= fdr_threshold) 
                                        & (eqtl_stat_df$pvalue <= eqtl_stat_df$egene_min_pvalue), , drop=F]
  if(nrow(egene_sig_eqtl_stat_df) > 0 ){
    egene_sig_eqtl_stat_df$tissue = tissue
    all_egene_sig_trans_eqtls = rbind(all_egene_sig_trans_eqtls, egene_sig_eqtl_stat_df)
  }
  # min gene-level pvalue
  min_egene_pvalue_eqtl_stat_df = eqtl_stat_df[eqtl_stat_df$pvalue <= eqtl_stat_df$egene_min_pvalue, , drop=F]
  if(nrow(min_egene_pvalue_eqtl_stat_df) > 0 ){
    min_egene_pvalue_eqtl_stat_df$tissue = tissue
    all_min_egene_trans_eqtls = rbind(all_min_egene_trans_eqtls, min_egene_pvalue_eqtl_stat_df)
  }
  
}

### add gene annotations
gene_annot_df = read_df(gene_annot_fn, header = T, row.names = F)
cn = colnames(gene_annot_df)
cn[cn == "gene_id"] = "gene"
cn[cn == "chr"] = "gene_chr"
cn[cn == "start_pos"] = "gene_start_pos"
cn[cn == "end_pos"] = "gene_end_pos"
cn[cn == "strand"] = "gene_strand"
cn[cn == "mappability"] = "gene_mappability"
colnames(gene_annot_df) = cn
# gene_annot_df = gene_annot_df[,c('gene', 'gene_name', 'gene_chr', 'gene_start_pos', 'gene_end_pos', 'gene_strand', 'gene_type', 'gene_mappability')]

all_sig_trans_eqtls = merge(all_sig_trans_eqtls, gene_annot_df, by = "gene", all.x = T, all.y = F)
all_egene_sig_trans_eqtls = merge(all_egene_sig_trans_eqtls, gene_annot_df, by = "gene", all.x = T, all.y = F)
all_min_egene_trans_eqtls = merge(all_min_egene_trans_eqtls, gene_annot_df, by = "gene", all.x = T, all.y = F)


### add snp annotations
snp_annot_df = read_df(snp_annot_fn, header = T, row.names = F)
cn = colnames(snp_annot_df)
cn[cn == "snp"] = "snps"
cn[cn == "chr"] = "snp_chr"
cn[cn == "pos"] = "snp_pos"
cn[cn == "variant1"] = "snp_allele1"
cn[cn == "variant2"] = "snp_allele2"
cn[cn == "mappability"] = "snp_mappability"
colnames(snp_annot_df) = cn
# snp_annot_df = snp_annot_df[,c('snps', 'snp_chr', 'snp_pos', 'snp_allele1', 'snp_allele2', 'snp_mappability')]

fdr_snp_annot_df = snp_annot_df[snp_annot_df$snps %in% all_sig_trans_eqtls$snps,,drop=F]
all_sig_trans_eqtls = merge(all_sig_trans_eqtls, fdr_snp_annot_df, by = "snps", all.x = T, all.y = F)

egene_fdr_snp_annot_df = snp_annot_df[snp_annot_df$snps %in% all_egene_sig_trans_eqtls$snps,,drop=F]
all_egene_sig_trans_eqtls = merge(all_egene_sig_trans_eqtls, egene_fdr_snp_annot_df, by = "snps", all.x = T, all.y = F)

egene_pvalue_snp_annot_df = snp_annot_df[snp_annot_df$snps %in% all_min_egene_trans_eqtls$snps,,drop=F]
all_min_egene_trans_eqtls = merge(all_min_egene_trans_eqtls, egene_pvalue_snp_annot_df, by = "snps", all.x = T, all.y = F)


### order eqtls by FDR/egene FDR/egene pvalue
all_sig_trans_eqtls = all_sig_trans_eqtls[order(all_sig_trans_eqtls$FDR, all_sig_trans_eqtls$pvalue), ,drop = F]
all_egene_sig_trans_eqtls = all_egene_sig_trans_eqtls[order(all_egene_sig_trans_eqtls$egene_FDR, all_egene_sig_trans_eqtls$pvalue), ,drop = F]
all_min_egene_trans_eqtls = all_min_egene_trans_eqtls[order(all_min_egene_trans_eqtls$egene_min_pvalue, all_min_egene_trans_eqtls$pvalue), ,drop = F]

### save
write_df(all_sig_trans_eqtls, file = sprintf("%s_fdr_%s.txt",fdr_combined_eqtl_stat_pfx, fdr_threshold), row.names = F, col.names = T)
write_df(all_egene_sig_trans_eqtls, file = sprintf("%s_egene_fdr_%s.txt",fdr_combined_eqtl_stat_pfx, fdr_threshold), row.names = F, col.names = T)
write_df(all_min_egene_trans_eqtls, file = sprintf("%s_min_pvalue_per_gene.txt",fdr_combined_eqtl_stat_pfx), row.names = F, col.names = T)
