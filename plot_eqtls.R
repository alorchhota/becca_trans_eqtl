library(ioutil)
library(miscutil)
library(genomicsutil)
library(ggplot2)
library(argparser)


args <- arg_parser("program");
args <- add_argument(args, "-eqtl", help="eqtl file", default="/work-zfs/abattle4/ashis/progres/becca_trans_eqtl/per_tissue_trans_eqtl_combined/combined_egene_fdr_1.txt")
args <- add_argument(args, "-color", help="gtex color definition file", default="data/gtex_colors.txt")
args <- add_argument(args, "-geno_pfx", help="genotype file name prefix with path (before chromosome in the name)", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_")
args <- add_argument(args, "-geno_sfx", help="genotype file name suffix (after chromosome in the name)", default="_dosage_MAF_05_not_in_repeat.RData")
args <- add_argument(args, "-geno_na", help="character for NA in the genotype file", default="-")
args <- add_argument(args, "-expr_pfx", help="expression file name prefix with path (before chromosome in the name)", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_expression/")
args <- add_argument(args, "-expr_sfx", help="expression file name suffix (after chromosome in the name)", default=".v8.normalized_expression.bed")
args <- add_argument(args, "-cov_pfx", help="covariate file name prefix with path (before chromosome in the name)", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/covariates/")
args <- add_argument(args, "-cov_sfx", help="covariate file name suffix (after chromosome in the name)", default=".v8.covariates.txt")
args <- add_argument(args, "-stat_pfx", help="eqtl stats file name prefix with path (before chromosome in the name)", default="")
args <- add_argument(args, "-stat_sfx", help="eqtl stats file name suffix (after chromosome in the name)", default="")
args <- add_argument(args, "-o", help="plot file (*.pdf)", default="results/eqtl_plot.pdf")

argv = parse_args(args)
significant_eqtl_fn = argv$eqtl
gtex_colors_fn = argv$color
geno_pfx = argv$geno_pfx
geno_sfx = argv$geno_sfx
geno_na = argv$geno_na
expr_pfx = argv$expr_pfx
expr_sfx = argv$expr_sfx
cov_pfx = argv$cov_pfx
cov_sfx = argv$cov_sfx
# eqtl_stats_pfx = "/work-zfs/abattle4/ashis/progres/covid19_eqtl/covid_related_trans_eqtl/"
# eqtl_stats_sfx = "_covid19_related_crossmap_filtered_trans_eqtls_p_0.001.txt.gz"
plt_fn = argv$o


### function to read genotype file
read_genotype_file <- function(geno_fn, geno_na = '-'){
  stopifnot(endsWith(geno_fn, ".RData"))
  load(geno_fn)
  colnames(genotype_mat_not_in_repeat) = gsub('\\.', '-', colnames(genotype_mat_not_in_repeat))
  genotype_mat_not_in_repeat[genotype_mat_not_in_repeat==geno_na] = NA
  for(cn in colnames(genotype_mat_not_in_repeat))
    genotype_mat_not_in_repeat[,cn] = as.numeric(genotype_mat_not_in_repeat[,cn])
  return(genotype_mat_not_in_repeat)
}

### gtex color
gtex_col = read.table(gtex_colors_fn, sep='\t', header = T, stringsAsFactors = F)
gtex_col$tissue_color_hex = paste0("#", gtex_col$tissue_color_hex)
rownames(gtex_col) = gtex_col$tissue_id

### read significant eqtls (one snp from each gene-tissue pair)
significant_eqtls = read_df(significant_eqtl_fn, header = T, row.names = F)
significant_eqtls = significant_eqtls[order(significant_eqtls$FDR),,drop=F]
isdupsnp = duplicated(significant_eqtls[,c('gene', 'tissue')])
significant_eqtls = significant_eqtls[!isdupsnp, , drop = F]

### plot eqtls
pdf(plt_fn)
for(si in seq_len(nrow(significant_eqtls))){
  verbose_print(sprintf("plotting eqtl %s of %s ...", si, nrow(significant_eqtls)))
  esnp = as.character(significant_eqtls[si,'snps'])
  egene = as.character(significant_eqtls[si, 'gene'])
  egene_name = as.character(significant_eqtls[si, 'gene_name'])
  tissue = as.character(significant_eqtls[si, 'tissue'])
  
  expr_fn = sprintf("%s%s%s", expr_pfx, tissue, expr_sfx)
  expr_df = read_df(expr_fn, header = T, row.names = F)
  rownames(expr_df) = expr_df$gene_id
  expr_df = expr_df[,5:ncol(expr_df), drop = F]
  stopifnot(egene %in% rownames(expr_df))
  
  ### read covariate
  cov_fn = sprintf("%s%s%s", cov_pfx, tissue, cov_sfx)
  cov_df = read_df(cov_fn, header = T, row.names = T)
  stopifnot(length(intersect(colnames(expr_df), colnames(cov_df))) == ncol(expr_df))
  cov_df = cov_df[,colnames(expr_df), drop = F]
  n_uniq = sapply(rownames(cov_df), function(idx)length(unique(as.numeric(cov_df[idx,]))))
  if(any(n_uniq<=1)){
    const_covariates = names(n_uniq[n_uniq<=1])
    cov_df = cov_df[-which(rownames(cov_df) %in% const_covariates), , drop=F]
  }
  cov_df = t(cov_df)
  
  ### remove covariates
  corrected_expr_df = rm_cov_from_expr(expr_df = expr_df, cov_df = cov_df)
  
  ### snps
  esnp_parts = strsplit(esnp, split = "_")[[1]]
  esnp_chr = esnp_parts[1]
  esnp_pos = as.numeric(esnp_parts[2])
  geno_fn = sprintf("%s%s%s", geno_pfx, esnp_chr, geno_sfx)
  geno_df = read_genotype_file(geno_fn, geno_na = geno_na)
  stopifnot(esnp %in% rownames(geno_df))
  
  common_samples = intersect(colnames(corrected_expr_df), colnames(geno_df))
  
  # eqtl_stat_fn = sprintf("%s%s%s", eqtl_stats_pfx, tissue, eqtl_stats_sfx)
  # eqtl_stat_df = read.table(file=gzfile(eqtl_stat_fn), 
  #                           header = T, 
  #                           sep = '\t', 
  #                           quote = "", 
  #                           comment.char = "", 
  #                           stringsAsFactors = F)
  # eqtl_stat_snps_parts = strsplit(eqtl_stat_df$snps, split = "_")
  # eqtl_stat_df$snp_chr = as.character(sapply(eqtl_stat_snps_parts, function(x)x[1]))
  # eqtl_stat_df$snp_pos = as.numeric(sapply(eqtl_stat_snps_parts, function(x)x[2]))
  # # eqtl_stat_df = eqtl_stat_df[(eqtl_stat_df$snp_chr == esnp_chr) & (eqtl_stat_df$snp_pos >= esnp_pos-1e6) & (eqtl_stat_df$snp_pos <= esnp_pos+1e6), , drop = F]
  # eqtl_stat_df = eqtl_stat_df[(eqtl_stat_df$snp_chr == esnp_chr) 
  #                             & (eqtl_stat_df$snp_pos >= esnp_pos-1e6) 
  #                             & (eqtl_stat_df$snp_pos <= esnp_pos+1e6)
  #                             & (eqtl_stat_df$gene == egene), , drop = F]
  
  ### boxplot with jitter
  jitter_plt_df = data.frame(geno = factor(geno_df[esnp, common_samples]),
                             expr = as.numeric(corrected_expr_df[egene, common_samples]),
                             stringsAsFactors = F)
  
  g1 <- ggplot(jitter_plt_df, aes(x=geno, y=expr)) + 
    geom_boxplot(outlier.colour = NA)+
    #geom_violin(trim = FALSE)+
    geom_jitter(position=position_jitter(width = 0.1, height = 0), color=gtex_col[tissue, 'tissue_color_hex']) +
    theme_bw() +
    xlab(sprintf("Genotype of %s", esnp)) +
    ylab(sprintf("Corrected Expression of %s", sprintf("%s (%s)", egene, egene_name))) +
    ggtitle(sprintf("%s: p<=%.3g fdr<=%.2g, beta=%.4f", tissue, as.numeric(significant_eqtls[si,"pvalue"]), as.numeric(significant_eqtls[si,"FDR"]), as.numeric(significant_eqtls[si, "beta"]) ))
  print(g1)
}

dev.off()
