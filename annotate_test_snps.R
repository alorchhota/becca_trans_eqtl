library(ioutil)
library(genomicsutil)
library(argparser)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-snps", help="target snps file", default="data/max_trans_eqtl_input.txt")
args <- add_argument(args, "-annot", help="snp annotation file (txt)", default="data/variant_annot_MAF_05_not_in_repeat.txt.gz")
args <- add_argument(args, "-th", help="minimum threshold for snp mappability", default=1)
args <- add_argument(args, "-oann", help="output file with snp annotations", default="results/becca_snps_annotated.txt")
args <- add_argument(args, "-o", help="output prefix", default="results/becca_snps")

argv = parse_args(args)

snps_fn = argv$snps
snp_annot_fn = argv$annot
min_snp_mappability = argv$th
out_pfx = argv$o
annotated_snps_fn = argv$oann

### read target snps
snps_df = read_df(snps_fn, header = T, row.names = F)
snps_df$chr = extend_chr(snps_df$chr)
print(sprintf("#unique snps: %s", length(unique(snps_df$pos))))  # 366

### read snp annotations
snp_annot_df = read_df(snp_annot_fn, sep = '\t', header = T, row.names = F)

### get target snp annotations
target_snps_annot = merge(snps_df, snp_annot_df, by = c('chr', 'pos'))
print(sprintf("#unique snps with MAF>=0.05 and not-in-repeat: %s", length(unique(target_snps_annot$pos))))  # 164

### save
out_annot = unique(target_snps_annot[, colnames(target_snps_annot)!="tissue", drop = F])
write_df(out_annot, file = annotated_snps_fn, row.names = F, col.names = T)

tissues = unique(target_snps_annot$tissue)
for(tis in tissues){
  filetered_snps = target_snps_annot$snp[(target_snps_annot$tissue == tis) & (target_snps_annot$mappability >= min_snp_mappability)]
  filtered_snps_fn = sprintf("%s_ids_%s.txt", out_pfx, tis)
  write_df(filetered_snps, file = filtered_snps_fn, row.names = F, col.names = F)
}

