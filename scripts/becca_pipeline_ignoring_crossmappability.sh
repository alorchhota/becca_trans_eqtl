### settings
gtex_expr_dir="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_expression"
gtex_cov_dir="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/covariates"
geno_pfx="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_"
geno_chromosomes="chr14"
geno_sfx="_dosage_MAF_05_not_in_repeat.RData"
geno_na="-"
gene_annot_fn="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/gencode.v26.annotation.gene.txt"
genemap_fn="/work-zfs/abattle4/lab_data/annotation/mappability_hg38_gencode26/hg38_gene_mappability.txt"
crossmap_fn="data/no_crossmap.txt"
min_gene_mappability=0.8
min_snp_mappability=1
snp_gene_dist_for_crossmap=1e6
p_th=1
root_out_dir="/work-zfs/abattle4/ashis/progres/becca_trans_eqtl_ignoring_crosmap"

eqtl_prog_dir="/work-zfs/abattle4/ashis/prog/becca_trans_eqtl"
input_test_genes="${eqtl_prog_dir}/data/becca_target_genes.txt"
target_genes="${eqtl_prog_dir}/data/becca_target_genes_ensemblid.txt"
target_genes_annot="${eqtl_prog_dir}/data/becca_target_genes_annot.txt"
input_test_snps="${eqtl_prog_dir}/data/max_trans_eqtl_input.txt"
bg_snps_annot="${eqtl_prog_dir}/data/variant_annot_MAF_05_not_in_repeat.txt.gz"
target_snps_pfx="${eqtl_prog_dir}/data/becca_snps"
target_snps_annot="${eqtl_prog_dir}/data/becca_snps_hg38_annot.txt"
gtex_col_fn="${eqtl_prog_dir}/data/gtex_colors.txt"
statdir="$root_out_dir/per_tissue_trans_eqtl"
statpfx=""
statsfx="_crossmap_filtered_trans_eqtls_p_1.txt.gz"
combined_stat_dir="$root_out_dir/per_tissue_trans_eqtl_combined"
  

tissues=( 'Artery_Aorta' 'Brain_Cerebellar_Hemisphere' 'Brain_Cerebellum' \
          'Cells_Cultured_fibroblasts' \
          'Esophagus_Gastroesophageal_Junction' 'Esophagus_Mucosa' 'Esophagus_Muscularis' \
          'Heart_Atrial_Appendage' 'Lung' 'Muscle_Skeletal' 'Nerve_Tibial' \
          'Skin_Not_Sun_Exposed_Suprapubic' 'Skin_Sun_Exposed_Lower_leg' \
          'Spleen' 'Thyroid' 'Whole_Blood')

# tissues=( 'Muscle_Skeletal' )      


annotate_test_genes=false
annotate_test_snps=false
generate_eqtl_slurm_jobs=true
run_eqtl_slurm_jobs=true
combine_eqtls=true
plot_eqtls=true

########### annotate test genes ###########################
if [[ $annotate_test_genes == true ]]
then
  cd $eqtl_prog_dir
  Rscript annotate_test_genes.R -genes "$input_test_genes" -annot "$gene_annot_fn" -map "$genemap_fn" -th $min_gene_mappability -oann "$target_genes_annot" -oensembl "$target_genes"
fi

########### annotate test snps ###########################
if [[ $annotate_test_snps == true ]]
then
  cd $eqtl_prog_dir
  Rscript annotate_test_snps.R -snps "$input_test_snps" -annot "$bg_snps_annot" -th $min_snp_mappability -oann "$target_snps_annot" -o "$target_snps_pfx"
fi

########### generate/run eqtl slurm jobs ###########################
if [[ $generate_eqtl_slurm_jobs == true ]]
then
  cd $eqtl_prog_dir
  source "scripts/slurm_util.sh"
  ntasks=5
  time=2:0:0
  mem="20GB"
  partition="shared"
  nodes=1
  
  out_dir="$statdir"
  script_dir="$out_dir/jobs"
  if [[ ! -d "$out_dir" ]]; then mkdir "$out_dir"; fi
  if [[ ! -d "$script_dir" ]]; then mkdir "$script_dir"; fi
  
  for tis in ${tissues[@]}
  do
    echo "[$(date +'%Y/%m/%d %H:%M:%S')] calling eqtls for ${tis} tissue ..."
    expr_fn="${gtex_expr_dir}/${tis}.v8.normalized_expression.bed"
    cov_fn="${gtex_cov_dir}/${tis}.v8.covariates.txt"
    out_pfx="${out_dir}/${tis}"
    target_snps="${target_snps_pfx}_ids_${tis}.txt"
    
    # check if file already exists
    check_fn="${out_pfx}_crossmap_filtered_trans_eqtls_p_1.txt.gz"  # TODO: make independent of p_th
    if [[ -f "$check_fn" ]] ; then 
      donithing=T
      continue
    fi
    
    echo "[$(date +"%Y-%m-%d %T")] writing slurm script for ${tis} tissue ..."
    script_fn="$script_dir/trans_eqtl_${tis}.sh"
    
    job_name="${tis: 0 : 6}"
    cmd_header=$(get_slurm_header ${partition} ${script_dir} ${time} ${nodes} ${ntasks} ${mem} ${job_name} "${script_fn}.%%j.out" "${script_fn}.%%j.err")
    cmd_body="module load R/3.5.1"
    cmd_body="$cmd_body\ncd '$eqtl_prog_dir'"
    cmd_body="$cmd_body\nRscript bulk_genes_trans_eqtl.R -expr '$expr_fn' -cov '$cov_fn' -geno_pfx '$geno_pfx' -geno_chr '$geno_chromosomes' -geno_sfx '$geno_sfx' -genes '$target_genes' -snps '$target_snps' -annot '$gene_annot_fn' -crossmap '$crossmap_fn' -d $snp_gene_dist_for_crossmap -p $p_th -o '$out_pfx'"
    cmd_body="$cmd_body\necho DONE"
    cmd_body="$cmd_body\n"
    
    printf "$cmd_header\n" > "$script_fn"
    printf "$cmd_body" >> "$script_fn"
    
    if [[ $run_eqtl_slurm_jobs == true ]]
    then
      echo "[$(date +"%Y-%m-%d %T")] running slurm script for ${tis} tissue ..."
      # sbatch "$script_fn"
      sh "$script_fn" 2>&1 | tee "$script_fn.manual.$(date +'%Y_%m_%d_%H_%M_%S').out"
    fi
  done
fi

### combine trans-eqtls
if [[ $combine_eqtls == true ]]
then
  cd $eqtl_prog_dir
  fdr=1
  combined_trans_eqtl_pfx="${combined_stat_dir}/combined"
  if [[ ! -d "$combined_stat_dir" ]]; then mkdir "$combined_stat_dir"; fi
  Rscript combine_bulk_genes_trans_eqtls.R  -statdir "$statdir" \
                                            -statpfx "$statpfx" \
                                            -statsfx "$statsfx" \
                                            -fdr $fdr \
                                            -geneannot "$gene_annot_fn" \
                                            -snpannot "$bg_snps_annot" \
                                            -o "$combined_trans_eqtl_pfx"
fi


if [[ $plot_eqtls == true ]]
then
  cd $eqtl_prog_dir
  
  eqtl_fn="${combined_stat_dir}/combined_egene_fdr_1.txt"
  expr_pfx="${gtex_expr_dir}/"
  expr_sfx=".v8.normalized_expression.bed"
  cov_pfx="${gtex_cov_dir}/"
  cov_sfx=".v8.covariates.txt"
  plt_fn="${combined_stat_dir}/combined_egene_fdr_1.pdf"
  Rscript plot_eqtls.R  -eqtl "$eqtl_fn" \
                        -color "$gtex_col_fn" \
                        -geno_pfx "$geno_pfx" \
                        -geno_sfx "$geno_sfx" \
                        -expr_pfx "$expr_pfx" \
                        -expr_sfx "$expr_sfx" \
                        -cov_pfx "$cov_pfx" \
                        -cov_sfx "$cov_sfx" \
                        -o "$plt_fn"
fi

