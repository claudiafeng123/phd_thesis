## mainly exploration
## adjusted p-values can be computed in the next section. this is to decide a cutoff
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggpointdensity))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(lme4))


# set relevent i/o paths
section_name <- "15c_calc_var_expl_by_target"
subsection_name <- "15c_05_calc_var_expl_by_target_check_pvals"

# set relevent i/o paths
date<- "2024-01-10"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
keep_singleton_lines <- FALSE
include_line_quality <- T
include_sex <- T
num_perms_to_beat <- 10
max_perms2compute <- 10^4
REDO <- T

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--date_new"){ date_new <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## linear model parameters
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_line_quality"){ include_line_quality <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_sex"){ include_sex <- as.logical(args[[ind + 1]]) }
  
  ## parameters
  if (arg == "--num_perms_to_beat"){ num_perms_to_beat <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_perms2compute"){ max_perms2compute <- as.numeric(args[[ind + 1]]) }
  
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))



setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)
sexes <- LineMetadata$sex; names(sexes) <- LineMetadata$name

outdir <- file.path(file.path(OutFolder, section_name, subsection_name,'/'))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_name,'/') 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- LoadData

fnm <- sort(grep(list.files(file.path(OutFolder, "7b_target_summary"), pattern = "_target_meta_data.csv", full.names = T), pattern = "params", invert = T, value = T), decreasing =  T)[1]
target_meta <- fread(fnm)

## loading this makes things faster for plotting later
most_recent <- unlist(lapply(strsplit(list.files(file.path(OutFolder, "15a_get_line_perms"), 
                                                 pattern = "target_line_mapping"), "_"), "[[", 1))
all_permutations <- readRDS(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = paste0(most_recent , '_line_set-1_'), full.names = T)[1])
fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = 'total_perms', full.names = T), decreasing = T)[1]
n_perms_per_line <- fread(fnm)
min_n_paired_lines <- n_perms_per_line %>%
  dplyr::filter(n_perms > max_perms2compute) %>%
  .$n_lines %>% min()
control_expression = get_control_expression()
line_ind_by_target <- fread(paste0(OutFolder, "15a_get_line_perms/", most_recent, "_target_line_mapping.tsv"))

## heritability
folder_name <- list.files(file.path(OutFolder, section_name), pattern = "15c_02")
most_recent <- sort(unlist(lapply(strsplit(list.files(file.path(OutFolder, section_name, folder_name)), "_"), "[[", 1)), decreasing = T)[1]
heritability_df <- fread(file.path(OutFolder, "15c_calc_var_expl_by_target", folder_name, paste0(most_recent, "_var_expl_combined_no_pval.tsv.gz")))

## add in p-values
## LFC
perm_status_fnms <- list.files(file.path(OutFolder, section_name, '15c_04_calc_var_expl_by_target_get_permutation_status', 'lfc'), pattern = "permutation_status.tsv", full.names = T, recursive = T)
lfc_pvals <- lapply(perm_status_fnms, FUN = fread) %>% bind_rows() %>% as.data.frame()
lfc_pvals <- left_join(lfc_pvals , heritability_df %>% 
                            dplyr::select(c("target", "downstream_gene_name", 'n_paired_lines')))
#post_kd_expr
most_recent <- sort(gsub(unlist(lapply(strsplit(
  list.files(file.path(OutFolder, section_name, '15c_04_calc_var_expl_by_target_get_permutation_status', 'post_kd_expr'), pattern = "permutation_status.tsv", recursive = T), "_"), "[[", 1)),
  pattern = '.*/', replacement = ""), decreasing = T)[1]
perm_status_fnms <- grep(list.files(file.path(OutFolder, section_name, '15c_04_calc_var_expl_by_target_get_permutation_status', 'post_kd_expr'), pattern = "permutation_status.tsv", full.names = T, recursive = T), pattern = most_recent, value = T)
post_kd_expr_pvals <- lapply(perm_status_fnms, FUN = fread) %>% bind_rows() %>% as.data.frame()
post_kd_expr_pvals  <- left_join(post_kd_expr_pvals  , heritability_df %>% 
                         dplyr::select(c("target", "downstream_gene_name", 'n_paired_lines')))

off_target <- fread(paste0(OutFolder, "/../Magpie/Preprocess_External_Datasets/off_target/off_target_pairs_magpie-", "2022-08-15", "-version.tsv.gz"), header = T) %>%
  dplyr::filter(target_gene != off_target_gene)
genes2exclude <- sort(unique(off_target$target_gene, off_target$off_target_gene))

## crossmapping
cross_mapped_pairs <- fread(paste0(ResourcesFolder, 'eQTLs/crossmapping_reads_iPSC_PairedEndData_with_gene_name.txt'))

## ---- DataPreprocessing

heritability_fnm_out <- paste0(outdir, '/', date, "_heritability_all_no_pval_adj.tsv.gz")


## add in target metadata
target_meta2add <- target_meta %>% dplyr::select(c('target' = 'gene',
                                                   'magpie_deg', 'n_downstream_excl_target',
                                                   contains('is_'),
                                                   'n_cells' = 'n_cells_for_analysis',
                                                   'n_paired_lines'))
heritability_df <- left_join(heritability_df,  target_meta2add)

## add a column for off-target
heritability_df$potential_off_target <- heritability_df$target %in% genes2exclude

## lfc and post-kd expression
heritability_df <- left_join(left_join(heritability_df,  
                                       lfc_pvals %>% dplyr::select(c("target", "downstream_gene_name",
                                                                     'lfc_dLL_donor_pval' = 'pval'))), 
                                       post_kd_expr_pvals %>% 
                                         dplyr::select(c("target", "downstream_gene_name",
                                                         'post_kd_expr_dLL_donor_pval' = 'pval')))
                             
## control expression
## add in wild-type expression
targets_by_donors_with_data <- split.data.frame(line_ind_by_target, f = line_ind_by_target$line_set_name)
wt_expr_df <- lapply(targets_by_donors_with_data, FUN = function(df){
  line_set_ind <- df$line_set_name[1]
  #print(line_set_ind)
  lines_in_set <- sort(unlist(strsplit(df$paired_lines[1], ';')))
  most_recent <- sort(unique(unlist(lapply(strsplit(list.files(file.path(OutFolder, "15b_downstream_heritability"), pattern = '.tsv'), "_"), "[[", 1))), decreasing = T)
  wt_downstream_meta <- fread(file.path(OutFolder, "15b_downstream_heritability", paste0(most_recent, '_', paste(lines_in_set, collapse = '_'), '.tsv.gz')))
  rtn <- left_join(heritability_df %>%
                     dplyr::filter(target %in% df$gene),
                   wt_downstream_meta %>% dplyr::select(c('downstream_gene_name', 
                                                          'control_expr_dLL_donor_pval' = 'pval', 'control_expr_dLL_donor_pval_adj' = 'pval_adj'))
  )
})
heritability_df <- as.data.frame(bind_rows(wt_expr_df))


## do a full list
fwrite(heritability_df, heritability_fnm_out, sep = '\t', compress = 'gzip')


## ---- Calculate Adjusted P-values
## do a shortened list

heritability_df <- fread(heritability_fnm_out)

pairs2compute <- heritability_df %>%
  dplyr::filter(!(is.na(control_expr_dLL_donor_pval)) & !(is.na(lfc_dLL_donor_pval)) & !(is.na(post_kd_expr_dLL_donor_pval)) & 
                  !(is.infinite(control_expr_dLL_donor_pval)) & !(is.infinite(lfc_dLL_donor_pval)) & !(is.infinite(post_kd_expr_dLL_donor_pval)) &
                  bulk_range_lfc > sig_abs_lfc_thresh &
                  n_paired_lines > min_n_paired_lines &
                  bulk_mean_control_expr > min_expr_thresh &
                  !(paste0(target, "_", downstream_gene_name) %in% paste0(cross_mapped_pairs$gene_1, "_", cross_mapped_pairs$gene_2)) & !(paste0(target, "_", downstream_gene_name) %in% paste0(cross_mapped_pairs$gene_2, "_", cross_mapped_pairs$gene_1)) &
                  !(potential_off_target)) 
##  
pairs2compute$lfc_dLL_donor_pval_adj <- p.adjust(pairs2compute$lfc_dLL_donor_pval, method = 'BH')
pairs2compute$post_kd_expr_dLL_donor_pval_adj <- p.adjust(pairs2compute$post_kd_expr_dLL_donor_pval, method = 'BH')
pairs2compute <- pairs2compute %>%
  dplyr::select(c("target", 'downstream_gene_name', 
                  "bulk_mean_control_expr", "bulk_control_expr_var", ## target meta
                  "bulk_mean_lfc", "bulk_lfc_var", ## bulk lfc meta
                  'control_expr_dLL_donor', 'control_expr_dLL_donor_pval', "control_expr_dLL_donor_pval_adj",## control expression
                  'lfc_dLL_donor', 'lfc_dLL_donor_pval', "lfc_dLL_donor_pval_adj", ## lfc
                  'post_kd_expr_dLL_donor', 'post_kd_expr_dLL_donor_pval', "post_kd_expr_dLL_donor_pval_adj", ## post knockdown,
                  contains('varExpl')
  ))
fwrite(pairs2compute, paste0(outdir, '/', date, "_heritability_with_pval_adj.tsv.gz"), sep = '\t', compress = 'gzip')

significant_pairs <- pairs2compute %>%
  dplyr::filter(lfc_dLL_donor_pval_adj < sig_pval_thresh)
fwrite(significant_pairs, paste0(outdir, '/', date, "_significant_pairs.tsv.gz"), sep = '\t', compress = 'gzip')


## ---- Check P-values

#heritability_df <- fread(heritability_fnm_out)
rmarkdown::render(file.path(CodeFolder, "Magpie", "pipeline", paste0(subsection_name, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, ".html")),
                  params = list(
                    date = date
                  ))



## ---- SessionInfo

sessionInfo()



## add in control expression heritability
fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = '_target_line_mapping', full.names = T), decreasing = T)[1]
lines_per_target <- fread(fnm)
fnm <- sort(list.files(file.path(OutFolder, "7b_target_summary"), pattern = '_target_meta_data', full.names = T), decreasing = T)[1]
target_meta <- fread(fnm)
## split by line set index
target_sets <- split.data.frame(lines_per_target, f = lines_per_target$line_set_name)
var_expl <- mapply(1:length(target_sets), FUN = function(line_set_ind){
  targets_in_set <- target_sets[[line_set_ind]]$gene
  n_paired_lines <- target_sets[[line_set_ind]]$n_paired_lines[1]
  
  ## control expression
  fnm <- sort(list.files(file.path(OutFolder, "15b_downstream_heritability"), pattern = paste0(paste(unlist(strsplit(target_sets[[line_set_ind]]$paired_lines[1], ';')), collapse = '_'),'.tsv.gz'), full.names = T), decreasing = T)[1]
  downstream_meta <- fread(fnm)
  df2merge_downstream_meta <- downstream_meta %>%
    mutate(n_paired_lines = n_paired_lines) %>%
    dplyr::select(c('downstream_gene_name',
                    n_paired_lines,
                    'bulk_mean_control_expr', 'bulk_control_expr_var', 
                    'control_expr_coef_sex' = 'coef_sex', 'control_expr_tval_sex' = 'tval_sex',
                    'control_expr_dLL_line' = 'lfc_dLL_line', 'control_expr_dLL_donor' = 'lfc_dLL_donor', 'control_expr_LL_m0' = 'LL_m0',
    ))
  
  df2merge_heritability_df <- var_expl %>%
    dplyr::filter(target %in% targets_in_set)
  rtn <- left_join(df2merge_heritability_df, df2merge_downstream_meta)
  
}, SIMPLIFY = F) %>% bind_rows() %>% as.data.frame()


## add in genomic coordinates
downstream_gene_entrez_ids <- select(org.Hs.eg.db, 
                                     keys = unique(var_expl$downstream_gene_name),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL") %>%
  dplyr::filter(!(is.na(ENTREZID)))
chr_meta <- as.list(org.Hs.egCHR[downstream_gene_entrez_ids$ENTREZID])
downstream_gene_entrez_ids$CHROM <- unlist(lapply(chr_meta[match(downstream_gene_entrez_ids$ENTREZID, names(chr_meta))], FUN = function(x){paste(x, collapse = "_")}))
var_expl$chrom <- downstream_gene_entrez_ids$CHROM[match(var_expl$downstream_gene_name, downstream_gene_entrez_ids$SYMBOL)]





mean_effect <- fread(file.path(OutFolder, "6b_calc_lfcs_transcriptome_wide_by_gene", paste0(date, "_combined_with_adj_pval.tsv.gz")))
lines_per_target_fnm <- sort(list.files(paste0(OutFolder, "/15a_get_line_perms/"), pattern = "_lines_per_target.tsv", full.names = T), decreasing = T)[1]
lines_per_target <- fread(lines_per_target_fnm)
target_meta_fnm <- sort(list.files(file.path(OutFolder, '7b_target_summary'), pattern = "_target_meta_data.csv", full.names = T), decreasing = T)[1]
target_meta <- fread(target_meta_fnm)
target_meta_by_line_fnm <- sort(list.files(file.path(OutFolder, '9b_target_gene_downregulation_by_line'), pattern = "_target_meta_by_line.csv", full.names = T), decreasing = T)[1]
target_meta_by_line <- fread(target_meta_by_line_fnm)

