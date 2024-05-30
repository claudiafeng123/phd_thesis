suppressMessages(library(umap))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library( org.Hs.eg.db))

# set relevent i/o paths
section_name <- "15c_calc_var_expl_by_target"
subsection_name <- "15c_02_calc_var_expl_by_target_variance_decomposition"
## combine the by-target results
## produce p-values

# set relevent i/o paths
date <- "2024-05-17"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
REDO <- F

args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--target_gene"){ target_gene <- args[[ind + 1]] }
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  if (arg == "--REDO"){ REDO <- as.logical(args[[ind + 1]]) }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## permuting
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_line_quality"){ include_line_quality <- as.logical(args[[ind + 1]]) }
  if (arg == "--include_sex"){ include_sex <- as.logical(args[[ind + 1]]) }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))
if(exists("donors2include")){print('including donors: '); print(donors2include)}



setwd(HomeFolder)
source(io_path)
source(utils_path)


outdir <- file.path(file.path(OutFolder, section_name, subsection_name))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_name) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)

## ---- Load Other Data

#no_fixed_effects <- fread(paste0(outdir, '/', date, '_var_expl_no_fixed_effects_combined_no_pval.tsv.gz'))
lfc_by_line <- fread(file.path(OutFolder, "6f_calc_lfcs_transcriptome_wide_by_gene_per_line", paste0(date, "_all_lines_with_adj_pval.tsv.gz")))


## ---- true values with fixed effects

fnm_out <- paste0(outdir, '/', date, '_var_expl_combined_no_pval.tsv.gz')
if (REDO == T | !(file.exists(fnm_out))){
  
  var_expl_fnms <- list.files(file.path(OutFolder, section_name, '15c_01_calc_var_expl_by_target_run_lme'), 
                              pattern = paste0(date, "_permute-cells-FALSE_permute-lines-FALSE_include-line-quality-TRUE_include-sex-FALSE"), 
                              full.names = T, recursive = T)
  var_expl <- lapply(var_expl_fnms, fread) %>%  bind_rows() %>% as.data.frame()
  var_expl <- var_expl %>% 
    dplyr::select(c('target', 'downstream_gene_name', 
                    "bulk_mean_lfc", "bulk_lfc_var", "max_abs_lfc", "bulk_range_lfc",
                    'lfc_dLL_donor', 'lfc_dLL_line', 'lfc_LL_m0',
                    'lfc_varExpl_cell_line', 'lfc_varExpl_donor', 'lfc_varExpl_guide', 'lfc_varExpl_on_target_expr', 'lfc_varExpl_Residual',
                    'lfc_coef_on_target_expr', 'lfc_tval_on_target_expr', 
                    'post_kd_expr_dLL_donor', 'post_kd_expr_dLL_line', 'post_kd_expr_LL_m0'))
  
  fwrite(var_expl, fnm_out, sep = '\t', compress = 'gzip')
} else {
  var_expl <- fread(fnm_out)
  
}

## ---- MakePlots

## load some more data

rmarkdown::render(file.path(CodeFolder, "Magpie", "pipeline", paste0(subsection_name, ".Rmd") ),
                  output_file = paste0(plotsdir, '/', date, '.html'),
                  params = list(
                    date =date,
                    section_name = section_name,
                    subsection_name = subsection_name
                  ))

## ---- SessionInfo

sessionInfo()

