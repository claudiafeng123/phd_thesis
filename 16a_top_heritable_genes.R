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
pica_target_meta <- fread(fnm)

## loading this makes things faster for plotting later
most_recent <- unlist(lapply(strsplit(list.files(file.path(OutFolder, "15a_get_line_perms"), 
                                                 pattern = "target_line_mapping"), "_"), "[[", 1))

## load lfcs
pica_lfcs <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/6f_calc_lfcs_transcriptome_wide_by_gene_per_line/2022-10-05_all_lines_with_adj_pval.tsv.gz")
heritability_d <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/2024-01-10_significant_pairs.tsv.gz")

heritability_df <- fread("/lustre/scratch123/hgi/projects/crispri_scrnaseq/outs/Pica/15c_calc_var_expl_by_target/15c_05_calc_var_expl_by_target_check_pvals/2024-01-10_heritability_with_pval_adj.tsv.gz")


## ---- Check P-values

#heritability_df <- fread(heritability_fnm_out)
rmarkdown::render(file.path(CodeFolder, "Magpie", "pipeline", paste0(subsection_name, ".Rmd")),
                  output_file = file.path(plotsdir, paste0(date, ".html")),
                  params = list(
                    date = date
                  ))



## ---- SessionInfo

sessionInfo()

