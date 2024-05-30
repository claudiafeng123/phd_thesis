suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(lme4))
suppressMessages(library(lme4qtl))
suppressMessages(library(lattice))
suppressMessages(library(variancePartition))
suppressMessages(library(gtools))
suppressMessages(library(lmtest))
suppressMessages(library(org.Hs.eg.db))

# set relevent i/o paths
section_name <- "15b_downstream_heritability"
subsection_long_name <- "15b_02_downstream_heritability_calc_control_heritability"

# set relevent i/o paths
date <- "2023-10-18"
HomeFolder <- "/lustre/scratch123/hgi/teams/parts/cf14/crispr_scrnaseq_hipsci/"
ExperimentName <- "Pica"

# options for testing
io_path <- paste0("scripts/io/", ExperimentName, "_io.R")
utils_path <- "scripts/io/Magpie_Utils.R"
control_tag <- c("unassigned", "NonTarget") # use all unassigned cells as control
keep_singleton_lines <- FALSE
line_set_ind <- 7
start_ind <- 4100
no_genes2test <- 100
perm_seed <- 0


args <- commandArgs(trailingOnly = TRUE)
for (ind in 1:length(args)){
  arg <- args[[ind]]
  ##needed arguments
  if (arg == "--date"){ date <- args[[ind + 1]] }
  if (arg == "--home_folder"){ HomeFolder <- args[[ind + 1]] }
  if (arg == "--no_genes2test"){ no_genes2test <- as.numeric(args[[ind + 1]]) }
  if (arg == "--start_ind"){ start_ind <- as.numeric(args[[ind + 1]]) }
  
  ## donors to include
  if (arg == "--donors2include"){ donors2include <- args[[ind + 1]] }
  if (arg == "--line_set_ind"){ line_set_ind <- as.numeric(args[[ind + 1]]) }
  
  ## not needed
  if (arg == "--section_name"){ section_name <- args[[ind + 1]] }
  if (arg == "--utils_path"){ utils_path <- args[[ind + 1]] }
  if (arg == "--io_path"){ io_path <- args[[ind + 1]] }
  
  ## permuting
  if (arg == "--permute_cells"){ permute_cells <- as.logical(args[[ind + 1]]) }
  if (arg == "--permute_lines"){ permute_lines <- as.logical(args[[ind + 1]]) }
  if (arg == "--rnd_seed_start"){ rnd_seed_start <- as.numeric(args[[ind + 1]]) }
  if (arg == "--max_seeds2test"){ max_seeds2test <- as.numeric(args[[ind + 1]]) }
  if (arg == "--n_perms"){ n_perms <- as.numeric(args[[ind + 1]]) }
  if (arg == "--keep_singleton_lines"){ keep_singleton_lines <- as.logical(args[[ind + 1]]) }
}

print(paste0(section_name))
print(paste("date:", date))
print(paste("home folder:", HomeFolder))
if(exists("donors2include")){print('including donors: '); print(donors2include)}



setwd(HomeFolder)
source(io_path)
source(utils_path)

GuideMetadata <- fread(GuideMetadataPath)
LineMetadata <- fread(LineMetadataPath)

fnm <- sort(list.files(file.path(OutFolder, "15a_get_line_perms"), pattern = "_lines_per_target.tsv", full.names = T), decreasing = T)[1]
donor_sets <- fread(fnm )
if (exists("donors2include")){
  d <- sort(unlist(strsplit(donors2include, '_')))
} else {
  d <- sort(unlist(strsplit(donor_sets$lines[which(donor_sets$line_set_name == line_set_ind)], ";")))
}

outdir <- file.path(file.path(OutFolder, section_name, subsection_long_name, paste0('line_set_ind_', line_set_ind)))
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
plotsdir <- file.path(HTMLFolder, "pipeline", section_name, subsection_long_name, paste0('line_set_ind_', line_set_ind)) 
if(!dir.exists(plotsdir)) dir.create(plotsdir, recursive = T)


## ---- load data

gene_expression <- get_control_expression() 
gene_expression$target <- NonTargetGeneName
gene_expression <- filter(gene_expression, donor %in% d)


n_cells <- dim(gene_expression)[1]
n_lines <- length(unique(gene_expression$cell_line))
n_donors <- length(unique(gene_expression$donor))
print("cells per line:")
table(gene_expression$cell_line)

## add guide

## cell metadata (for guides)
most_recent <- sort(gsub(
  list.files(file.path(OutFolder, "5a_combine"), pattern = "control_cells_metadata.tsv", full.names = T),
  pattern = paste0(".*5a_combine[/]|_.*"),
  replacement = ""),
  decreasing = T)[1]
cell_metadata_fnms <- grep(list.files(file.path(OutFolder, "5a_combine"), pattern = "control_cells_metadata.tsv", full.names = T),
                           pattern = most_recent, value = T)
cell_metadata <- lapply(cell_metadata_fnms, FUN = fread) %>%
  bind_rows() %>% as.data.frame() %>%
  #dplyr::filter(target_gene_name == target_gene) %>%
  dplyr::select(c("V1", "guide" = 'feature_call'))
## add to gene expression df
gene_expression <- left_join(gene_expression,cell_metadata, by = c("id" = "V1"))


## ---- Throw out singleton lines

if (keep_singleton_lines == F){
  print("removing singleton lines...")
  lines2keep <- data.frame(
    line = unique(gene_expression$cell_line),
    donor = unlist(lapply(strsplit(unique(gene_expression$cell_line), "_"), "[[", 1)),
    n_cells = as.numeric(table(gene_expression$cell_line))
  ) %>% filter(n_cells > min_cells_per_gene_per_donor)
  donors2keep <- table(lines2keep$donor)
  donors2keep  <- names(donors2keep)[which(donors2keep  == 2)]
  lines2keep <- lines2keep %>% filter(donor %in% donors2keep)
  print(paste0("keeping ", paste(lines2keep$line, collapse = ', ')))
  gene_expression <- filter(gene_expression, cell_line %in% lines2keep$line)
  print("cells per line:")
  print(table(gene_expression$cell_line))
} else {
  lines2keep <- unique(gene_expression$cell_line)
}
n_cells_used <- dim(gene_expression)[1]
n_lines_used <- length(unique(gene_expression$cell_line))


## ---- Calculate Heritability

### some functions
most_recent <- sort(gsub(list.files(paste0(OutFolder, "6a_run_control_lm")), pattern = "_.*", replace = ""), decreasing = T)[1]
downstream_gene_fnms <- list.files(paste0(OutFolder, "6a_run_control_lm"), pattern = most_recent )
downstream_gene_fnms <- downstream_gene_fnms[(start_ind+1):(start_ind + no_genes2test)]

fnm <- sort(list.files(file.path(OutFolder, '15a_get_line_perms'), 
                       pattern = paste0("_line_set-", line_set_ind, "_", length(lines2keep$line)/2, "_perms_rnd-seed-", perm_seed, ".RDS"), full.names = T),
                       decreasing = T)[1]
all_permutations <- readRDS(fnm)

registerDoParallel(50)
res <- foreach(downstream_gene_fnm = downstream_gene_fnms) %dopar% {
  
  downstream_gene_ind <- which(gsub(colnames(gene_expression), pattern = ":", replacement = '-') == gsub(downstream_gene_fnm, pattern = paste0(most_recent, "_gene-|[.]RDS"), replacement = ""))
  downstream_gene_name <- unlist(lapply(strsplit(colnames(gene_expression)[downstream_gene_ind], ":"), "[[", 2))
  
  
  df4lmm_control <- get_df4lmm(kd_expr = gene_expression %>% mutate(target = NonTargetGeneName),
                                 downstream_gene_name = downstream_gene_name,
                                 lines2include = lines2keep$line,
                                 include_line_quality = F,
                                 include_sex = T)
  
  
  ## with fixed effects
  models_to_test <- gsub(var_expl_models, pattern = '[+] on_target_expr', replacement = '')
  names(models_to_test) <- names(var_expl_models)
  
  df4lmm = df4lmm_control
  df2write <- run_perm(df4lmm = df4lmm_control, 
                       val2calc = 'post_kd_expr',
                       models_to_test = models_to_test)
  fwrite(df2write, paste0(outdir, '/', date, "_line-set-ind_", line_set_ind, '_', downstream_gene_name, "-", dim(df2write)[1] - 1, ".tsv.gz"), sep = '\t')
  
}

##LL_m2 <- logLik(lmer(formula =  as.formula(paste0('post_kd_expr ', models_to_test['m2'])), 
#                      data = df4lmm_control))
#LL_m1 <- logLik(lmer(formula =  as.formula(paste0('post_kd_expr ', models_to_test['m1'])), 
#                     data = df4lmm_control))
## ---- SessionInfo

sessionInfo()

## AARS had error at 840
## EXOSC2 had error at 2690

#fit1 <- lmer(lfc ~ on_target_expr + (1|cell_line) + (1|donor), df4lmm)
#fit2 <- lmer(lfc ~ on_target_expr + (1|cell_line), df4lmm)
#lrtest(fit1, fit2)

