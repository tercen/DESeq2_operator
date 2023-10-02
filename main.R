suppressPackageStartupMessages({
  library(tercen)
  library(tercenApi)
  library(dplyr)
  library(DESeq2)
  library(ashr)
  library(apeglm)
  library(reshape2)
})

ctx = tercenCtx()

if(length(ctx$yAxis) == 0) stop("y axis is missing.")
if(length(ctx$rnames) == 0) stop("row information is missing.")
if(length(ctx$cnames) == 0) stop("column information is missing.")
if(length(ctx$colors) == 0) stop("color information is missing.")

alpha <- ctx$op.value('alpha', as.double, 0.1)
LFC_shrinkage <- ctx$op.value('LFC_shrinkage', as.logical, TRUE)
shrinkage_type <- ctx$op.value('shrinkage_type', as.character, "apeglm")
size_factor_type <- ctx$op.value('size_factor_type', as.character, "ratio")
reference.index <- ctx$op.value("reference.index", as.double, 1)

all_data <- ctx$select(c(".ci", ".ri", ".y", ctx$colors))

dups <- all_data %>% select(.ci, .ri) %>% duplicated %>% any
stopifnot("Some cells contain multiple values." = !dups)

count_matrix <- acast(
  all_data,
  .ri ~ .ci,
  value.var = ".y",
  fun.aggregate = mean
)

df_tmp <- filter(all_data, .ri == 0) %>%
  tidyr::unite(col = "condition", unlist(ctx$colors)) %>%
  select(.ci, condition)

# Reorder factors
lev <- unique(df_tmp$condition)
lev_idx <- c(seq_along(lev)[reference.index], seq_along(lev)[-reference.index]) 
df_tmp$condition <- factor(
  df_tmp$condition, levels = lev[lev_idx]
)

colData <- df_tmp %>% unique()

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition
)

dds <- DESeq(dds, sfType = size_factor_type, quiet = TRUE)

effects <- resultsNames(dds)
effects <- effects[!effects %in% "Intercept"]

res <- list()
for(x in effects) {
  dds_results <- results(dds, alpha = alpha, name = x)
  
  if(LFC_shrinkage) {
    dds_results <- lfcShrink(
      dds,
      type = shrinkage_type,
      res = dds_results,
      coef = x,
      quiet = TRUE
    )
  }
  
  dds_results$comparison <- x
  dds_results$.ri <- as.integer(rownames(dds_results))
  res[[x]] <- dds_results
  
}

res_out <- do.call(rbind, res) %>%
  as_tibble() %>%
  mutate(neg_log10_padj = -log10(padj)) %>%
  ctx$addNamespace()

ctx$save(res_out)
