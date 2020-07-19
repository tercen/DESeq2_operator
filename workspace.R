library(tercen)
library(tidyverse)
library(DESeq2)

options("tercen.workflowId" = "0c98153037f24682e45289f1760072cb")
options("tercen.stepId"     = "8bf22ba4-978e-4f33-b26c-34ae548f1737")

getOption("tercen.workflowId")
getOption("tercen.stepId")

ctx = tercenCtx()

if(inherits(try(ctx$select(".x")), 'try-error')) stop("x axis is missing.")
if(inherits(try(ctx$select(".y")), 'try-error')) stop("y axis is missing.")
if(inherits(try(ctx$select(".colorLevels")), 'try-error')) stop("color information is missing.")

alpha <- as.double(ctx$op.value('alpha'))
LFC_shrinkage <- as.logical(ctx$op.value('LFC_shrinkage'))
shrinkage_type <- as.character(ctx$op.value('shrinkage_type'))

all_data <- ctx  %>% 
  select(.x, .y, .ci, .ri,
         .colorLevels)

count_matrix <- all_data %>%
  select(.ri, .x, .y) %>%
  spread(key = ".x", value = ".y")

row_indexes <- count_matrix %>%
  select(.ri)

count_matrix <- count_matrix %>%
  select(-.ri)

colData <- select(all_data,
                  .x, condition = .colorLevels) %>%
  unique()

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

dds_results <- results(dds, alpha = alpha)

if(LFC_shrinkage) {
  dds_results <- lfcShrink(dds, type=shrinkage_type,
                      res = dds_results,
                      coef = "condition")
}

dds_results <- tibble(row_indexes, as_tibble(dds_results))

ctx$addNamespace(dds_results) %>%
  ctx$save()
