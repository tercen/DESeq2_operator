library(tercen)
library(tidyverse)
library(DESeq2)
library(ashr)
library(apeglm)

options("tercen.workflowId" = "0c98153037f24682e45289f1760072cb")
options("tercen.stepId"     = "6c765e2e-417a-4223-9122-05624c8d6e89")

getOption("tercen.workflowId")
getOption("tercen.stepId")

ctx = tercenCtx()

if(length(ctx$yAxis) == 0) stop("y axis is missing.")
if(length(ctx$rnames) == 0) stop("row information is missing.")
if(length(ctx$cnames) == 0) stop("column information is missing.")
if(length(ctx$colors) == 0) stop("color information is missing.")

alpha <- as.double(ctx$op.value('alpha'))
LFC_shrinkage <- as.logical(ctx$op.value('LFC_shrinkage'))
shrinkage_type <- as.character(ctx$op.value('shrinkage_type'))

all_data <- ctx %>% select(.ci, .ri,.colorLevels)

count_matrix <- ctx$as.matrix()

colData <- filter(all_data, .ri == 0) %>%
  select(.ci, condition = .colorLevels)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

dds_results <- results(dds, alpha = alpha)

if(LFC_shrinkage) dds_results <- lfcShrink(dds, type=shrinkage_type,
                                           res = dds_results,
                                           coef = "condition")

dds_results <- tibble(filter(all_data, .ci == 0) %>% select(.ri),
                      as_tibble(dds_results) %>%
                        mutate(minus_log10_padj = -log10(padj)))

ctx$addNamespace(dds_results) %>%
  ctx$save()
