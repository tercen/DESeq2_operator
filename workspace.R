library(tercen)
library(tidyverse)
library(DESeq2)
library(ashr)
library(apeglm)

options("tercen.workflowId" = "581ff0ca91a059dff254ff6579007bbd")
options("tercen.stepId"     = "e984055c-47ca-4560-b54a-cbb930c24a87")

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

alpha <- 0.1
LFC_shrinkage <- FALSE
shrinkage_type <- "local"


all_data <- ctx %>% select(.ci, .ri,.colorLevels, .y)

count_matrix <- acast(all_data, .ri~.ci, value.var=".y")

colData <- filter(all_data, .ri == 0) %>%
  select(.ci, condition = .colorLevels)

colData$condition<-factor(colData$condition)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

dds_results <- results(dds, alpha = alpha)

if(LFC_shrinkage) dds_results <- lfcShrink(dds, type=shrinkage_type,
                                           res = dds_results,
                                           coef = "condition")

res_out_tmp <- tibble(filter(all_data, .ci == 0) %>% select(.ri),
                      as_tibble(dds_results) %>%
                      mutate(minus_log10_padj = -log10(padj)))
                     
res_out<-res_out_tmp %>%
  ctx$addNamespace()

res_out %>%
  ctx$save()


tim::build_test_data(res_table = res_out, ctx = ctx, test_name = "test1")
