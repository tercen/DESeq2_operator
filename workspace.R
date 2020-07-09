library(tercen)
library(tidyverse)
library(DESeq2)

options("tercen.workflowId" = "0c98153037f24682e45289f1760072cb")
options("tercen.stepId"     = "8bf22ba4-978e-4f33-b26c-34ae548f1737")

getOption("tercen.workflowId")
getOption("tercen.stepId")

all_data <- (ctx = tercenCtx())  %>% 
  select(.x, .y, .ci, .ri, js0.condition)

count_matrix <- all_data %>%
  select(.ri, .x, .y) %>%
  pivot_wider(names_from = ".x", values_from = ".y") %>%
  select(-.ri)

colData <- select(all_data,
                  .x, condition = js0.condition) %>%
  unique()

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

dds_results <- results(dds)
