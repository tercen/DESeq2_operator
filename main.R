library(tercen)
library(tidyverse)
library(DESeq2)

all_data <- (ctx = tercenCtx())  %>% 
  select(.x, .y, .ci, .ri, js0.condition)

count_matrix <- all_data %>%
  select(.ri, .x, .y) %>%
  pivot_wider(names_from = ".x", values_from = ".y")

row_indexes <- count_matrix %>%
  select(.ri)

count_matrix <- count_matrix %>%
  select(-.ri)

colData <- select(all_data,
                  .x, condition = js0.condition) %>%
  unique()

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

dds_results <- results(dds)

dds_results <- tibble(row_indexes, as_tibble(dds_results))

ctx$addNamespace(dds_results) %>%
  ctx$save()
