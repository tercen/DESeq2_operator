library(tercen)
library(dplyr)

options("tercen.workflowId" = "0c98153037f24682e45289f1760072cb")
options("tercen.stepId"     = "8bf22ba4-978e-4f33-b26c-34ae548f1737")

getOption("tercen.workflowId")
getOption("tercen.stepId")

(ctx = tercenCtx())  %>% 
  select(.x, .y, .ci, .ri, js0.condition)

(ctx = tercenCtx())  %>% 
  select(.y, .ci, .ri) %>% 
  group_by(.ci, .ri) %>%
  summarise(median = median(.y)) %>%
  ctx$addNamespace() %>%
  ctx$save()
