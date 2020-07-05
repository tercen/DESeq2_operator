# DESeq2_two_conditions operator

## Description

`DESeq2_two_conditions` tests for differential gene expression in samples from two conditions using the `DESeq2` package from BioConductor (Love, et al, Genome Biology, 2014).

## Usage

| Input projection | Description                      |
| ---------------- | -------------------------------- |
| `color`          | Represents the groups to compare |
| `y-axis`         | Sequence counts                  |

| Input parameters | Description                                                                              |
| -----------------| ---------------------------------------------------------------------------------------- |
| `alpha`          | Numeric, adjusted _p value_ cutoff for independent filtering (default = 0.1)             |
| `LFC_shrinkage`  | Logical, whether the returned log fold-change values should be shrinked (default = TRUE) |
| `shrinkage_type` | "normal", "apeglm" or "ashr", Type of shrinkage estimator to use (default = "normal")    |

| Output relations | Description                                                               |
| ---------------- | ------------------------------------------------------------------------- |
| `pvalue`         | numeric, p-value calculated per gene                                      |
| `padj`           | numeric, p-value calculated per gene after adjusting for multiple testing |
| `baseMean`       | numeric, mean of normalized counts for all the samples                    |
| `log2FoldChange` | numeric, shrunken log2 fold-change between the two groups to compare      |

## Details

The operator uses the `DESeq2` package from BioConductor.

## References

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

See ["Analyzing RNA-seq data with DESeq2"](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for further information on DESeq2 by Love, et al.

