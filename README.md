# DESeq2

##### Description

The `DESeq2` operator tests for differential gene expression in samples from two conditions using the `DESeq2` package from BioConductor (Love, et al, Genome Biology, 2014).

##### Usage

| Input projection | Description                      |
| ---------------- | -------------------------------- |
| `row`            | Gene name/identifier             |
| `column`         | Sample name/identifier           |
| `color`          | Represents the groups to compare |
| `y-axis`         | Sequence counts                  |

| Input parameters | Description                                                                              |
| -----------------| ---------------------------------------------------------------------------------------- |
| `alpha`          | Numeric, adjusted _p value_ cutoff for independent filtering (default = 0.1)             |
| `LFC_shrinkage`  | Logical, whether the returned log fold-change values should be shrinked (default = TRUE) |
| `shrinkage_type` | "normal", "apeglm" or "ashr", Type of shrinkage estimator to use (default = "normal")    |
| `size_factor_type` | Method for size factor estimation: either "ratio", "poscounts", or "iterate". "ratio" uses the standard median ratio method introduced in DESeq. The size factor is the median ratio of the sample over a "pseudosample": for each gene, the geometric mean of all samples. "poscounts" and "iterate" offer alternative estimators, which can be used even when all genes contain a sample with a zero. The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts. This evolved out of use cases with Paul McMurdie's phyloseq package for metagenomic samples. The "iterate" estimator iterates between estimating the dispersion with a design of ~1, and finding a size factor vector by numerically optimizing the likelihood of the ~1 model. |
| `reference.index` | Index of the reference category to be used. |

|  Output relations  | Description                                                                |
| ------------------ | -------------------------------------------------------------------------- |
| `pvalue`           | numeric, p-value calculated per gene                                       |
| `padj`             | numeric, p-value calculated per gene after adjusting for multiple testing  |
| `baseMean`         | numeric, mean of normalized counts for all the samples                     |
| `log2FoldChange`   | numeric, shrunken log2 fold-change between the two groups to compare       |
| `minus_log10_padj` | numeric, negative log10 transformation of padj for more intuitive plotting |

##### Details

The operator uses the `DESeq2` package from BioConductor.

##### References

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

See ["Analyzing RNA-seq data with DESeq2"](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for further information on DESeq2 by Love, et al.

