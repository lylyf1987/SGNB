# SGNB
This is a splicing graph based negative binomial model for detecting differential expression by using RNA-seq data


## Introduction

The splicing graph-based negative binomial model is designed for detecting differentially expressed(DE) genes under
two conditions, e.g. normal and study.

For doing such an analysis, we need three steps.

1. modify gene annotation file with modify_ann, which will generate block_ann list, gene_range list and gene_size list for later usage.

1. summarize read in sam format, which is the output of TopHat, into read type data.frame with summarize_read_single_end or summarize_read_paired_end.

1. fit model with fit_SGNB_exact, which will provide p-value and expression level for each gene.

## Example 1

Here is an example.

Suppose we want to detect DE genes for human and our working directory is './'. We need to have a gene annotation file for human, e.g. './human.gtf'. And we need to use Tophat and samtools to map and filter the RNA-seq data(in .fastq format) to get the .sam files. Assume that these .sam files are saved in two folders, e.g. './group0', './group1'.

Then we can do the analysis as shown below.

```r
#STEP 1. Load the required R libraries------------------------------------------------------------
require(SGNB) 
require(data.table)
require(Rsubread)

#STEP 2. Summarize gene info from the annotation file-------------------------------------
# Create 'block_ann.RData', 'gene_range.RData', and 'gene_size.RData' based on 
# the gene annotation file 'human.gtf'.
# Run ?modify_ann to obtain the R documentation of the modify_ann function in 
# SGNB package if needed.
modify_ann(gene_ann_path='./human.gtf', 
           line_skip=5, 
           block_ann_path='./block_ann.RData', 
           gene_range_path='./gene_range.RData', 
           gene_size_path='./gene_size.RData')
           
# STEP 3. Load summarized gene info files--------------------------------------
load('./block_ann.RData')
load('./gene_range.RData')
load('./gene_size.RData')

# STEP 4. Summarize mapped reads into count table----------------------------------------------------------
# Summarize the mapped reads for two sets of .sam files saved in folders 
# “group0” and “group1” based on two conditions.
# Run ?summarize_read_single_end to obtain the R documentation of the 
# summarize_read_single_end function in SGNB package if needed.
read_summarized <- summarize_read_single_end('./group0', 
                                             './group1', 
                                             block_ann, 
                                             gene_range, 
                                             run.parrallel = TRUE, 
                                             core.num = 2)

# STEP 5：Analyze genes based on SGNB method using fit_SGNB_exact function.
res <- fit_SGNB_exact(read_summarized, gene_size, tol = 0.001, times = 100)

# Optional STEP 6. Report adjusted p-values using FDR  -------------------------------------------
res$p_value_ad <- p.adjust(res$p_value, method = 'BH') 
```

## Example 2

In this package, we already saved block_ann.RData, gene_range.RData, and gene_size.RData as internal datasets. So for human RNA-seq analysis, we don't need to run "modify_ann" to generate them again. And there is also a real_data.rda file which is the summarized real data after running "summarize_read_single_end" function. So we can use it directly to run "fit_SGNB_exact" as below.

```r
res <- fit_SGNB_exact(real_data)
```
