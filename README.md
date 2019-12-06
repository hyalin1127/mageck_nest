Integrative CRISPR analysis using MAGeCK-NEST
====================================================================================
Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout (MAGeCK_NEST) is a computational tool to identify important genes from the genome-scale CRISPR-Cas9 knockout screens technology. Based on the previous generation, MAGeCK_VISPR, MAGeCK_NEST includes some new features, including:
```
* Outliers removal.
* Protein-protein interaction information integration.
* QC:
  1. Histogram of beta scores (modified log fold changes)
  2. QQ-plot of z-statistics
  3. Gene set enrichment analysis (GSEA) using positive control genes as quality control.
  4. Correlation coefficients between interacting genes.
  5. QC report
```

# Prerequisites #
The input of the MAGeCK-NEST workflow are read count files. The formats are sgRNA, gene, readcounts
Ex: 
```
sgRNA         gene     sample1_readcount     sample2_readcount...
gene1_gRNA_1  gene1    557                   421
gene1_gRNA_2  gene1    295                   128
     .          .       .
     .          .       .
     .          .       .
gene2_gRNA_1  gene2    173                   68
gene2_gRNA_2  gene2    85                    38
     .          .       .
     .          .       .
     .          .       .
```
# Installation #

```
python setup.py install
```
# Usage #

```
mageck_nest.py [-h] -k COUNT_TABLE -d DESIGN_MATRIX
                           [-n OUTPUT_PREFIX] [-i INCLUDE_SAMPLES]
                           [-b BETA_LABELS]
                           [--norm-method {none,median,total,control}]
                           [-e NEGATIVE_CONTROL]
                           [--genes-varmodeling GENES_VARMODELING]
                           [--adjust-method {fdr,holm,pounds}] [-o] [-p] [-q]
```

# Arguments #
```
optional arguments:
  -h, --help            show this help message and exit

Required arguments:

  -k COUNT_TABLE, --count_table COUNT_TABLE
                        Provide a tab-separated count table. Each line in the
                        table should include sgRNA name (1st column), target
                        gene (2nd column) and read counts in each sample.
  -d DESIGN_MATRIX, --design_matrix DESIGN_MATRIX
                        Provide a design matrix, either a quoted string of the
                        design matrix or a file name. If using quoted string,
                        for instance, "1,0;1,1" can be used for 2-sample
                        conditions like 0_day and 4_weeks, and --include-
                        samples should be specified as "0_day,4_weeks". If a
                        file is given, the row of the design matrix must match
                        the order of the samples in the count table, or the
                        order of the samples by the --include-samples option.

Optional arguments for input and output:

  -n OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        The prefix of the output file(s).
  -i INCLUDE_SAMPLES, --include-samples INCLUDE_SAMPLES
                        Specify the sample labels if the design matrix is not
                        given by file in the --design-matrix option. Sample
                        labels are separated by ",", and must match the labels
                        in the count table.
  -b BETA_LABELS, --beta-labels BETA_LABELS
                        Specify the labels of the variables (i.e., beta), if
                        the design matrix is not given by file in the
                        --design-matrix option. Should be separated by ",",
                        and the number of labels must equal to (# columns of
                        design matrix), including baseline labels. Default
                        value: "bata_0,beta_1,beta_2,...".

Optional arguments for normalization.:
  "--norm-method control -e non_essential_genes_control" is recommended,
  which you have to specify negative control names, such as AAVS1

  --norm-method {none,median,total,control}
                        Method for normalization, including "none"
                        (nonormalization), "median" (median normalization,
                        default), "total" (normalization by total read
                        counts), "control" (normalization by control sgRNAs
                        specified by the --control-sgrna option). Defautl is
                        median.
  -e NEGATIVE_CONTROL, --negative_control NEGATIVE_CONTROL
                        The name of negative controls, such as AAVS1.
  --genes-varmodeling GENES_VARMODELING
                        The number of genes for mean-variance modeling.
                        Default is 2000.
  --adjust-method {fdr,holm,pounds}
                        Method for sgrna-level p-value adjustment, including
                        false discovery rate (fdr), holm's method (holm), or
                        pounds's method (pounds). Defautl is FDR.

Optional arguments for PPI incorporation and outliers removal:

  -o, --outliers_removal
                        Speicify whehter you want to remove outliers and
                        recalculate..
  -p, --PPI_prior       Specify whether you want to incorporate PPI as prior
  -q, --QC_metric       Specify whether you want to derive quality control
                        metrics

```
# Outputs #
* Log file
* gene_summary.txt
```
 Gene: gene name
 sgRNA: number of targeting sgRNA 
 beta_1|beta: beta score 
 beta_1|z: z-statistic of beta score
 beta_1|wald-p-value: p-value
 beta_1|wald-fdr: false discovery rate
```
* sgrna_summary.txt
```
 Gene: gene name	
 sgRNA: sgrna name	
 eff: whether or not the sgRNA is included in calculation. 1: included; 0: not included	
 mu: mean	
 k: read count  
 loglikelihood  
```
* Histogram of beta score distribution
![Histogram of beta score](https://github.com/hyalin1127/mageck_nest/blob/master/Histogram_of_beta_scores_demo.png)

* The distribution of p value against uniform distribution
![QQplot of wald_p value](https://github.com/hyalin1127/mageck_nest/blob/master/QQplot_of_pvalues.png)

# Demonstration #
* QC only
* Use "AAVS1" as negative controls
```
mageck_nest nest -n mageck_nest_cell_line_A -i day_0,week_4 -k readcount_table.txt --norm-method control -e AAVS1 -d "1,0;1,1" -q
```
* QC
* PPI-integration
* Use "AAVS1" as negative controls
```
mageck_nest nest -n mageck_nest_cell_line_A -i day_0,week_4 -k readcount_table.txt --norm-method control -e AAVS1 -d "1,0;1,1" -q -p
```

# Contact #
Chen-Hao Chen (chen-hao_chen@dfci.harvard.edu)

# Citation #
Chen, C. H. et al. Improved design and analysis of CRISPR knockout screens. Bioinformatics 34, 4095–4101 (2018).
