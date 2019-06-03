# BAGSE: Bayesian Analysis of Gene Set Enrichment

This repository contains the software package BAGSE designed for gene set enrichment analysis. BAGSE performs both enrichment (hypothesis) testing and quantification. It requires gene-level association evidence (in forms of either z-scores or estimated effect sizes with corresponding standard errors) and pre-defined gene set annotations as input.

The current release is version 1.1.

## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.

## Types of gene set annotation

When multiple gene set annotations are available, there are two possible ways to formulate gene set enrichment analysis involving multiple gene sets. 

Note that when only a single gene set is used for analysis, both approaches are applicable  but the second approach makes stronger assumption (in requiring the *same* effect size distributions of the associated genes for annotated and un-annotated genes). Thus, we recommend the first approach for analyzing a single gene set annotation.  



### 1. Single mutually exclusive gene set annotation

In this approach, we define the combination of multiple gene set annotation as a new annotation. Consider two potentially overlapping gene set annotations, there are 4 possible combination of annotations depending on the presence and absence of a gene in each gene set, i.e.,

```
0 0  ---> 0
0 1  ---> 1
1 0  ---> 2
1 1  ---> 3
```
For example ``0 1`` denotes a gene annotated in the second gene set but not in the first set, and we denote this combination as category ``1`` in the combined annotation. 
More generally, this specific input format should be used if every gene is annotated by one of K mutually exclusive categories. The summary statistics and annotation of the genes should be contained in a single text file. 
Importantly, we require that **the annotation for the baseline category is coded by 0**, other categories can be coded by arbitrary strings or integers.

Under this formulation, BAGSE allows that the effect size distribution under the alternative model (i.e., when gene-level association is genuine) is category-specific. 

Note, although this approach is the most general, it does not scale computationally for a large number of gene sets.  


BAGSE prefers estimated association effect size (b-hat) and its corresponding standard error, se(b-hat), as the summary statistics for each gene. The information should be organized in a single text file with the following format

``` 
gene-name  b-hat se(b-hat)  annotation
```

BAGSE also accepts gene-level z-scores as input. In such case, the expected format for the input text file is

```
gene-name  z-score  annotation
```
Use command option ``--load_zval`` to inform BAGSE that the z-score input is used. 

Additionally, p-values can be used in place of z-scores when the command line option ``--load_pval`` is specified.

#### 1.1  Usage 

```
 bagse  -d input_data [--load_zval | --load_pval ] [-fdr_level alpha]  [-fdr_out fdr_output_file]
```


#### 1.2 Sample data

A sample data set (``sample.combination_annotation.dat``) can be found in the ``sample_data`` folder. To estimate the enrichment parameter run

```
bagse -d sample.combination_annotation.dat --load_zval 
```


### 2. Multiple gene set annotations assuming additivity

An alternative approach is to assume that multiple gene set annotation is additively affecting the odds of a given gene being associated. This is a simplifying assumption in comparison to the previous approach, but has the advantage in computational efficiency, which allows to consider many gene sets simultaneously. 
Additionally, this approach assumes a single distribution of effects under the alternative model regardless of types of annotations, which is more restrictive. 
 
The input data for this approach are separated into two text files: one contains gene-level association statistics and the other contains gene set annotation information.  
The summary statistics file has the following format:

```
gene-name  b-hat se(b-hat)
```

Z-scores are also accepted as input by using the command line option ``--load_zval`` and the following format:

```
gene-name z-score
```
Additionally, p-values can be used in place of z-scores when the command line option ``--load_pval`` is specified.



The annotation file has the following format:

```
gene-name annotation1 annotation2 annotation3 ....
```
Importantly, a header starting with the keyword "Gene" is required:

```
Gene    set1-name set2-name set3-name ...
```



#### 2.1  Usage

```
 bagse  -d summary_data -a annot_dat [--load_zval|--load_pval] [-fdr_level alpha]  [-fdr_out fdr_output_file]
```

**The presence of the annotation file** ``annot_dat`` **and the**``-a`` **flag notify  BAGSE to switch to the algorithm using the additive prior.** 


#### 2.2 Sample data

A set of sample data in this format (``sample.additive_summary.dat`` and ``sample.additive_annot.dat``) can be found in the ``sample_data`` folder. To estimate the enrichment parameter run

```
bagse -d sample.additive_summary.dat -a  sample.additive_annot.dat  --load_zval
```




## Output 

### Enrichment estimates

Under the formulation of single mutually exclusive gene set annotation, the output of enrichment estimates from analyzing the sample data is given by 

```
       Baseline            0    -0.973     -1.075  -0.871
        annot.1   category_1     0.841      0.724   0.958
        annot.2   category_2     1.019      0.904   1.134
        annot.3   category_3     2.261      2.106   2.416
```

+ The first and second columns denote the names of the parameter and their corresponding annotations provided by the user. 
+ The third column shows the maximum likelihood estimates (MLEs) of the enrichment parameters. For baseline level, the estimate provides the log-odds ratio for category "0", which is served as the contrast for all other mutually exclusive categories. The estimates of the other categories are the estimated log-odds ratios *relative to* the baseline level, i.e., representing the relative enrichment in comparison to the baseline level. 
For example, in the above example, the absolute log-odds ratio for ``category_1`` can be computed by -0.973 + 0.841 = -0.132.
+ The last two columns jointly represents the 95% confidence intervals for the estimated parameters.



Under the formulation of multiple gene set annotations assuming additivity, the output is similar. The only difference is that the first two columns from the previous output are consolidated into one parameter column. The output from running the sample data is given by 

```
      Intercept   -1.065      -1.132  -0.997
         set1.1    1.001       0.905   1.097
         set2.1    1.207       1.112   1.301
```



### FDR control results


When the command line option ``-fdr_level alpha`` is specified, BAGSE perform FDR control to identify associated genes incorporating estimated enrichment information at ``alpha`` (e.g., 0.05) level.
If an output file is *not* specified through ``-fdr_out`` option, a default ``fdr.out`` file will be created and used to record the FDR control results.

The output of the FDR control has the following format:

```
          Gene1    4.787e-05   1
          Gene2    8.006e-01   0
          Gene3    7.301e-01   0
          Gene4    7.514e-01   0
          Gene5    3.472e-07   1

```

+ Column 1: gene name
+ Column 2: local fdr or posterior false discovery probability
+ Column 3: rejection at the defined alpha level (1 means reject)




