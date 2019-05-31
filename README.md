# BAGSE: Bayesian Analysis of Gene Set Enrichment

This repository contains the software package BAGSE designed for gene set enrichment analysis. BAGSE performs both enrichment (hypothesis) testing and quantification. It requires gene-level association evidence (in forms of either z-scores or estimated effect sizes with corresponding standard errors) and pre-defined gene set annotations as input.

The current release is version 1.1.

## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.

## Types of gene set annotation

When multiple gene set annotations are available, there are two possible ways to formulate gene set enrichment analysis involving multiple gene sets. When only a single gene set is used for analysis, both approaches apply but they yield slightly different results.



### 1. Single mutually exclusive gene set annotation

In this approach, we define the combination of multiple gene set annotation as a new annotation. Consider two potentially overlapping gene set annotations, there are 4 possible combination of annotations depending on the presence and absence of a gene in each gene set, i.e.,

```
0 0  ---> 0
1 0  ---> 1
0 1  ---> 2
1 1  ---> 3
```
For example ``0 1`` denotes a gene annotated in the second gene set but not in the first set, and we denote this combination as category ``3`` in the combined annotation. 
More generally, this specific input format should be used if every gene is annotated by one of K mutually exclusive categories. The summary statistics and annotation of the genes should be contained in a single text file. 
Importantly, we require that **the annotation for the baseline category is coded by 0**, other categories can be coded by arbitrary strings or integers.

Under this formulation, BAGSE allows the effect size distribution under the alternative model (i.e., when gene-level association is geniune) is category-specific. 

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


#### 1.1  Usage 

```
 bagse  -d input_data [--load_zval] [-fdr_level alpha]  [-fdr_out fdr_output_file]
```


#### 1.2 Sample data

A sample data set (``sample.combination_annotation.dat``) can be found in the ``sample_data`` folder. To estimate the enrichment parameter run

```
bagse -d sample.combination_annotation.dat --load_zval 
```


### 2.Multiple gene set annotations assuming additivity

#### 2.1  Usage

```
 bagse  -d summary_data -a annot_dat [--load_zval] [-fdr_level alpha]  [-fdr_out fdr_output_file]
```

The presence of the annotation file ``annot_dat`` and the ``-a`` flag notify  BAGSE to switch to the algorithm using the additive prior. 


#### 2.2 Sample data

A set of sample data in this format (``sample.additive_summary.dat`` and ``sample.additive_annot.dat``) can be found in the ``sample_data`` folder. To estimate the enrichment parameter run

```
bagse -d sample.additive_summary.dat -a  sample.additive_annot.dat  --load_zval
```




## Output 

### Enrichment estimates

### FDR control results




