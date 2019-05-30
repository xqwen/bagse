# BAGSE: Bayesian Analysis of Gene Set Enrichment

This repository contains the software package BAGSE designed for gene set enrichment analysis. BAGSE performs both enrichment (hypothesis) testing and quantification. It requires gene-level association evidence (in forms of either z-scores or estimated effect sizes with corresponding standard errors) and pre-defined gene set annotations as input.

The current release is version 1.1.

## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.

## Type of gene set annotation

### 1. Single mutually exclusive gene set annotation

This specific input format should be used if every gene is annotated by one of K mutually exclusive categories. The summary statistics and annotation of the genes should be contained in a single text file. 
Importantly, we require that **the annotation for the baseline category is coded by 0**, other categories can be coded by arbitrary strings or integers.


BAGSE prefers estimated association effect size (b-hat) and its corresponding standard error (sde) as the summary statistics for each gene. The information should be organized in a single text file with the following format

``` 
gene-name  b-hat sde  annotation
```

BAGSE also accepts gene-level z-scores as input. In such case, the expected format for the input text file is

```
gene-name z-score annotation
```
Use command option ``--load_zval`` to inform BAGSE that the z-score input is used. 


#### 1.1  Usage 

```
 bagse  -d input_data [--load_zval] [-fdr_level alpha]  [-fdr_out fdr_output_file]
```


#### 1.2 Sample data

A set of sample data (``sample.bagse.dat``) can be found in the ``src`` folder. To estimate the enrichment parameter run

```
bagse -d sample.bagse.dat --load_zval 
```


### 2.Multiple gene set annotations assuming additivity

#### 2.1  Usage

```
 bagse  -d summary_data -a annot_dat [--load_zval] [-fdr_level alpha]  [-fdr_out fdr_output_file]
```

The presence of the annotation file ``annot_dat`` notifies BAGSE to switch to the algorithm using the additive prior. 


#### 2.2 Sample data

A set of sample data (``sample.dat`` and ``sample.annot.dat``) can be found in the ``src`` folder. To estimate the enrichment parameter run

```
bagse -d sample.dat -a  sample.annot.dat  --load_zval
```




## Output 

### Enrichment estimates

### FDR control results




