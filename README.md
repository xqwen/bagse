# BAGSE: Bayesian Analysis of Gene Set Enrichment

This repository contains the software package BAGSE designed for gene set enrichment analysis. BAGSE performs both enrichment (hypothesis) testing and quantification. It requires gene-level association evidence (in forms of either z-scores or estimated effect sizes with corresponding standard errors) and pre-defined gene set annotations as input.


## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.

## Usage 

```
 bagse  -d input_data -annot gene_set_annotation [--load_zval] [-lfdr lfdr_output_file]
```

## Input data format

### Summary statistics input

The summary statistics input contains individual gene association evidence. BAGSE prefers using estimated association effect size (b-hat) and its corresponding standard error (sde) for each gene. The information should be organized in a single text file with the following format

``` 
gene-name  b-hat sde 
```

BAGSE also accepts gene-level z-scores as input. In such case, the expected format for the input text file is

```
gene-name z-score
```


### Gene set annotation

The gene set annotation should be organized in a text file, which requires the following header
```
Gene  pathway-name
```
The pathway-name can  be replaced with any actual gene set name. For data input, list each gene in the input data file. Use ``1`` to denote a member gene in the gene set and ``0`` otherwise.


### Sample data

A set of sample data can be found in the ``src`` folder.




