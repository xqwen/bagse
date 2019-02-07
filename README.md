# BAGSE: Bayesian Analysis of Gene Set Enrichment

This repository contains the software package BAGSE designed for gene set enrichment analysis. BAGSE performs both enrichment (hypothesis) testing and quantification. It requires gene-level association evidence (in forms of either z-scores or estimated effect sizes with corresponding standard errors) and pre-defined gene set annotations as input.


## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.

## Usage 

```
 bagse  -d input_data -annot gene_set_annotation [--load_bf | --load_zval] [-lfdr lfdr_output_file]
```

