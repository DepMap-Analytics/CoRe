# CoRe
![alt text](https://github.com/DepMap-Analytics/CoRe/blob/master/web/coRe_logo.jpg)

The CoRe package implements two main methods for identification of essential genes using genome-wide CRISPR-Cas9 screens. The main functions of this package can be summarized as follows:

1) Percentile gene score ranking method which aims to identify the basic intuition that if a gene is involved in essential processes in cells, it should fall in the top N depleted genes in at least 90% of cell lines. There are three versions of this method implemented in CoRe. First version calculates the rank distribution of genes at their least dependent 90-th percentile cell line and the second one calculates the average rank distribution of the least dependent 90-th percentile cell lines. The third version, instead, employs all the cell lines to fit a linear model to generate a gene score rank distribution for each gene. The slope distribution of the genes score ranks form a bimodal distribution similar to the percentile gene score rank methods of the first two versions. Using this distribution it is possible to determine the point of minimum density between two peaks and predict the essential genes. For each version, the resulting bimodal distribution is used for the identification of essential genes under the considered cell lines [1]. 

2) Adaptive Daisy Model which is a semi-supervised algorithm for computing a fuzzy-intersection of non-fuzzy sets by adaptively determining the minimal number of sets to which an element should belong in order to be a member of the fuzzy-intersection (the membership threshold). This threshold maximises the deviance from expectation of the cardinality of the resulting fuzzy-intersection, as well as the coverage of predefined elements. This method can be used to identify the minimal number of cell lines from a given tissue in which the inactivation of a gene (for example via CRISPR-Cas9 targeting) should exert a reduction of viability (or fitness effect) in order for that gene to be considered a core-fitness essential gene for the tissue under consideration [2].

Contributors: Emre Karakoc, Clare Pacini, Alessandro Vinceti & Francesco Iorio


Install
--

Install the R package via devtools:

```
install.packages("devtools")
library(devtools)

install_github("DepMap-Analytics/CoRe")
```

**References**

[1]    J. M. Dempster et al., “Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets.,” Nat. Commun., vol. 10, no. 1, p. 5817, 2019, doi: 10.1038/s41467-019-13805-y.

[2]  Behan FM & Iorio F & Picco G et al., Prioritisation of cancer therapeutic targets using CRISPR-Cas9 screens. Nature, In press.
