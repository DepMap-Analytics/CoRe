# CoRe
![alt text](https://github.com/DepMap-Analytics/CoRe/blob/master/web/coRe_logo.jpg)

The CoRe package implements two methods for identification of core fitness genes (at two level of stringency) from joint analyses of multiple genome-wide CRISPR-Cas9 screens: 

1) The percentile ranking method builds on the starting basic intuition that if a gene is involved in essential processes in every cell, it should be among the top essential genes in all the analysed screens, even those performed on the cell lines whose viability is lowly impacted by the knock-out of that gene, i.e. the least dependent cell lines. The CoRe package implements four versions of this method, which verify if the starting assumption is satisfied for (considering each gene in turn) in different ways. The first variant, named 'fixed', accounts for the distribution of the rank position of the genes in their least 90-th percentile dependent cell lines, obtained by sorting all genes based on their essentiality in that cell line. The second, named 'average', considers the distribution of the average gene rank positions in all the cell lines falling at least in 90-th percentile cell lines of least dependent ones.
The third version, named 'slope', fits a linear model on the curve obtained when sorting all cell lines based on how much they are dependent on a given gene then considering the essentiality rank position of that gene across all cell lines (in the obtained order). The fourth version, named 'AUC', the Area Under the Curve resulting from considering the sequence of gene fitness-rank-positions across all cell lines, sorted according to their dependency on that gene, is used in the subsequent step Percentile method.
The density of the considered distributions is then estimated using a Gaussian kernel with width 0.1 and
the point of minimum density is then identified. Genes whose rank position, average rank position, slope or AUC (according to the chosen variant) is predicted to be generated by the first distribution are called common essential genes, i.e. core fitness genes at a lower stringency [1]. 

2) The Adaptive Daisy Model (ADaM) is a semi-supervised algorithm for computing a fuzzy-intersection of non-fuzzy sets by adaptively determining the minimal number of sets to which an element should belong in order to be a member of the fuzzy-intersection (the membership threshold). This threshold maximises the deviance from expectation of the cardinality of the resulting fuzzy-intersection, as well as the coverage of predefined elements.
This method can be used to identify the minimal number of cell lines from a given tissue in which the inactivation of a gene (via CRISPR-Cas9 targeting) should exert a reduction of viabilty (or fitness effect) in order for that gene to be considered a core-fitness essential gene for the tissue under consideration. Iterating this method and considering cancer-type specific core fitness genes as sets to be fuzzy-intersected, then this method can be used to estimate a set of pan-cancer core fitness genes.
This method is used to discriminate between core-fitness and context-specific essential genes in a study describing a large scale genome-wide CRISPR-Cas9 pooled drop-out screening [1] (a detailed description of the algorithm is included in the Supplemental Information of [2]). ADaM was inspired by the Daisy Model method introduced in [3]

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

[1] J. M. Dempster et al., Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets., Nat. Commun., vol. 10, no. 1, p. 5817, 2019, doi: 10.1038/s41467-019-13805-y.

[2] Behan FM & Iorio F & Picco G et al., Prioritisation of cancer therapeutic targets using CRISPR-Cas9 screens. Nature, In press.

[3] Hart T et al., High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015;163:1515–26.
