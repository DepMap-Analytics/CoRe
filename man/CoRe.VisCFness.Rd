\name{CoRe.VisCFness}
\alias{CoRe.VisCFness}
\title{Common-essentiality tendency visualisation}
\usage{CoRe.VisCFness(depMat, gene,
                      percentile=0.9,
                      posControl='RPL8',
                      negControl='MAP2K1',
                      method='fixed')}
\arguments{
\item{depMat}{Quantitative Dependency Matrix containing Pan-cancer or tissue/cancer-types specific gene fitness/dependency scores across cell-lines/samples. The value in position \eqn{[i,j]} of such matrix quantifies the fitness/dependency score of the \eqn{i}-th gene in the \eqn{j}-th cell line.}
\item{gene}{Character, name of the gene whose common-essentiality tendency should be visualised.}

\item{method}{Character, default is 'fixed'. This parameter specifies which variant of the Fitness Percentile method should be used. Admissimble values are:\cr
  - fixed: a distribution of gene fitness-rank-positions in their least dependent \code{n}-th (determined by the \code{percentile} parameter) percentile cell line is used in the subsequent step of the Fitness Percentile method.\cr
  - average: a distribution of gene average fitness-rank-position across all cell lines at or over the \code{n}-th percentile of least dependent cell lines (where \eqn{n} is determined by the \code{percentile} parameter) is used in the subsequent step of the Fitness Percentile method.\cr
  - slope: for each gene, a linear model is fit on the sequence of gene fitness-rank-positions across all cell lines sorted according to their dependency on that gene, then a distrubution of models' slopes is used in the subsequent step of the  Fitness Percentile method.\cr
  - AUC: for each gene, the area under the curve resulting from considering the sequence of gene fitness-rank-positions across all cell lines sorted according to their dependency on that gene is used in the subsequent step of the Fitness Percentile method.}

\item{percentile}{Numerical value in range [0,1], default is 0.9. Percentile to be used as a comparative threshold.}
\item{posControl}{Character, name of a gene used as positive control for the visualization.}
\item{negControl}{Character, name of gene used as negative control for the visualization.}
}
\description{This function visualises the tendency of a given gene to be common essential (in the dependency matrix provided in input according the the Fitness Percentile method described in [1]) and compares this tendency to that of a positive and a negative control. It also visualises the distribution of scores computed by the Fitness percentile methods (according to the selected variant) with the discriminative threshold and rank position of the input gene as well as the control genes.}
\examples{
## Downloading and scaling quantitative dependency matrix from project score [2,3]
DepMat<-CoRe.download_DepMatrix()

## CFness visualisation of BRAF (a context specific essential gene)
CoRe.VisCFness(DepMat,'BRAF')

## CFness visualisation of RPL22 (a common essential gene)
CoRe.VisCFness(DepMat,'RPL22')
}
\references{
[1] Dempster, J.M., Pacini, C., Pantel, S. et al. Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets. Nat Commun 10, 5817 (2019).

[2] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[3] Dwane L, Behan FM, Gonçalves E, Lightfoot H, Yang W, van der Meer D, Shepherd R, Pignatelli M, Iorio F, Garnett MJ. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.
}

\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
