\name{CoRe.ADaM}
\alias{CoRe.ADaM}
\title{Adaptive Daisy Model to compute core fitness genes}
\usage{
CoRe.ADaM(depMat,
          display=TRUE,
          main_suffix='fitness genes in at least 1 cell line',
          xlab='n. dependent cell lines',
          ntrials=1000,
          verbose=TRUE,
          TruePositives)
}
\arguments{
\item{depMat}{Binary dependency matrix, rows are genes and columns are samples (screens, cell-cell lines). A 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.}

\item{display}{Boolean, default is TRUE. Should bars indicating dependency profiles and boxes for estimated null models be plotte.}

\item{main_suffix}{If display=TRUE, title suffix to be give to the plots.}

\item{xlab}{label to be used in the x-axis of the plots, default is 'n. cell lines'.}

\item{ntrials}{Integer, default =1000. How many times to randomly perturb dependency matrix to generate null distributions of number of genes called essentials in fixed number of cell lines.}
\item{verbose}{Boolean, default is TRUE. Should the computation progress be monitored.}
\item{TruePositives}{Vector of gene symbols to be used as reference prior known essential genes.}
}
\value{
\item{coreFitnessGenes}{A vector of strings with estimated Core Fitness Genes' symbols.}
}
\description{
This function estimates Core Fitness essential genes using the Adaptive Daisy Model [1] starting from a binary gene dependency matrix.
}
\details{
This function identifies Core Fitness essential genes from the joint analysis of multiple CRISPR-Cas9 viability screens performed on different cell-lines / models. It works with binary gene x cell-line essantial/non-essential matrices and it estimates the minimal number \eqn{n} of cell-lines in which a gene should be called as essential in order to be considered as a core-fitness essential gene for the tissue of origin of the screened cell-lines. This threshold is computed in a semi-supervised way and it is defined as that maximising the deviance from expectation of the number of genes that are essential in \eqn{n} cell-lines and their
true positive rates computed with respect to a set of prior known core-fitness essential genes (to be provided in input).
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMat<-CoRe.extract_tissueType_SubMatrix(BinDepMat)

## Loading a reference set of essential genes from
## from the CRISPRcleanR package, derived from [1] and [2]
data(curated_BAGEL_essential)

## Computing lung cancer core-fitness genes with ADaM
cfgenes <- CoRe.ADaM(LungDepMat, TruePositives = curated_BAGEL_essential)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, Mis M, Zimmermann M, Fradet-Turcotte A, Sun S, Mero P, Dirks P, Sidhu S, Roth FP, Rissland OS, Durocher D, Angers S, Moffat J. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26. doi: 10.1016/j.cell.2015.11.015. Epub 2015 Nov 25. PMID: 26627737.
}
\seealso{
\code{\link{CoRe.panessprofile}}
\code{\link{CoRe.generateNullModel}}
\code{\link{CoRe.empiricalOdds}}
\code{\link{CoRe.truePositiveRate}}
\code{\link{CoRe.tradeoffEO_TPR}}
\code{\link{CoRe.coreFitnessGenes}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
