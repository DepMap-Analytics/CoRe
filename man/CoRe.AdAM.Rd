\name{CoRe.AdAM}
\alias{CoRe.AdAM}
\title{Adaptive Daisy Model to compute core fitness genes}
\usage{
CoRe.AdAM(depMat,
          display=TRUE,
          main_suffix='fitness genes in at least 1 cell line',
          xlab='n. dependent cell lines',
          ntrials=1000,
          verbose=TRUE)
}
\arguments{
\item{depMat}{Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.}

\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted}

\item{main_suffix}{If display=TRUE, title suffix to give to plot of number of genes depleted in a give number of cell lines, default is 'genes depleted in at least 1 cell line'}

\item{xlab}{label to give to x-axis of the plots, default is 'n. cell lines'}

\item{ntrials}{Integer, default =1000. How many times to randomly perturb dependency matrix to generate the null distributions.}
\item{verbose}{Boolean, default is TRUE. Should the computation progress be monitored.}
}
\value{
A vector of string with Core Fitness Genes' symbols:
}
\description{
This function identifies the Core Fitness genes using the Adaptive Daisy Model [1] starting from a binary dependency matrix.
}
\details{
This function calculates the Core Fitness essential genes based on the calculated minimum number of cell lines that optimizes the True positive rates with log10 odds ratios. log10 odd ratios are calculated of observed vs. expected profiles of cumulative number of fitness genes in fixed number of cell lines. Expected values are the mean of those observed across randomised version of the observed binary matrix.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat)

## Loading a reference set of essential genes from
## from the CRISPRcleanR package, derived from [2]
data(BAGEL_essential)

## Computing lung cancer core-fitness genes with AdAM
cfgenes <- CoRe.AdAM(LungDepMap)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, Mis M, Zimmermann M, Fradet-Turcotte A, Sun S, Mero P, Dirks P, Sidhu S, Roth FP, Rissland OS, Durocher D, Angers S, Moffat J. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26. doi: 10.1016/j.cell.2015.11.015. Epub 2015 Nov 25. PMID: 26627737.
}
\author{
C. Pacini, E. Karakoc & F. Iorio
}
\keyword{functions}
