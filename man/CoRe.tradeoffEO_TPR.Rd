\name{CoRe.tradeoffEO_TPR}
\alias{CoRe.tradeoffEO_TPR}
\title{Calculate AdAM threshold}
\usage{
CoRe.tradeoffEO_TPR(EO,
                    TPR,
                    test_set_name,
                    display = TRUE)
}
\arguments{
\item{EO}{Profile of empirical odds values. Computed with the \code{CoRe.empiricalOdds} function.}

\item{TPR}{Profile of True positive rates for across number of cell line. Computed with the \code{CoRe.truePositiveRate} function.}

\item{test_set_name}{Name to give to the analysis, used for plotting titles.}
\item{display}{Boolean, default is TRUE. Should AdAM tradeoff strategy be plotted.}
}
\value{
ADAM model threshold:
\item{point}{Number of cell lines for which a gene needs to be a fitness gene in order to be predicted as core-fitness gene.}
}
\description{
This function finds the minimum number of cell lines in which a gene needs to be fitness in order to be called core-fitness
for all the considered cell lines. This is defined as the \emph{n} providing the best trade-off between i) coverage of priori-known essential genes in the resulting set of predicted core-fitness genes, i.e. fitness in at least \emph{n} cell lines, and ii) deviance from expectation of the number of fitness genes in \emph{n} cell lines.
}
\details{
Compare and plot the log10 odds ratios with the true positive rates to find the cross over point where the true positive rate falls below the odds ratio.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMat<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat,tissue_type="Non-Small Cell Lung Carcinoma")

## Loading a reference set of essential genes from
## from the CRISPRcleanR package, derived from [2]
data(BAGEL_essential)

# Generate the profiles of number of fitness genes across number of cell lines from
# observed data and corresponding comulative sums.
pprofile<-CoRe.panessprofile(depMat=LungDepMat)

# Generate a set of random profiles of number of genes depleted for a number of cell lines
# and corresponding cumulative sums by perturbing observed data.
nullmodel<-CoRe.generateNullModel(depMat=LungDepMat,ntrials = 1000)

# Calculate log10 odd ratios of observed/expected profiles of cumulative number of fitness
# genes in fixed number of cell lines.
# Observed values are from the ADAM.panessprofile function and expected are the average of
# random set from CoRe.generateNullModel
EO<-CoRe.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )

# Calculate True positive rates for fitness genes in at least n cell lines in the observed
# dependency matrix, with positive cases from a reference set of essential genes
TPR<-CoRe.truePositiveRate(LungDepMat,BAGEL_essential)

# Calculate minimum number of cell lines a gene needs to be a fitness gene in order to
# be considered as a core-fitness gene
crossoverpoint<-CoRe.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'BAGEL essential')
crossoverpoint
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, Mis M, Zimmermann M, Fradet-Turcotte A, Sun S, Mero P, Dirks P, Sidhu S, Roth FP, Rissland OS, Durocher D, Angers S, Moffat J. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26. doi: 10.1016/j.cell.2015.11.015. Epub 2015 Nov 25. PMID: 26627737.
}
\seealso{
\code{\link{CoRe.empiricalOdds}},
\code{\link{CoRe.truePositiveRate}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
