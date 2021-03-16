\name{BAGEL_essential}
\alias{BAGEL_essential}
\docType{data}
\title{
Reference Core fitness essential genes
}
\description{
A list of reference core fitness essential genes assembled from multiple RNAi studies used as classification template
by the BAGEL algorithm to call gene depletion significance [1].
}
\usage{data(BAGEL_essential)}
\format{
 A vector of strings containing HGNC symbols of reference core fitness essential genes.
}
\references{
[1] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, Mis M, Zimmermann M, Fradet-Turcotte A, Sun S, Mero P, Dirks P, Sidhu S, Roth FP, Rissland OS, Durocher D, Angers S, Moffat J. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26.
}
\examples{
data(BAGEL_essential)
head(BAGEL_essential)
}
\seealso{
    \code{\link{BAGEL_nonEssential}}
}
\keyword{datasets}
