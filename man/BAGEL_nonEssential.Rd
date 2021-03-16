\name{BAGEL_nonEssential}
\alias{BAGEL_nonEssential}
\docType{data}
\title{
Reference set of non essential genes
}
\description{
A list of reference non essential genes assembled from multiple RNAi studies used as classification template
by the BAGEL algorithm to call gene depletion significance [1].
}
\usage{data(BAGEL_nonEssential)}
\format{
 A vector of strings containing HGNC symbols of reference non essential genes.
}
\references{
[1] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, Mis M, Zimmermann M, Fradet-Turcotte A, Sun S, Mero P, Dirks P, Sidhu S, Roth FP, Rissland OS, Durocher D, Angers S, Moffat J. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26
}
\examples{
data(BAGEL_nonEssential)
head(BAGEL_nonEssential)
}
\seealso{
    \code{\link{BAGEL_essential}}
}
\keyword{datasets}
