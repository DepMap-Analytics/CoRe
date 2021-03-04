\name{curated_BAGEL_essential}
\alias{curated_BAGEL_essential}
\docType{data}
\title{
Curated Reference Core fitness essential genes
}
\description{
A list of reference core fitness essential genes assembled from multiple RNAi studies used as classification template
by the BAGEL algorithm to call gene depletion significance [1] and filtered as shown in [2].
}
\usage{data(curated_BAGEL_essential)}
\format{
 A vector of strings containing HGNC symbols of curated reference core fitness essential genes.
}
\references{
[1] BAGEL: a computational framework for identifying essential genes from pooled library screens. Traver Hart and Jason Moffat. BMC Bioinformatics, 2016 vol. 17 p. 164.

[2] Behan, F.M., Iorio, F., Picco, G. et al. Prioritization of cancer therapeutic targets using CRISPR–Cas9 screens. Nature 568, 511–516 (2019).
}
\examples{
data(curated_BAGEL_essential)
head(curated_BAGEL_essential)
}
\seealso{
    \code{\link{curated_BAGEL_nonEssential}}
}
\keyword{datasets}
