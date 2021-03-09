\name{CoRe.AssembleFPs}
\alias{CoRe.AssembleFPs}
\title{Assembling expression-based false positives}
\usage{CoRe.AssembleFPs(URL='https://ndownloader.figshare.com/files/26261476')}
\arguments{
\item{URL}{URL of the CCLE gene expression dataset.}
}
\value{
\item{LowlyExpr}{A vector of strings with symbols of genes that are lowly expressed across ~1,300 cell lines to be used as false positive by the benchmarking function.}
}
\description{Assembling a set of lowly expressed genes cancer cell lines from the CCLE [1].}
\details{
This function download CCLE gene expression data from DepMap portal [1]  then it estimates a set of overall lowly expressed genes as detailed in [2].
}
\examples{
FPs<-CoRe.AssembleFPs()
head(FPs)
}
\references{
[1] Barretina, J., Caponigro, G., Stransky, N. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603–607 (2012).

[2] Pacini, Dempster et al, Integrated cross-study datasets of genetic dependencies in cancer. https://doi.org/10.1101/2020.05.22.110247
.
}
\seealso{
\code{\link{CoRe.CF_Benchmark}}
}
\author{
A. Vinceti & F. Iorio
}
\keyword{functions}
