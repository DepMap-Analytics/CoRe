\name{CoRe.AssembleFPs}
\alias{CoRe.AssembleFPs}
\title{Assemble expression-based false positives}
\usage{CoRe.AssembleFPs(URL='https://ndownloader.figshare.com/files/26261476')}
\arguments{
\item{URL}{url link to download CCLE gene expression data.}
}
\value{
\item{LowlyExpr}{Vector of lowly expressed genes used as false positive for benchmarcking function.}
}
\description{Download CCLE gene expression data [1].}
\details{
Download CCLE gene expression data from DepMap portal [1].
}
\examples{
FPs<-CoRe.AssembleFPs(URL='https://ndownloader.figshare.com/files/26261476')
head(FPs)
}
\references{
[1] Barretina, J., Caponigro, G., Stransky, N. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603–607 (2012).
}
\seealso{
\code{\link{CoRe.CF_Benchmark}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
