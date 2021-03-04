\name{CoRe.download_BinaryDepMatrix}
\alias{CoRe.download_BinaryDepMatrix}
\title{Download Binary Dependency Matrix}
\usage{CoRe.download_BinaryDepMatrix(
URL='https://cog.sanger.ac.uk/cmp/download/binaryDepScores.tsv.zip')}
\arguments{
\item{URL}{url link to download binary matrix.}
}
\value{
Binary matrix
\item{X}{Binary Dependency Matrix introduced in Behan 2019 from Project Score.}
}
\description{Downloading Binary Dependency Matrix introduced in Behan 2019 from Project Score [1].}
\details{
If URL link exists download Binary Dependency Matrix.
}
\examples{
## Downloading Binary Dependency Matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

head(BinDepMat)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
