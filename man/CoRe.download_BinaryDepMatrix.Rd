\name{CoRe.download_BinaryDepMatrix}
\alias{CoRe.download_BinaryDepMatrix}
\title{Download Binary Dependency Matrix}
\usage{CoRe.download_BinaryDepMatrix(
URL='https://cog.sanger.ac.uk/cmp/download/binaryDepScores.tsv.zip')}
\arguments{
\item{URL}{URL pointing to the online location of a zipped folder containing a binary dependency matrix.
By default this will point to an entry on the data download page of the Project Score portal with data from [1].}
}
\value{
A binary Dependency Matrix (from [1]) where rows are genes and columns are cell lines. A 1 in position [i,j] indicates that the inactivation of the i-th gene exerts a significant loss of fitness in the j-th sample, 0 otherwise.
}
\description{Downloading Binary Dependency Matrix introduced in [1] from the Project Score portal [2].}
\details{
If the URL points to a valid online location this function downloads a Binary Dependency Matrix.
}
\examples{
## Downloading Binary Dependency Matrix
## for > 300 cancer cell lines from [1,2]
BinDepMat<-CoRe.download_BinaryDepMatrix()

head(BinDepMat)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Dwane L, Behan FM, Gonçalves E, et al. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
