\name{CoRe.randomisedepMat}
\alias{CoRe.randomisedepMat}
\title{Column-wise matrix randomisation}
\usage{
CoRe.randomisedepMat(depMat)
}
\arguments{
\item{depMat}{A numeric matrix}
}
\value{
The matrix given in input with entries shuffled column wisely.
}
\description{
This function takes in input a matrix and shuffles its entries column wisely. Then matrix resulting from this shuffling will have the same column marginal totals of the inpputted one.
}
\examples{

## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat)

## Randomising the colums of the submatrix
rnd_exampleDepMat<-CoRe.randomisedepMat(LungDepMap)

}
\author{
C. Pacini, E. Karakoc & F. Iorio
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\keyword{functions}
