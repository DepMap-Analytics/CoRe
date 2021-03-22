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
This function takes in input a matrix and shuffles its entries column-wise. Then matrix resulting from this shuffling will have the same column marginal totals of the inputted one.
}
\examples{

## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1,2]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMat<-CoRe.extract_tissueType_SubMatrix(BinDepMat)

## Randomising the colums of the submatrix
rnd_exampleDepMat<-CoRe.randomisedepMat(LungDepMat)

}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Dwane L, Behan FM, Gonçalves E, Lightfoot H, Yang W, van der Meer D, Shepherd R, Pignatelli M, Iorio F, Garnett MJ. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.

}
\keyword{functions}
