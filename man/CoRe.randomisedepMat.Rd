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
This function takes in input a matrix and shuffles its entries column wisely. If the matrix is binary then then matrix resulting from this shuffling will have the same column marginal totals of the inpputted one.
}
\examples{
data(exampleDepMat)
rnd_exampleDepMat<-CoRe.randomisedepMat(exampleDepMat)
}
\author{
C. Pacini, E. Karakoc & F. Iorio
}
\keyword{functions}
