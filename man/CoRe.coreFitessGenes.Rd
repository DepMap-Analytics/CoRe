\name{CoRe.coreFitnessGenes}
\alias{CoRe.coreFitnessGenes}
\title{Calculate the Core Fitness genes given the binary dependency matrix and the minimal number of cell line threshold.}
\usage{
CoRe.coreFitnessGenes(depMat,
                      crossoverpoint)
}
\arguments{
\item{depMat}{Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.}

\item{crossoverpoint}{minimum number of cell lines in which a gene needs to be fitness in order to be called core-fitness}
}
\value{
A vector that containing the Core Fitness Genes:
}
\description{
This function identifies as Core Fitness the genes that are fitness in a number of cell lines at least equal to the inputted trheshold, this is computed by the \code{CoRe.tradeoffEO_TPR} function.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat)

## Compute as core-fitness genes those that are fitness
## in at least 20 lung cance cell lines
cfgenes <- CoRe.coreFitnessGenes(depMat=exampleDepMat,crossoverpoint=3800)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\seealso{
\code{\link{CoRe.tradeoffEO_TPR}},
}
\author{
C. Pacini, E. Karakoc & F. Iorio
}
\keyword{functions}
