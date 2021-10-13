\name{CoRe.coreFitnessGenes}
\alias{CoRe.coreFitnessGenes}
\title{Determining Core Fitness from a binary dependency matrix and required minimal number of dependent cell lines.}
\usage{
CoRe.coreFitnessGenes(depMat,
                      crossoverpoint)
}
\arguments{
\item{depMat}{Binary dependency matrix where rows are genes and columns are cell-lines/samples. A 1 in position \emph{[i,j]} indicates that the inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, i.e. that gene is a fitness gene for that cell line}

\item{crossoverpoint}{The estimated minimum number of cell lines in which a gene should be a significant fitness gene in order to be called a core-fitness gene.}
}
\value{
A vector of string containing the predicted core fitness genes.
}
\description{
This function identifies as Core Fitness those genes that are fitness genes in at least \emph{n} cell lines (defined in input), according to the binary dependency matrix defined in input. This minimal \emph{n} is estimated through the ADaM method [1] by the \code{CoRe.tradeoffEO_TPR} function.
}
\examples{
## Downloading a binary dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungdepMat<-CoRe.extract_tissueType_SubMatrix(BinDepMat)

## Compute as core-fitness genes those that are fitness
## in at least 20 Non-Small Cell Lung Carcinoma cell lines
cfgenes <- CoRe.coreFitnessGenes(depMat=LungdepMat,crossoverpoint=20)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\seealso{
\code{\link{CoRe.tradeoffEO_TPR}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
