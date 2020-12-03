\name{CoRe.truePositiveRate}
\alias{CoRe.truePositiveRate}
\title{Profile of True Positive Rates}
\usage{
CoRe.truePositiveRate(depMat,
                      essentialGeneSet)
}
\arguments{
\item{depMat}{Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, i.e. the \emph{i}-th gene is a fitness gene for the \emph{j}-th cell line, 0 otherwise.}

\item{essentialGeneSet}{Reference set of predefined essential genes. This is used to define positive cases.}
}
\value{
A list of the following vectors:
\item{P}{Vector of number of genes that are fitness genes in a given number of cell lines.}
\item{TP}{Vector of number of genes that are fitness genes in a given number of cell lines and are true positives, i.e. in the reference set of essential genes provided in input.}
\item{TPR}{\code{TP} divided by number of genes in set reference set of essential genes provided in input.}
}
\description{
This function calculates a profile of True Positive Rates for fitness genes in at least \emph{n} cell lines, with positive cases from a reference set of essential genes.
}
\details{
This function calculates true positive rates for fitness genes in at least {n} cell lines (for each {n}). First, this function calculates the number of cell lines for which each gene is a fitness gene. Second, for a given number of cell lines, the set of genes that are fitness genes in at least that number of cell lines is determined. Finally, this set of genes is then compared to the reference set of essential genes to calculate a true positive rate.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat)

## Loading a reference set of essential genes from
## from the CRISPRcleanR package, derived from [2]
data(BAGEL_essential)

TPR<-CoRe.truePositiveRate(LungDepMap,BAGEL_essential)
}
\author{
C. Pacini, E. Karakoc & F. Iorio
}
\keyword{functions}
