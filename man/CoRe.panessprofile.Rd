\name{CoRe.panessprofile}
\alias{CoRe.panessprofile}
\title{Profile of number of fitness genes across fixed numbers of cell lines and its cumulative sums}
\usage{
CoRe.panessprofile(depMat,
                   display=TRUE,
                   main_suffix='fitness genes in at least 1 cell line',
                   xlab='n. dependent cell lines')
}
\arguments{
\item{depMat}{A binary dependency matrix, i.e. a binary matrix with genes on rows and samples on columns. A 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.}

\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted}

\item{main_suffix}{If display=TRUE, title suffix to give to plot showing number of genes depleted across fixed number of cell lines, default is 'genes depleted in at least 1 cell line'}

\item{xlab}{If display=TRUE, label to be given to the x-axis of the plots, default is 'n. cell lines'}
}
\value{
A list with the following two named vectors:
\item{panessprof}{Number of genes that are depleted for a number of cell lines}
\item{CUMsums}{Cumulative number of genes depleted in at least x cell lines}
}
\description{
This function calculates the numbers (and cumulative numbers) of genes whose inactivation exerts a fitness effect in \emph{n} cell lines, varying \emph{n} from 1 to the number of cell lines in the dataset in input.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat)

## Computing number of fitness genes across fixed numbers
## of cell lines and its cumulative sums
pprofile <- CoRe.panessprofile(depMat = LungDepMap)
head(pprofile)
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\keyword{functions}
