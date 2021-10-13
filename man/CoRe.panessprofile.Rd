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
\item{depMat}{A binary dependency matrix, i.e. a binary matrix with genes on rows and cell-lines/samples on columns. A 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.}

\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted}

\item{main_suffix}{If display=TRUE, title suffix to give to plot showing number of genes depleted across fixed number of cell lines, default is 'genes depleted in at least 1 cell line'}

\item{xlab}{If display=TRUE, label to be given to the x-axis of the plots, default is 'n. cell lines'}
}
\value{
A list with the following two named vectors:
\item{panessprof}{Number of genes that are depleted across fixed number of cell lines in the inputted dependency matrix}
\item{CUMsums}{Cumulative number of genes depleted in at least \eqn{n} cell lines in the inputted dependency matrix}
}
\description{
This function calculates the numbers (and cumulative numbers) of genes whose inactivation exerts a fitness effect in \emph{n} cell lines, varying \emph{n} from 1 to the number of cell lines in the dependency map given in input.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1,2]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_SubMatrix(BinDepMat)

## Computing number of fitness genes across fixed numbers
## of cell lines and its cumulative sums
pprofile <- CoRe.panessprofile(depMat = LungDepMat)
head(pprofile)
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Dwane L, Behan FM, Gonçalves E, Lightfoot H, Yang W, van der Meer D, Shepherd R, Pignatelli M, Iorio F, Garnett MJ. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.
}
\keyword{functions}
