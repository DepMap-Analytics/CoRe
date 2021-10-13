\name{CoRe.generateNullModel}
\alias{CoRe.generateNullModel}
\title{Null model of number of fitness genes across numbers of cell lines and their cumulative sums}
\usage{
CoRe.generateNullModel(depMat,
                       ntrials=1000,
                       display=TRUE,
                       verbose=TRUE)
}
\arguments{
\item{depMat}{A binary dependency matrix where rows are genes and columns are cell-lines/samples with a 1 in position \emph{[i,j]} indicating that the inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th cell-line/sample.}

\item{ntrials}{Integer, default = 1000. How many randomly permuted versions of the inputted dependency matrix should be generated
to simulate null distributions of number of fitenss genes across number of cell lines and their cumulative sums.}

\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted.}
\item{verbose}{Boolean, default is TRUE. Should a progress bar be displayed during the execution of the function.}
}
\value{
A list with the following two named vectors:
\item{nullProf}{Matrix of number of fitness genes across fixed number of cell lines from. Each rows of matrix corresponds to a randomised version of the inputted matrix.}
\item{nullCumSum}{Matrix of profile of cumulative number of fitness genes in fixed number of cell lines. Each rows of matrix corresponds to a randomised version of the inputted matrix.}
}
\description{
This function randomly perturbs the binary dependency matrix to generate a null distribution of profiles of fitness genes across fixed number of cell lines, and corresponding null distribution of cumulative sums.
}
\details{
For a number of trials specified in (\code{ntrials}) the inputted binary dependency matrix is randomised, keeping its column marginal sums. The profiles of fitness genes across fixed number of cell lines, and corresponding cumulative sums, are returned for each executed random perturbation.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1,2]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_SubMatrix(BinDepMat)

## Simulating Null model for the number of fitness genes
## across numbers of cell lines and their cumulative sums
pprofile <- CoRe.generateNullModel(depMat = LungDepMap)

## Inspecting output
head(pprofile$nullProf)

## Inspecting output
head(pprofile$nullCumSUM)

}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Dwane L, Behan FM, Gonçalves E, Lightfoot H, Yang W, van der Meer D, Shepherd R, Pignatelli M, Iorio F, Garnett MJ. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.
}
\seealso{
\code{\link{CoRe.randomisedepMat}}
\code{\link{CoRe.panessprofile}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
