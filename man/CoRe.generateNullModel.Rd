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
\item{depMat}{Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.}

\item{ntrials}{Integer, default = 1000. How many times to randomly perturb dependency matrix to generate the null distributions.}

\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted}
\item{verbose}{Boolean, default is TRUE. Should a progress bar be displayed}
}
\value{
A list with the following two named vectors:
\item{nullProf}{Matrix of number of fitness genes for fixed number of cell lines from. Each rows of matrix corresponds to a randomisation trial.}
\item{nullCumSum}{Matrix of profile of cumulative number of fitness genes in fixed number of cell lines. Each row of matrix is one randomisation trial.}
}
\description{
This function randomly perturbs the binary dependency matrix to generate a null distribution of profiles of fitness genes across fixed number of cell lines, and corresponding null distribution of cumulative sums.
}
\details{
For a number of trials specified in (\code{ntrials}) the inputted binary dependency matrix is randomised, keeping its column marginal sums. The profiles of fitness genes across fixed number of cell lines, and corresponding cumulative sums, are returned for each random perturbation.
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat)

## Simulating Null model for the number of fitness genes
## across numbers of cell lines and their cumulative sums
pprofile <- CoRe.generateNullModel(depMat = LungDepMap)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
C. Pacini, E. Karakoc & F. Iorio
}
\keyword{functions}
