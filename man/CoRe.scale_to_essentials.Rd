\name{CoRe.scale_to_essentials}
\alias{CoRe.scale_to_essentials}
\title{Scale matrix cell-wise using positive and negative reference gene sets}
\usage{
CoRe.scale_to_essentials(ge_fit,
                         ess_genes,
                         noness_genes)
}
\arguments{
\item{ge_fit}{Dependency matrix, rows are genes and columns are samples (screens, cell-cell lines).}

\item{ess_genes}{Character vector, list of reference essential genes.}

\item{noness_genes}{Character vector, list of reference nonessential genes.}
}
\value{
\item{scaled_ge_fit}{Scaled dependency matrix.}
}
\description{
This function takes as input a continuous dependency matrix and perform a cell-wise (i.e. column-wise) scaling by employing two reference sets of essential and nonessential genes.
}
\details{
This function takes as input a continuous dependency matrix and perform a cell-wise (i.e. column-wise) scaling by setting, in every cell line, the median of reference essential genes to -1 and the median of reference non-essential genes to 0.
}
\examples{
data('curated_BAGEL_essential')
data('curated_BAGEL_nonEssential')

DepMat<-CoRe.download_DepMatrix()
scaled_DepMat<-CoRe.scale_to_essentials(DepMat,
                                        curated_BAGEL_essential,
                                        curated_BAGEL_nonessential)
\references{
[1] J. M. Dempster et al., Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets., Nat. Commun., vol. 10, no. 1, p. 5817, 2019, doi: 10.1038/s41467-019-13805-y.

[2] Clare Pacini et al., Integrated cross-study datasets of genetic dependencies in cancer., Nat. Commun., vol. 12, p. 1661, 2021, doi: 10.1038/s41467-021-21898-7.
}
\seealso{
\code{\link{CoRe.download_DepMatrix}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
