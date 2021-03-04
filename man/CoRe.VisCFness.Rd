\name{CoRe.VisCFness}
\alias{CoRe.VisCFness}
\title{Visualization of CFness of a gene}
\usage{CoRe.VisCFness(depMat,
                         gene,
                         percentile=0.9,
                         posControl='RPL12',
                         negControl='MAP2K1')}
\arguments{
\item{depMat}{Quantitative Dependency Matrix containing all cell models.}
\item{gene}{Character, name of gene of interest.}
\item{percentile}{Numerical value in range [0,1], default is 0.9. Percentile to be used as threshold.}
\item{posControl}{Name of gene used as positive control for the visualization.}
\item{negControl}{Name of gene used as negative control for the visualization.}
}
\description{Visualization of CFness of a gene and comparison to positive and negative control gene.}
\examples{
## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1]
DepMat<-CoRe.download_DepMatrix()

gene<-"A1BG"

## CFness visualization
CoRe.VisCFness(DepMat,
               gene,
               percentile=0.9,
               posControl='RPL12',
               negControl='MAP2K1')
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
