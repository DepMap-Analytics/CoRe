\name{CoRe.download_DepMatrix}
\alias{CoRe.download_DepMatrix}
\title{Download Quantitative Dependency Matrix}
\usage{CoRe.download_DepMatrix(
                     URL='https://cog.sanger.ac.uk/cmp/download/essentiality_matrices.zip',
                     scaled=FALSE,
                     ess,
                     noness)}
\arguments{
\item{URL}{url link to download quantitative dependency matrix.}
\item{scaled}{Boolean, default is FALSE. Should the Quantitative Dependency matrix be scaled using CERES strategy.}
\item{ess}{Vector of gene symbols to be used as prior known essential genes.}
\item{noness}{Vector of gene symbols to be used as prior known non-essential genes.}
}
\value{
Quantitative Dependency Matrix
\item{X}{Quantitative Dependency Matrix introduced in Behan 2019 from Project Score.}
}
\description{Downloading Quantitative Dependency Matrix introduced in Behan 2019 from Project Score [1].}
\details{
If URL link exists download Quantitative Dependency Matrix.
}
\examples{
## Downloading Quantitative Dependency Matrix
## for > 300 cancer cell lines from [1]
data(curated_BAGEL_essential)
data(curated_BAGEL_nonEssential)

DepMat<-CoRe.download_DepMatrix(ess = curated_BAGEL_essential, noness = curated_BAGEL_nonEssential)

head(DepMat)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
