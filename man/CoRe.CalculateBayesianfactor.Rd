\name{CoRe.CalculateBayesianfactor}
\alias{CoRe.CalculateBayesianfactor}
\title{Compute Bayesian Factors}
\usage{
CoRe.CalculateBayesianfactor(RankDistribution,
                             display=TRUE)
}
\arguments{
  \item{RankDistribution}{Dataframe of gene scores rank distribution.}
  \item{display}{Boolean, default is TRUE. Should probability density distribution derived from the normal mixture model be plotted.}
}
\description{
Compute Bayesian Factors based on gene score rank distribution.
}
\details{
Compute Bayesian Factors based on gene score rank distribution using a normal mixture model.
}
\value{
  \item{bak}{Numeric vector of Bayesian Factors derived from mixture model.}
}
\examples{
## Compute Bayesian Factor scores
data(curated_BAGEL_essential)
data(curated_BAGEL_nonEssential)

depMat<-CoRe.download_DepMatrix(scaled = TRUE, ess = curated_BAGEL_essential,
  noness = curated_BAGEL_nonEssential) ## Quantitative Dependency Matrix introduced in [1]

RankDistribution <- data.frame(rankings=sample(1:nrow(depMat),
  size = nrow(depMat), replace = TRUE))
rownames(RankDistribution) <- rownames(depMat)

cfBFs<-CoRe.CalculateBayesianfactor(RankDistribution,display=FALSE)
names(cfBFs) <- rownames(RankDistribution)
head(cfBFs)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.
}
\seealso{
\code{\link{CoRe.PercentileCF}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
