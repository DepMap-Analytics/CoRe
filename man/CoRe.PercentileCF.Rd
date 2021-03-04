\name{CoRe.PercentileCF}
\alias{CoRe.PercentileCF}
\title{Execute 90-th percentile method}
\usage{CoRe.PercentileCF(depMat,
                  display=TRUE,
                  percentile=0.9,
                  method='fixed',
                  thresholding='localMin')}
\arguments{
\item{depMat}{Quantitative Dependency Matrix containing Pan-cancer or tissue/cancer types specific models.}
\item{display}{Boolean, default is TRUE. Should gene score rank distribution of the dependency be plotted.}
\item{percentile}{Numerical value in range [0,1], default is 0.9. Percentile to be used as threshold.}
\item{method}{Character, default is 'fixed'. Specify which version of the 90-th percentile use, options are:\cr
  - fixed: calculates the rank distribution of genes at their least dependent 90-th percentile cell line.\cr
  - average: calculates the average rank distribution of the least dependent 90-th percentile cell lines.\cr
  - slope: employs all the cell lines to fit a linear model to generate a gene score rank distribution for each gene. The slope distribution of the genes score ranks form a bimodal distribution similar to the percentile gene score rank methods of the first two versions. Using this distribution it is possible to determine the point of minimum density between two peaks and predict the essential genes.\cr\cr
  For each version, the resulting bimodal distribution is used for the identification of essential genes under the considered cell lines.}

\item{thresholding}{Character, default is 'localMin'. Specify the thresholding strategy to adopt, options are: \cr
      - localMin: use local minimum for setting threshold. \cr
      - second option set threshold using Bayesian Factor (BF) scores. Genes having a BF score >= 10 are considered as core fitness.}
}
\value{
List of the following items:
\item{cfgenes}{A vector of string with Core Fitness Genes' symbols for the tissue/cancer type of interest.}
\item{geneRanks}{Dataframe containing rank scores for each gene.}
\item{LocalMinRank}{Discriminative threshold estimated using the 'localMin' thresholding option.}
\item{cfBFs}{Numeric vector containing Bayesian Factor scores for each gene based on gene score rank distribution.}
}
\description{Execute 90-th percentile method [1] on Quantitative Dependency Matrix.}
\details{
Calculate the Core Fitness genes using the 90th-percentile [1] least dependent cell line from Quantative knockout screen dependency matrix containing Pan-cancer or tissue/cancer types specific models.
}
\examples{
## Execute all possible combination of 90-th percentile method on Quantitative
## knockout screen dependency
data(curated_BAGEL_essential)
data(curated_BAGEL_nonEssential)

depMat<-CoRe.download_DepMatrix(scaled = TRUE, ess = curated_BAGEL_essential, noness = curated_BAGEL_nonEssential)

CFgenes<-CoRe.PercentileCF(depMat,method = 'fixed',thresholding='localMin')
CFgenesAVG<-CoRe.PercentileCF(depMat,method = 'average',thresholding='localMin')
CFgenesSLOPE<-CoRe.PercentileCF(depMat,method = 'slope',thresholding='localMin')

CFgenes_BFs<-CoRe.PercentileCF(depMat,method = 'fixed',thresholding='BFs')
CFgenesAVG_BFs<-CoRe.PercentileCF(depMat,method = 'average',thresholding='BFs')
CFgenesSLOPE_BFs<-CoRe.PercentileCF(depMat,method = 'slope',thresholding='BFs')
}
\references{
[1] Dempster, J.M., Pacini, C., Pantel, S. et al. Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets. Nat Commun 10, 5817 (2019).
}
\seealso{
\code{\link{CoRe.CalculateBayesianfactor}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
