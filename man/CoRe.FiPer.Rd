\name{CoRe.FiPer}
\alias{CoRe.FiPer}
\title{Fitness Percentile method to estimate common essential genes}
\usage{CoRe.FiPer(depMat,
                  display=TRUE,
                  percentile=0.9,
                  method='fixed')}
\arguments{
\item{depMat}{Quantitative Dependency Matrix containing Pan-cancer or tissue/cancer-types specific gene fitness/dependency scores across cell-lines/samples. The value in position \eqn{[i,j]} of such matrix quantifies the fitness/dependency score of the \eqn{i}-th gene in the \eqn{j}-th cell line.}
\item{display}{Boolean, default is TRUE. Should gene score rank distributions of dependency scores be plotted.}
\item{percentile}{Numerical value in range [0,1], default is 0.9. Percentile to be used as a threshold.}
\item{method}{Character, default is 'fixed'. This parameter specifies which variant of the Fitness Percentile method should be used. Admissimble values are:\cr
  - fixed: a distribution of gene fitness-rank-positions in their most dependent \code{n}-th (determined by the \code{percentile} parameter) percentile cell line is used in the subsequent step of the Fitness Percentile method.\cr
  - average: a distribution of gene average fitness-rank-position across all cell lines at or over the \code{n}-th percentile of most dependent cell lines (where \eqn{n} is determined by the \code{percentile} parameter) is used in the subsequent step of the Fitness Percentile method.\cr
  - slope: for each gene, a linear model is fit on the sequence of gene fitness-rank-positions across all cell lines sorted according to their dependency on that gene, then a distrubution of models' slopes is used in the subsequent step of the Fitness Percentile method.\cr
  - AUC: for each gene, the area under the curve resulting from considering the sequence of gene fitness-rank-positions across all cell lines sorted according to their dependency on that gene is used in the subsequent step of the Fitness Percentile method.\cr
  \cr}
}
\value{
List of the following items:
\item{cfgenes}{A vector of strings with common essential genes' symbols for the tissue/cancer type of interest.}
\item{geneRanks}{Dataframe containing rank scores for each gene.}
\item{LocalMinRank}{Discriminative rank threshold.}
}
\description{Fitness percentile method to identify common essential genes from the joint analysis of multiple gene dependency profiles.}
\details{
This function implements the Fitness Percentile method to estimate common essential genes from multiple gene dependency screens introduced in [1]. For each gene in the quantitative dependency matrix provided in input a score is computed using three possible variants of the method. This can ben its fitness/essentiality rank when considering all gene essentiality scores in its \code{n}-th percentile least dependent cell line, or the average rank when considering all cell lines falling within or over the
\code{n}-th percentile of least dependent cell lines, or the slope of the curve obtained when fitting a linear model on the sequence of fitness-rank-positions across cell lines sorted according to the dependency on that gene (where both variant and \code{n} are user defined.)

The density of these gene scores (generally bimodal) is then estimated using a Gaussian kernel and
the central point of minimum density is identified. Genes whose score falls
below the point of minimum density are classified as common essential.
This method is fully detailed in [1]. Observed histrogram, estimated density distributions and discriminative threshold are also plotted.
}
\examples{
## Downloading a set of priori known essential genes to be used as true positives from [2] and manually
## curated as detailed in [3]
data(curated_BAGEL_essential)
data(curated_BAGEL_nonEssential)

## Downloading and scaling quantitative dependency matrix from project score [3,4]
depMat<-CoRe.download_DepMatrix(scaled = TRUE, ess = curated_BAGEL_essential, noness = curated_BAGEL_nonEssential)

## Executing the three variants of the Fitness percentile method
CFgenes<-CoRe.FiPer(depMat,method = 'fixed')
CFgenesAVG<-CoRe.FiPer(depMat,method = 'average')
CFgenesSLOPE<-CoRe.FiPer(depMat,method = 'slope')
CFgenesAUC<-CoRe.FiPer(depMat,method = 'AUC')

## Inspect the identified common essential genes
CFgenes$cfgenes
CFgenesAVG$cfgenes
CFgenesSLOPE$cfgenes
CFgenesAUC$cfgenes
}
\references{
[1] Dempster, J.M., Pacini, C., Pantel, S. et al. Agreement between two large pan-cancer CRISPR-Cas9 gene dependency data sets. Nat Commun 10, 5817 (2019).

[2] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, Mis M, Zimmermann M, Fradet-Turcotte A, Sun S, Mero P, Dirks P, Sidhu S, Roth FP, Rissland OS, Durocher D, Angers S, Moffat J. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26. doi: 10.1016/j.cell.2015.11.015. Epub 2015 Nov 25. PMID: 26627737.

[3] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[4] Dwane L, Behan FM, Gonçalves E, Lightfoot H, Yang W, van der Meer D, Shepherd R, Pignatelli M, Iorio F, Garnett MJ. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.

}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
