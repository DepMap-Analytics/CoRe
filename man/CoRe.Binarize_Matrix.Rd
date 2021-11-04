\name{CoRe.Binarize_Matrix}
\alias{CoRe.Binarize_Matrix}
\title{Binarize Quantitative Dependency Matrix}
\usage{CoRe.Binarize_Matrix(depMat,
                            method=c('fdr', 'thr'),
                            ess_genes=NULL,
                            noness_genes=NULL,
                            scaled=FALSE,
                            Bayes_Factor=FALSE,
                            FDRth=0.05,
                            threshold=NULL)}
\arguments{
\item{depMat}{Quantitative Dependency Matrix containing Pan-cancer or tissue/cancer-types specific gene fitness/dependency scores across cell-lines/samples. The value in position \eqn{[i,j]} of such matrix quantifies the fitness/dependency score of the \eqn{i}-th gene in the \eqn{j}-th cell line.}
\item{method}{A character string indicating which method is to be used for the computation of the binary matrix. One of "fdr" or "thr".\cr
  - fdr: For every cell line, the function computes the threshold at fixed FDR computed on the two distributions derived from sets of reference essential and non-essential genes. Genes with a logFC score less than the FDR threshold will be converted to 1, 0 otherwise. Different percentages of FDR can be specified by the user, according to the level of stringency (default is 0.05).\cr
  - thr: For every cell line, any gene with a logFC score less than \code{threshold} will be converted to 1, 0 otherwise.\cr}
\item{ess_genes}{Vector of gene symbols to be used as reference essential genes. Ignored if \code{scaled} is set to FALSE.}
\item{noness_genes}{Vector of gene symbols to be used as reference non-essential genes. Ignored if \code{scaled} is set to FALSE.}
\item{scaled}{Boolean, default is FALSE. Should the Quantitative Dependency matrix be scaled using reference set of essential and non essential genes (provided in input), as detailed in [2]. We always recommend scaling as a normalization procedure.}
\item{Bayes_Factor}{Boolean, default is FALSE. It indicates whether the Quantitative Dependency matrix was obtained from BAGEL2 [3].}
\item{FDRth}{a numerical value > 0 and < 1 specifying the false discovery rate threshold at which the matrix binarization will be computed, using two reference sets of essential and non-essential genes.}
\item{threshold}{default is NULL. A numerical value specifying a hard logFC score threshold. For every cell line, any gene with a logFC score less than \code{threshold} will be converted to 1, 0 otherwise. Scaling is warmly recommended for this method, in order to normalize essentiality profiles across cell lines.}
}
\value{
A binarized version of the Dependency Matrix where rows are genes and columns are cell lines. The entry in position [i,j] indicates whether the i-th gene is essential for the viability of the j-th cell line.
}
\description{Binarazation of Quantitative Dependency Matrix.}
\details{
This function implements several strategies for the binarization of a quantitative Dependency Matrix.
}
\examples{

## loading reference set of essential/non-essential genes
data(curated_BAGEL_essential)
data(curated_BAGEL_nonEssential)

## Downloading and scaling Quantitative Dependency Matrix
## for > 300 cancer cell lines from [1]
depMat<-CoRe.download_DepMatrix()

## Binarization of logFC matrix using the FDR method, and
## scaling using sets of reference essential/non-essential genes
BinDepMat<-CoRe.Binarize_Matrix(depMat,
                                method='fdr',
                                ess_genes=curated_BAGEL_essential,
                                noness_genes=curated_BAGEL_nonEssential,
                                scaled=TRUE)
head(BinDepMat)

## Binarization of logFC matrix using the threshold method, and
## scaling using sets of reference essential/non-essential genes
BinDepMat<-CoRe.Binarize_Matrix(depMat,
                                method='thr',
                                ess_genes=curated_BAGEL_essential,
                                noness_genes=curated_BAGEL_nonEssential,
                                scaled=TRUE,
                                threshold = -0.5)
head(BinDepMat)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Meyers RM, Bryan JG, McFarland JM, et al. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nat Genet. 2017 Dec;49(12):1779-1784. doi: 10.1038/ng.3984.

[3] Eiru Kim & Traver Hart. Improved analysis of CRISPR fitness screens and reduced off-target effects with the BAGEL2 gene essentiality classifier. Genome Medicine 13, 2 (2021). doi: 10.1186/s13073-020-00809-3.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
