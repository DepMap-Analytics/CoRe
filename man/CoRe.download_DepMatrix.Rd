\name{CoRe.download_DepMatrix}
\alias{CoRe.download_DepMatrix}
\title{Download Quantitative Dependency Matrix}
\usage{CoRe.download_DepMatrix(
                     URL='https://cog.sanger.ac.uk/cmp/download/essentiality_matrices.zip',
                     scaled=FALSE,
                     ess=NULL,
                     noness=NULL)}
\arguments{
\item{URL}{URL pointing to the online location of a zipped folder containing a quantitative dependency matrix.
By default this will point to an entry on the data download page of the Project Score portal with data from [1].}
\item{scaled}{Boolean, default is FALSE. Should the Quantitative Dependency matrix be scaled using reference set of essential and non essential genes (provided in input), as detailed in [2].}
\item{ess}{Vector of gene symbols to be used as reference essential genes. Ignored if \code{scaled} is set to FALSE.}
\item{noness}{Vector of gene symbols to be used as reference non-essential genes. Ignored if \code{scaled} is set to FALSE.}
}
\value{
A Dependency Matrix (from [1,3]) where rows are genes and columns are cell lines. The entry in position [i,j] quantifies the effect of the inactivation of the i-th gene on the survival of the j-th cell line (its fitness).
}
\description{Downloading Quantitative Dependency Matrix introduced in Behan 2019 from Project Score [1,3].}
\details{
If the URL points to a valid online location this function downloads a quantitative Dependency Matrix.
}
\examples{

## loading reference set of essential/non-essential genes
data(curated_BAGEL_essential)
data(curated_BAGEL_nonEssential)

## Downloading unscaled Quantitative Dependency Matrix
## for > 300 cancer cell lines from [1]
DepMat<-CoRe.download_DepMatrix()

head(DepMat)

## Downloading and scaling Quantitative Dependency Matrix
DepMat_scaled<-CoRe.download_DepMatrix(scaled=TRUE, ess = curated_BAGEL_essential, noness = curated_BAGEL_nonEssential)

head(DepMat_scaled)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Meyers RM, Bryan JG, McFarland JM, et al. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nat Genet. 2017 Dec;49(12):1779-1784. doi: 10.1038/ng.3984.

[3] Dwane L, Behan FM, Gonçalves E, et al. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
