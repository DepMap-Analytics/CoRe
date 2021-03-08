\name{CoRe.extract_tissueType_SubMatrix}
\alias{CoRe.extract_tissueType_SubMatrix}
\title{Subset Dependency Matrix based on tissue/cancer type of interest}
\usage{CoRe.extract_tissueType_SubMatrix(fullDepMat,
                                         tissue_type="Non-Small Cell Lung Carcinoma")}
\arguments{
\item{fullDepMat}{Dependency Matrix containing all cell models.}
\item{tissue_type}{Character, name of tissue/cancer type of interest.}
}
\value{
Dependency SubMatrix
\item{cs_depmat}{Subset of Dependency Matrix containing cell lines of interest.}
}
\description{Extracting Dependency SubMatrix for a given tissue or cancer type.}
\details{
Extracting Dependency SubMatrix for a given tissue or cancer type, among those included in the latest model annotation file on the cell model passports [1][2].
}
\examples{
## Downloading dependency matrix
## for > 300 cancer cell lines from [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Subset Binary Dependency Matrix based on "Non-Small Cell Lung Carcinoma" cancer type
LungDepMat<-CoRe.extract_tissueType_SubMatrix(BinDepMat,tissue_type="Non-Small Cell Lung Carcinoma")

head(LungDepMat)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports—a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Research 2019; 47(D1):D923–D929. doi:10.1093/nar/gky872.
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
