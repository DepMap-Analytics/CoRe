\name{CoRe.PanCancer_AdAM}
\alias{CoRe.PanCancer_AdAM}
\title{Execute AdAM at the Pan-cancer level}
\usage{CoRe.PanCancer_AdAM(pancan_depMat,
                          tissues_ctypes,
                          clannotation = NULL,
                          display=TRUE,
                          ntrials=1000,
                          verbose=TRUE,
                          TruePositives)}
\arguments{
\item{pancan_depMat}{Binary Dependency Matrix containing all cell models.}
\item{tissues_ctypes}{Vector of string with tissue/cancer type names of interest.}
\item{clannotation}{Cancer cell lines latest model annotation file on the cell model passports.}
\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted.}
\item{ntrials}{Integer, default =1000. How many times to randomly perturb dependency matrix to generate the null distributions.}
\item{verbose}{Boolean, default is TRUE. Should the computation progress be monitored.}
\item{TruePositives}{Vector of gene symbols to be used as prior known essential genes.}
}
\value{
\item{PanCancer_CF_genes}{A vector of string with Core Fitness Genes' symbols for the tissue/cancer type of interest.}
}
\description{Execute AdAM on tissue or cancer type specific dependency submatrix.}
\details{
Execute sequentially the whole AdAM pipeline on every tissue and identify Pan-cancer core fitness genes.
}
\examples{
# Identifying pan-cancer core-fitness genes with the AdAM model, as
# described in Behan et al 2019, i.e. performing analyses at individual
# tissues/cancer-type level then collapsing results at pan-cancer level

## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

tissues_ctypes<-c("Haematopoietic and Lymphoid",
                  "Ovary",
                  "Peripheral Nervous System",
                  "Central Nervous System",
                  "Pancreas",
                  "Head and Neck",
                  "Bone",
                  "Lung",
                  "Large Intestine",
                  "Esophagus",
                  "Endometrium",
                  "Stomach",
                  "Breast")

clannotation<-
  CoRe.download_AnnotationModel('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz') ## dataset from [2]

data(curated_BAGEL_essential)
PanCancer_CF_genes<-
  CoRe.PanCancer_AdAM(pancan_depMat = BinDepMat,
                      tissues_ctypes = tissues_ctypes,
                      clannotation = clannotation,
                      TruePositives = curated_BAGEL_essential,
                      display = FALSE)

}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.
}
\seealso{
\code{\link{CoRe.CS_AdAM}}
\code{\link{CoRe.AdAM}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
