\name{CoRe.CS_AdAM}
\alias{CoRe.CS_AdAM}
\title{Execute AdAM on a specific tissue/cancer type}
\usage{CoRe.CS_AdAM(pancan_depMat,
                     tissue_ctype = 'Non-Small Cell Lung Carcinoma',
                     clannotation = NULL,
                     display=TRUE,
                     main_suffix='fitness genes in at least 1 cell line',
                     xlab='n. dependent cell lines',
                     ntrials=1000,
                     verbose=TRUE,
                     TruePositives)}
\arguments{
\item{pancan_depMat}{Binary Dependency Matrix containing all cell models.}
\item{tissue_ctype}{Character, name of tissue/cancer type of interest.}
\item{clannotation}{Cancer cell lines latest model annotation file on the cell model passports.}
\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted.}
\item{main_suffix}{If display=TRUE, title suffix to give to plot of number of genes depleted in a give number of cell lines, default is 'genes depleted in at least 1 cell line'.}
\item{xlab}{label to give to x-axis of the plots, default is 'n. cell lines'.}
\item{ntrials}{Integer, default =1000. How many times to randomly perturb dependency matrix to generate the null distributions.}
\item{verbose}{Boolean, default is TRUE. Should the computation progress be monitored.}
\item{TruePositives}{Vector of gene symbols to be used as prior known essential genes.}
}
\value{
\item{coreFitnessGenes}{A vector of string with Core Fitness Genes' symbols for the tissue/cancer type of interest.}
}
\description{Execute AdAM on tissue or cancer type specific dependency submatrix.}
\details{
Execute sequentially the whole AdAM pipeline on tissue or cancer type specific dependency submatrix.
}
\examples{
## Downloading Quantitative Dependency Matrix
## for > 300 cancer cell lines from [1]

data(curated_BAGEL_essential)

## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Perform all the analyses but on different tissues or cancer-types
clannotation<-
  CoRe.download_AnnotationModel('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz') ## dataset from [2]

SNCLC_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Non-Small Cell Lung Carcinoma',
                             clannotation = clannotation,
                             TruePositives = curated_BAGEL_essential)
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.
}
\seealso{
\code{\link{CoRe.AdAM}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
