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
\item{tissue_ctype}{A string specifying the tissue/cancer type of interest, this must be compliant with the Cell Model Passports annotation [1].}
\item{clannotation}{Cancer cell line models' annotation from the cell model passports. This can be downloaded using the \code{\link{CoRe.download_AnnotationModel}} function}
\item{display}{Boolean, default is TRUE. Should bar plots of dependency profiles and boxplots of estimated empirical distribution be visualised.}
\item{main_suffix}{If display=TRUE, title suffix to be given to plots of number of genes that are essential/fitness in a give number of cell lines, default is 'genes depleted in at least 1 cell line'.}
\item{xlab}{x-axis label of the plots, default is 'n. cell lines'.}
\item{ntrials}{Integer, default =1000. How many times the dependency matrix shouldd be suffled in order to generate null distributions of number of genes that are essential in fixed numbers of cell lines}
\item{verbose}{Boolean, default is TRUE. Should the computation progress be monitored.}
\item{TruePositives}{Vector of gene symbols to be used as prior known essential genes.}
}
\value{
\item{coreFitnessGenes}{A vector of strings with estimated Core Fitness Genes' symbols for the tissue/cancer type of interest.}
}
\description{Execute AdAM on a tissue- or cancer-type-specific binary dependency submatrix.}
\details{
Execute sequentially the whole AdAM pipeline on a tissue or cancer type specific dependency submatrix.
}
\examples{
## downloading a reference set of prior known essential genes from [3]
## curated as detailed in [2]

data(curated_BAGEL_essential)

## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [2]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Perform all the analyses but on different tissues or cancer-types
clannotation<-
  CoRe.download_AnnotationModel('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz') ## dataset from [2]

SNCLC_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Non-Small Cell Lung Carcinoma',
                             clannotation = clannotation,
                             TruePositives = curated_BAGEL_essential)
}
\references{
[1] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.

[2] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[3] Hart T, Chandrashekhar M, Aregger M, et al. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26.
}
\seealso{
\code{\link{CoRe.AdAM}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
