\name{CoRe.PanCancer_ADaM}
\alias{CoRe.PanCancer_ADaM}
\title{Execute ADaM at the Pan-cancer level}
\usage{CoRe.PanCancer_ADaM(pancan_depMat,
                          tissues_ctypes,
                          clannotation = NULL,
                          display=TRUE,
                          ntrials=1000,
                          verbose=TRUE,
                          TruePositives)}
\arguments{
\item{pancan_depMat}{A binary dependency matrix derived from screening (ideally 100s of) cell-lines from multiple tissue lineages and where rows are genes and columns are cell-lines/samples, with a 1 in position \emph{[i,j]} indicating that the inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th cell-line/sample.}
\item{tissues_ctypes}{Vector of strings with tissue/cancer type names of interest. These should be compatible with the cell model annotations of the Cell Model Passports [2] (downloadable through the function \code{CoRe.download_AnnotationModel}).}
\item{clannotation}{Data frame containing the Cancer cell lines' annotations, derived from the cell model passports [2] (downloadable through the function \code{CoRe.download_AnnotationModel}).}
\item{display}{Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted.}
\item{ntrials}{Integer, default =1000. How many times to randomly perturb the dependency matrix in order to generate null distributions of number of fitness genes across fixed number of cell lines.}
\item{verbose}{Boolean, default is TRUE. Should the computation progress be monitored.}
\item{TruePositives}{Vector of gene symbols to be used as prior known essential genes by the ADaM algorithm.}
}
\value{
\item{PanCancer_CF_genes}{A vector of string with predicted PanCancer Core Fitness Genes' symbols.}
}
\description{Execute ADaM at PanCancer level.}
\details{
This function executes ADaM on every tissue in cascade to identify Cancer Type specific Core Fitness genes, then iterates the procedure as detailed in [1] to identify a set of Pan-cancer core fitness genes.
}
\examples{
# Identifying pan-cancer core-fitness genes with the ADaM model, as
# described in Behan et al 2019, i.e. performing analyses at individual
# tissues/cancer-type level then iterating the proceudre at pan-cancer level

## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1,4]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Defining tissues/cancer-types that should be considered in the
## first phase of ADaM executions
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

## Downloading cell line annotations from the Cell Model Passports [2]
clannotation<-
  CoRe.download_AnnotationModel('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz') ## dataset from [2]

## Downloading a set of priori known essential genes to be used as true positives from [3] and manually
## curated as detailed in [1]
data(curated_BAGEL_essential)

## Execute ADaM at the pancancer level
PanCancer_CF_genes<-
  CoRe.PanCancer_ADaM(pancan_depMat = BinDepMat,
                      tissues_ctypes = tissues_ctypes,
                      clannotation = clannotation,
                      TruePositives = curated_BAGEL_essential,
                      display = FALSE)


## Inspect output
PanCancer_CF_genes
}
\references{
[1] Behan FM, Iorio F, Picco G, Gonçalves E, Beaver CM, Migliardi G, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019;568:511–6.

[2] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.

[3] Hart T, Chandrashekhar M, Aregger M, Steinhart Z, Brown KR, MacLeod G, Mis M, Zimmermann M, Fradet-Turcotte A, Sun S, Mero P, Dirks P, Sidhu S, Roth FP, Rissland OS, Durocher D, Angers S, Moffat J. High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities. Cell. 2015 Dec 3;163(6):1515-26. doi: 10.1016/j.cell.2015.11.015. Epub 2015 Nov 25. PMID: 26627737.

[4] Dwane L, Behan FM, Gonçalves E, Lightfoot H, Yang W, van der Meer D, Shepherd R, Pignatelli M, Iorio F, Garnett MJ. Project Score database: a resource for investigating cancer cell dependencies and prioritizing therapeutic targets. Nucleic Acids Res. 2021 Jan 8;49(D1):D1365-D1372.

}
\seealso{
\code{\link{CoRe.CS_ADaM}}
\code{\link{CoRe.ADaM}}
\code{\link{CoRe.download_AnnotationModel}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
