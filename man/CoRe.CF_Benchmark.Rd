\name{CoRe.CF_Benchmark}
\alias{CoRe.CF_Benchmark}
\title{Computes recall and other ROC indicators for identified core fitness genes}
\usage{CoRe.CF_Benchmark(testedGenes,
                          background,
                          priorKnownSignatures,
                          falsePositives,
                          displayBar=FALSE)}
\arguments{
\item{testedGenes}{Vector of gene symbols that have been identified as tissue/Pan-cancer core fitness.}
\item{background}{Vector of gene symbols that are part of the Binary Dependency Matrix.}
\item{priorKnownSignatures}{List, signatures of independent known essential genes [1].}
\item{falsePositives}{Lowly expressed genes used as false positive during the benchmarking [2].}
\item{displayBar}{Boolean, default is TRUE. Should mutual exclusivity pattern of signatures be plotted.}
}
\value{
\item{TPRs}{Dataframe listing Recall and p-values (obtained from hypergeometric distribution) associated
with each signature of independent known essential genes.}
\item{PPV}{Positive predicted value obtained by comparing pooled signatures against inputted tested genes.}
\item{FPR}{False positive rate obtained by comparing pooled signatures against inputted false positive genes.}
}
\description{Computes recall and other ROC indicators for identified core fitness genes.}
\details{
Computes recall and other ROC indicators for identified core fitness genes with respect to pre-defined signatures of essential [1] and
false positive genes [2].
}
\examples{
# Benchmarking the identified PanCancer Core fitness genes against
# prior known essential genes [1]

data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.HISTONES)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)
data(EssGenes.ribosomalProteins)

data("curated_BAGEL_essential")

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

signatures<-list(DNA_REPLICATION=EssGenes.DNA_REPLICATION_cons,
                 HISTONES=EssGenes.HISTONES,
                 RNA_POLYMERASE=EssGenes.KEGG_rna_polymerase,
                 PROTEASOME=EssGenes.PROTEASOME_cons,
                 SPLICEOSOME=EssGenes.SPLICEOSOME_cons,
                 RIBOSOMAL_PROTS=EssGenes.ribosomalProteins)

clannotation<-
  CoRe.download_AnnotationModel(
  'https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz') ## dataset from [3]

PanCancer_CF_genes<-
  CoRe.PanCancer_AdAM(pancan_depMat = BinDepMat,
                      tissues_ctypes = tissues_ctypes,
                      clannotation = clannotation,
                      TruePositives = curated_BAGEL_essential,
                      display = FALSE)

FPs<-CoRe.AssembleFPs()
AdAMperf<-CoRe.CF_Benchmark(PanCancer_CF_genes,
  background = rownames(BinDepMat),priorKnownSignatures =
  signatures,falsePositives=FPs)
}
\references{
[1] Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102:15545.

[2] Barretina, J., Caponigro, G., Stransky, N. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603–607 (2012).

[3] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.
}
\seealso{
\code{\link{CoRe.AssembleFPs}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
