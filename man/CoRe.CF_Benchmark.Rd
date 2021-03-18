\name{CoRe.CF_Benchmark}
\alias{CoRe.CF_Benchmark}
\title{Recall of known essential genes and ROC indicators}
\usage{CoRe.CF_Benchmark(testedGenes,
                          background,
                          priorKnownSignatures,
                          falsePositives,
                          displayBar=FALSE)}
\arguments{
\item{testedGenes}{Vector of gene symbols that have been identified as tissue-specific or Pan-cancer core fitness genes.}
\item{background}{Vector of gene symbols included in the Dependency Matrix used to make the prediction (the background population).}
\item{priorKnownSignatures}{A List of string vectors containg each a signature of prior known essential genes (their symbol)[1].}
\item{falsePositives}{Genes to be used to compute false positive rates, this can be for example lowly expressed genes from the CCLE [2], assembled through the \code{CoRe.AssembleFPs} function.}
\item{displayBar}{Boolean, default is TRUE. Should a heatmap of the signatures' coverage be plotted.}
}
\value{
\item{TPRs}{Dataframe listing Recall and enrichment p-values (obtained from hypergeometric distribution) associated with each signature of prior known essential genes.}
\item{PPV}{Positive predicted value obtained pooling all inputed signatures together and using them as positive cases.}
\item{FPR}{False positive rate of the inputed false positive genes.}
}
\description{This function assesses the set of predicted core fitness genes by computing the recall (and other ROC indicators) of prior known essential genes and false positives.}
\details{
Computes recall and other ROC indicators for identified core fitness genes with respect to pre-defined signatures of essential and false positive genes defined in input.
}
\examples{
# Benchmarking the identified PanCancer Core fitness genes against
# prior known essential genes [1]

# loading signtures of prior known essential genes
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.HISTONES)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)
data(EssGenes.ribosomalProteins)
data(curated_BAGEL_essential)

signatures<-list(DNA_REPLICATION=EssGenes.DNA_REPLICATION_cons,
                 HISTONES=EssGenes.HISTONES,
                 RNA_POLYMERASE=EssGenes.KEGG_rna_polymerase,
                 PROTEASOME=EssGenes.PROTEASOME_cons,
                 SPLICEOSOME=EssGenes.SPLICEOSOME_cons,
                 RIBOSOMAL_PROTS=EssGenes.ribosomalProteins)

# downloading binary dependency matrix from project Score [3]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Running ADaM [3] to identify Pan-Cancer core fitness genes

## defining the cell line tissues to be used in the first step of ADaM
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

## Downloading cell line model annotations from the Cell Model Passports [3]
clannotation<-
  CoRe.download_AnnotationModel(
  'https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz')

## Running ADaM [2]
PanCancer_CF_genes<-
  CoRe.PanCancer_ADaM(pancan_depMat = BinDepMat,
                      tissues_ctypes = tissues_ctypes,
                      clannotation = clannotation,
                      TruePositives = curated_BAGEL_essential,
                      display = FALSE)

## Assemling lowly expressed genes from the CCLE [2]
FPs<-CoRe.AssembleFPs()

## benchmarking the core fitness genes predicted by ADaM
## plotting a heatmap highlighting the recalled prior known essential genes
## with barplots and enrichhment pvalues
ADaMperf<-CoRe.CF_Benchmark(PanCancer_CF_genes,
  background = rownames(BinDepMat),priorKnownSignatures =
  signatures,falsePositives=FPs)

## Inspecting TPRs, PPV and FPR
ADaMper$

}
\references{
[1] Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102:15545.

[2] Barretina, J., Caponigro, G., Stransky, N. et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603–607 (2012).

[3] Behan FM, Iorio F, Picco G, Gonçalves E, et al. Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens. Nature. 2019 Apr;568(7753):511-516.

[4] Van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929.
}
\seealso{
\code{\link{CoRe.AssembleFPs}}
}
\author{
C. Pacini, E. Karakoc, A. Vinceti & F. Iorio
}
\keyword{functions}
