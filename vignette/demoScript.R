library(CRISPRcleanR)
library(CELLector)
library(MutExMatSorting)
library(pheatmap)
library(magrittr)
library(mixdist)


## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat,tissue_type="Non-Small Cell Lung Carcinoma")

# Generate the profiles of number of fitness genes across number of cell lines from observed data and
# corresponding comulative sums.
pprofile<-CoRe.panessprofile(depMat=LungDepMap)

# Generate a set of random profiles of number of genes depleted for a number of cell lines and corresponding
# cumulative sums by perturbing observed data.
nullmodel<-CoRe.generateNullModel(depMat=LungDepMap,ntrials = 1000)

# Calculate log10 odd ratios of observed/expected profiles of cumulative number of fitness genes in fixed number of cell lines
EO<-CoRe.empiricalOdds(observedCumSum = pprofile$CUMsums,
                       simulatedCumSum =nullmodel$nullCumSUM)

#load a reference set of essential genes
data(curated_BAGEL_essential)

# Calculate True positive rates for fitness genes in at least n cell lines in the observed dependency matrix,
# with positive cases from a reference set of essential genes
TPR<-CoRe.truePositiveRate(LungDepMap,curated_BAGEL_essential)

# Calculate minimum number of cell lines a gene needs to be a fitness gene in order to be considered
# as a core-fitness gene
crossoverpoint<-CoRe.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'BAGEL essential')

#coreFitnessGenes is the list of genes predicted as core-fitness by AdAM.
coreFitnessGenes<-CoRe.coreFitnessGenes(LungDepMap,crossoverpoint)

#======================================================================
#Reperform all the analyses but with a single call to the AdAM wrapper
coreFitnessGenes<-CoRe.AdAM(LungDepMap,TruePositives=BAGEL_essential)

#======================================================================
#Reperform all the analyses but on different tissues or cancer-types
clannotation<-CELLector.CMPs_getModelAnnotation()

SNCLC_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Non-Small Cell Lung Carcinoma',
                             clannotation = clannotation,
                             TruePositives = curated_BAGEL_essential)
BRCA_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Breast',
                            clannotation = clannotation,
                            TruePositives = curated_BAGEL_essential)
CRC_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Large Intestine',
                           clannotation = clannotation,
                           TruePositives = curated_BAGEL_essential)

#======================================================================
# Identifying pan-cancer core-fitness genes with the AdAM model, as
# described in Behan et al 2019, i.e. performing analyses at individual
# tissues/cancer-type level then collapsing results at pan-cancer level

## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1]
library(CELLector)

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
  CELLector.CMPs_getModelAnnotation('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz')

data(curated_BAGEL_essential)
PanCancer_CF_genes<-
  CoRe.PanCancer_AdAM(pancan_depMat = BinDepMat,
                      tissues_ctypes = tissues_ctypes,
                      clannotation = clannotation,
                      TruePositives = curated_BAGEL_essential,
                      display = FALSE)

#======================================================================
# Benchmarking the identified PanCancer Core fitness genes against
# prior known essential genes

library(CRISPRcleanR)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.HISTONES)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)
data(EssGenes.ribosomalProteins)

data(BAGEL_essential)

signatures<-list(DNA_REPLICATION=EssGenes.DNA_REPLICATION_cons,
                 HISTONES=EssGenes.HISTONES,
                 RNA_POLYMERASE=EssGenes.KEGG_rna_polymerase,
                 PROTEASOME=EssGenes.PROTEASOME_cons,
                 SPLICEOSOME=EssGenes.SPLICEOSOME_cons,
                 RIBOSOMAL_PROTS=EssGenes.ribosomalProteins)

FPs<-CoRe.AssembleFPs()
AdAMperf<-CoRe.CF_Benchmark(PanCancer_CF_genes,background = rownames(BinDepMat),priorKnownSignatures = signatures,falsePositives=FPs)


#======================================================================
# Percentile Method

## Downloading quantitative dependency matrix
## from Project Score [1]
depMat<-CoRe.download_DepMatrix(scaled = TRUE)

CFgenes<-CoRe.PercentileCF(depMat,method = 'fixed',thresholding='localMin')
CFgenesAVG<-CoRe.PercentileCF(depMat,method = 'average',thresholding='localMin')
CFgenesSLOPE<-CoRe.PercentileCF(depMat,method = 'slope',thresholding='localMin')

CFgenes_BFs<-CoRe.PercentileCF(depMat,method = 'fixed',thresholding='BFs')
CFgenesAVG_BFs<-CoRe.PercentileCF(depMat,method = 'average',thresholding='BFs')
CFgenesSLOPE_BFs<-CoRe.PercentileCF(depMat,method = 'slope',thresholding='BFs')


HARTperf<-CoRe.CF_Benchmark(testedGenes = BAGEL_essential,
                              background = rownames(depMat),
                              priorKnownSignatures = signatures,
                              falsePositives=FPs)

Perc90perf<-CoRe.CF_Benchmark(testedGenes = CFgenes$cfgenes,
                  background = rownames(depMat),
                  priorKnownSignatures = signatures,
                  falsePositives=FPs)

Perc90AVGperf<-CoRe.CF_Benchmark(testedGenes = CFgenesAVG$cfgenes,
                              background = rownames(depMat),
                              priorKnownSignatures = signatures,
                              falsePositives=FPs)

Perc90SLOPEperf<-CoRe.CF_Benchmark(testedGenes = CFgenesSLOPE$cfgenes,
                                 background = rownames(depMat),
                                 priorKnownSignatures = signatures,
                                 falsePositives=FPs)

Perc90_BFs_perf<-CoRe.CF_Benchmark(testedGenes = CFgenes_BFs$cfgenes,
                              background = rownames(depMat),
                              priorKnownSignatures = signatures,
                              falsePositives=FPs)

Perc90AVG_BFs_perf<-CoRe.CF_Benchmark(testedGenes = CFgenesAVG_BFs$cfgenes,
                                 background = rownames(depMat),
                                 priorKnownSignatures = signatures,
                                 falsePositives=FPs)

Perc90SLOPE_BFs_perf<-CoRe.CF_Benchmark(testedGenes = CFgenesSLOPE_BFs$cfgenes,
                                   background = rownames(depMat),
                                   priorKnownSignatures = signatures,
                                   falsePositives=FPs)


barplot(rbind(length(BAGEL_essential),
              length(PanCancer_CF_genes),
              length(CFgenes$cfgenes),
              length(CFgenes_BFs$cfgenes),
              length(CFgenesAVG$cfgenes),
              length(CFgenesAVG_BFs$cfgenes),
              length(CFgenesSLOPE$cfgenes),
              length(CFgenesSLOPE_BFs$cfgenes)),
              beside = TRUE,
              ylab='n. genes',col=c('black',
                                                '#D81B60',
                                                '#1E88E5','#1EDAE5',
                                                '#FF8A07','#FFC107',
                                                '#004D40','#009C40'),border = FALSE)

barplot(rbind(HARTperf$TPRs$Recall,
  AdAMperf$TPRs$Recall,
              Perc90perf$TPRs$Recall,
              Perc90_BFs_perf$TPRs$Recall,
              Perc90AVGperf$TPRs$Recall,
              Perc90AVG_BFs_perf$TPRs$Recall,
              Perc90SLOPEperf$TPRs$Recall,
              Perc90SLOPE_BFs_perf$TPRs$Recall),beside = TRUE,ylab='Recall',
        names.arg = names(signatures),las=2,ylim=c(0,1),col=c('black','#D81B60',
                                                              '#1E88E5','#1EDAE5',
                                                              '#FF8A07','#FFC107',
                                                              '#004D40','#009C40'),border = FALSE)

barplot(rbind(1-HARTperf$PPV,
              1-AdAMperf$PPV,
              1-Perc90perf$PPV,
              1-Perc90_BFs_perf$PPV,
              1-Perc90AVGperf$PPV,
              1-Perc90AVG_BFs_perf$PPV,
              1-Perc90SLOPEperf$PPV,
              1-Perc90SLOPE_BFs_perf$PPV),beside = TRUE,
        ylab='ratio of potential novel hits (1-PPV)',
        las=2,ylim=c(0,1),col=c('black',
          '#D81B60',
                                '#1E88E5','#1EDAE5',
                                '#FF8A07','#FFC107',
                                '#004D40','#009C40'),border = FALSE)

barplot(rbind(HARTperf$FPR,
  AdAMperf$FPR,
              Perc90perf$FPR,
              Perc90_BFs_perf$FPR,
              Perc90AVGperf$FPR,
              Perc90AVG_BFs_perf$FPR,
              Perc90SLOPEperf$FPR,
              Perc90SLOPE_BFs_perf$FPR),beside = TRUE,
        ylab='not expressed genes (FPR)',
        las=2,col=c('black','#D81B60',
                    '#1E88E5','#1EDAE5',
                    '#FF8A07','#FFC107',
                    '#004D40','#009C40'),border = FALSE)

plot(rbind(1-HARTperf$PPV,
           1-AdAMperf$PPV,
            1-Perc90perf$PPV,
            1-Perc90_BFs_perf$PPV,
            1-Perc90AVGperf$PPV,
            1-Perc90AVG_BFs_perf$PPV,
            1-Perc90SLOPEperf$PPV,
            1-Perc90SLOPE_BFs_perf$PPV),
      rbind(HARTperf$FPR,
            AdAMperf$FPR,
            Perc90perf$FPR,
            Perc90_BFs_perf$FPR,
            Perc90AVGperf$FPR,
            Perc90AVG_BFs_perf$FPR,
            Perc90SLOPEperf$FPR,
            Perc90SLOPE_BFs_perf$FPR),
     xlab='n. potential novel hits (1 - PPV)',
     ylab='Recall of known core fitness genes',
     xlim=c(0,1),ylim=c(0,1))

plot(rbind(1-HARTperf$PPV,
           1-AdAMperf$PPV,
           1-Perc90perf$PPV,
           1-Perc90_BFs_perf$PPV,
           1-Perc90AVGperf$PPV,
           1-Perc90AVG_BFs_perf$PPV,
           1-Perc90SLOPEperf$PPV,
           1-Perc90SLOPE_BFs_perf$PPV),
     rbind(mean(HARTperf$TPRs$Recall),
           mean(AdAMperf$TPRs$Recall),
           mean(Perc90perf$TPRs$Recall),
           mean(Perc90_BFs_perf$TPRs$Recall),
           mean(Perc90AVGperf$TPRs$Recall),
           mean(Perc90AVG_BFs_perf$TPRs$Recall),
           mean(Perc90SLOPEperf$TPRs$Recall),
           mean(Perc90SLOPE_BFs_perf$TPRs$Recall)),
     xlab='n. potential novel hits (1 - PPV)',
     ylab='Recall of known core fitness genes',
     xlim=c(0,1),ylim=c(0,1))

abline(0,1)













