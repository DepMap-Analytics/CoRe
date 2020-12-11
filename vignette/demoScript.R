library(CRISPRcleanR)
library(CELLector)

## Downloading binary dependency matrix
## for > 300 cancer cell lines from Project Score [1]
BinDepMat<-CoRe.download_BinaryDepMatrix()

## Extracting dependency submatrix for
## Non-Small Cell Lung Carcinoma cell lines only
LungDepMap<-CoRe.extract_tissueType_BinDepMatrix(BinDepMat)

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
data(BAGEL_essential)

# Calculate True positive rates for fitness genes in at least n cell lines in the observed dependency matrix,
# with positive cases from a reference set of essential genes
TPR<-CoRe.truePositiveRate(LungDepMap,BAGEL_essential)

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

SNCLC_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Non-Small Cell Lung Carcinoma',clannotation)
BRCA_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Breast',clannotation)
CRC_cf_genes<-CoRe.CS_AdAM(BinDepMat,tissue_ctype = 'Large Intestine',clannotation)

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
PanCacer_CF_genes<-
  CoRe.PanCancer_AdAM(pancan_depMat = BinDepMat,
                      tissues_ctypes = tissues_ctypes,
                      clannotation = clannotation,
                      TruePositives = curated_BAGEL_essential,
                      display = FALSE)

RES<-systematic_CS_AdAM_res<-lapply(tissues_ctypes[8],function(x){
  cls<-clannotation$model_name[clannotation$tissue==x | clannotation$cancer_type==x]
  cls<-intersect(colnames(pancan_depMat),cls)
  cs_depmat<-pancan_depMat[,cls]

  res<-lapply(10:length(cls),function(x){
    print(x)
    do.call(rbind,lapply(1:100,function(y){
      print(y)
      ss<-sample(cls,x)
      CFs<-CoRe.AdAM(cs_depmat[,ss],display = FALSE,verbose = FALSE,
                     TruePositives = curated_BAGEL_essential)

      atLeastOne<-names(which(rowSums(cs_depmat[,ss])>0))
      CS<-length(setdiff(atLeastOne,CFs))

      return(data.frame(nCF=length(CFs),nCS=CS))
      }))
  })
  CSspec<-setdiff(names(which(rowSums(cs_depmat)>0)),CFs)
  return(list(ncls=length(cls),CFs=CFs,CSs=CSspec))
})


