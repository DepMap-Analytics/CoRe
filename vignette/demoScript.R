library(CRISPRcleanR)

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
TPR<-ADAM.truePositiveRate(exampleDepMat,curated_BAGEL_essential)


############################################################################################################


# Calculate True positive rates for fitness genes in at least n cell lines in the observed dependency matrix,
# with positive cases from a reference set of essential genes
TPR<-ADAM.truePositiveRate(exampleDepMat,curated_BAGEL_essential)


# Calculate minimum number of cell lines a gene needs to be a fitness gene in order to be considered
# as a core-fitness gene
crossoverpoint<-ADAM.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'curated BAGEL essential')

#coreFitnessGenes is the list of genes predicted as core-fitness by AdAM.
coreFitnessGenes<-rownames(exampleDepMat)[rowSums(exampleDepMat)>=crossoverpoint]
