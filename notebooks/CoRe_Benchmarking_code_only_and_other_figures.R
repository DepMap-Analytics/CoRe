## Working directory needs to be /CoRe/notebooks

options(warn=-1)

library(tidyverse)
library(pheatmap)
library(CoRe)
library(magrittr)
library(nVennR)
library(limma)


data('curated_BAGEL_essential')
data('curated_BAGEL_nonEssential')

url <- 'https://www.depmap.org/broad-sanger/integrated_Sanger_Broad_essentiality_matrices_20201201.zip'
temp <- tempfile()
download.file(url, temp, mode="wb")
unzip(temp, exdir = 'integrated_dataset')
unlink(temp)

depFC <- read.table('integrated_dataset/integrated_Sanger_Broad_essentiality_matrices_20201201/CERES_FC.txt',
                    row.names = 1, sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
system('rm -r integrated_dataset')

print('CERES processed joint dataset of cancer dependencies correctly downloaded')


clannotation<-CoRe.download_AnnotationModel(URL='https://cog.sanger.ac.uk/cmp/download/model_list_20210326.csv')

## convert model IDs to model names
depFC <- depFC[,-which(is.na(colSums(depFC)))]
cells <- colnames(depFC)

dim(depFC)

tot_id <- clannotation %>% filter(model_id %in% cells | BROAD_ID %in% cells) %>%
  select(model_name,model_id,BROAD_ID) %>%
  pivot_longer(cols = c("model_id","BROAD_ID"), names_to = "institute", values_to = "model_id") %>%
  select(model_name,model_id) %>% filter(model_id %in% cells) %>%
  arrange(factor(model_id, levels = cells))

depFC <- depFC[,which(cells %in% tot_id$model_id)]
colnames(depFC) <- tot_id$model_name
scaled_depFC <- CoRe.scale_to_essentials(depFC,curated_BAGEL_essential,curated_BAGEL_nonEssential)

bdep <- apply(scaled_depFC, 2, function(x){
  x[which(x >= -0.5)] <- 0
  x[which(x < -0.5)] <- 1
  x
})

print('Done: quantitative and binary dependency matrices have been created')
print(paste('(accounting for',dim(bdep)[1],'genes and',dim(bdep)[2],'cell lines'))


## determining tissues with at least 15 cell lines in the dependency matrix
tissues_ctypes<-
  names(which(summary(as.factor(clannotation$tissue[match(colnames(bdep),clannotation$model_name)]))>=15))

print('ADaM will be executed at the tissue type level on the following tissue lineages:')
print(tissues_ctypes)

## The procomputed ADaM CFGs are also available and can be loaded uncommenting
## and executing this line instead of the the following one:

load('data/preComputed/ADaM.RData')

#ADaM <- CoRe.PanCancer_ADaM(bdep,
#                            tissues_ctypes,
#                            clannotation = clannotation,
#                            display=FALSE,
#                            ntrials=1000,
#                            verbose=TRUE,
#                            TruePositives = curated_BAGEL_essential)

print('Done')

load('data/preComputed/FiPer_outputs.RData')

## FiPer with fixed strategy
#Perc_fixed <- CoRe.FiPer(scaled_depFC,
#                         display=FALSE,
#                         method = 'fixed')$cfgenes

## FiPer with average strategy
#Perc_avg <- CoRe.FiPer(scaled_depFC,
#                       display=FALSE,
#                       method = 'average')$cfgenes

## FiPer with slope strategy
#Perc_slope <- CoRe.FiPer(scaled_depFC,
#                         display=FALSE,
#                         method = 'slope')$cfgenes

## Fiper with AUC strategy
#Perc_AUC <- CoRe.FiPer(scaled_depFC,
#                       display=FALSE,
#                       method = 'AUC')$cfgenes

## consensual core-fitness genes across the three less stringent variants
#Perc_Consensus <- Reduce(intersect,list(Perc_fixed,Perc_slope,Perc_AUC))

#print('Done')

## Hart2014 (built-in the CoRe package)
data(BAGEL_essential)
Hart_2014<-BAGEL_essential

## Hart2017
load('data/BAGEL_v2_Essentials.RData')
Hart_2017<-BAGEL_essential

## Behan2019
load('data/ADaM_CFs_Behan_et_Al_2019.RData')
Behan_2019<-PanCancerCoreFitnessGenes

## Sharma2020
Sharma_2020<-read.table(file = 'data/CenTools_essential_Sharma_et_al.txt',sep='\t',stringsAsFactors = FALSE)$V1

print(paste('Loaded',length(Hart_2014),'CFGs from Hart2014'))
print(paste('Loaded',length(Hart_2017),'CFGs from Hart2017'))
print(paste('Loaded',length(Behan_2019),'CFGs from Behan2019'))
print(paste('Loaded',length(Sharma_2020),'CFGs from Sharma2020'))

## Executing a logistic regression model part of the CEN-tools suite  using curated BAGEL essential
## and curated BAGEL never-essential genes
## This is followed by a k-means clustering and a prediction phase based on the clusters' silhouettes
## identifying core essential genes.

load('data/preComputed/CENtools.RData')

# write.csv(scaled_depFC, 'CENtools/data/curated_data/CERES_scaled_depFC.csv', quote = FALSE)
#
# system('pip3 install -r CENtools/requirements.txt')
#
# system('mkdir -p CENtools/data/objects')
# system('python3 CENtools/LR.py', wait = TRUE)
#
# source('CENtools/clustering.R')
#
# project_path = 'CENtools/prediction_output/INTEGRATED/'
# CENtools <- ClusterEssentiality(Chosen_project= 'INTEGRATED',
#                                 binPath= paste0(project_path, 'INTEGRATED_Histogram_DF_20_BIN.txt'),
#                                 resultPath = project_path)
#
# system('rm -r CENtools/prediction_output/')
# system('rm CENtools/data/curated_data/CERES_scaled_depFC.csv')

print(paste('Computed',length(CENtools),'CFGs via the logistic regression based method (part of CEN-tools)'))

## Assembling all CFGs in a unique list and creating a vector of colors to be used in the plots
CFs_sets<-list(Hart_2014,
               Hart_2017,
               Behan_2019,
               Sharma_2020,
               CENtools,
               ADaM,
               Perc_avg,
               Perc_Consensus,
               Perc_slope,
               Perc_AUC,
               Perc_fixed)

names(CFs_sets)<-c('Hart 2014',
                   'Hart 2017',
                   'ADaM (Behan 2019)',
                   'CEN-tools (Sharma 2020)',
                   'CEN-tools',
                   'ADaM',
                   'FiPer Average',
                   'FiPer Consensual',
                   'FiPer Slope',
                   'FiPer AUC',
                   'FiPer Fixed')

col=c("#03B2C8","#034DD9","#E6AB02","#1B9E77","#337100","#FC8D62",
      '#F4CAE4','#800080','#CAB2D6','#BEBADA','#777892')
names(col)<-names(CFs_sets)

## Assembling a set of genes included in one of the sets used as training set by at least one method
TrainingSets<-unique(c(BAGEL_essential,
                       BAGEL_nonEssential,
                       curated_BAGEL_essential,
                       curated_BAGEL_nonEssential))

print(paste(length(TrainingSets),'genes used for training by at least one method'))

## Creating a list of putative novel CFGs by removing the training genes from the sets of predicted CFGs
novelCFs_sets<-lapply(CFs_sets,function(x){
  setdiff(x,TrainingSets)
})

## Comparing and plotting CFG sets' sizes.
GeneLengths<-unlist(lapply(CFs_sets,length))
novelGeneLengths<-unlist(lapply(novelCFs_sets,length))

print('novel hits:')
print(sort(novelGeneLengths,decreasing=TRUE))

print('Unsupervised method all hits:')
print(GeneLengths[c('FiPer Average','FiPer Slope','FiPer AUC','FiPer Fixed')])
print(median(GeneLengths[c('FiPer Average','FiPer Slope','FiPer AUC','FiPer Fixed')]))

print('Unsupervised method novel hits:')
print(novelGeneLengths[c('FiPer Average','FiPer Slope','FiPer AUC','FiPer Fixed')])
print(median(novelGeneLengths[c('FiPer Average','FiPer Slope','FiPer AUC','FiPer Fixed')]))

print(paste('FiPer consensual n.genes (total): ',length(CFs_sets$`FiPer Consensual`)))
print(paste('FiPer consensual n.genes (novel hits):',length(novelCFs_sets$`FiPer Consensual`)))

par(mar=c(4,12,2,2))
barplot(rbind(novelGeneLengths,GeneLengths-novelGeneLengths),
        las=2,border = FALSE,horiz = TRUE,xlim=c(0,2300),
        xlab='n. genes')

legend('bottomright',legend = c('Training sets','(Hart2017 + curated Hart2014)'),
       fill=c('lightgray',NA),border=NA,bty = 'n')

myNV <- plotVenn(CFs_sets[c(7,9:11)],outFile='SuppFig1A.svg')
myNV <- plotVenn(novelCFs_sets[c(7,9:11)],outFile='SuppFig1B.svg')

par(mar=c(4,12,2,2))
barplot(rbind(novelGeneLengths,GeneLengths-novelGeneLengths),las=2,border = FALSE,horiz = TRUE,xlim=c(0,2300))


## Computing Recall rates of state-of-the-art sets of CFGs.

Recall_Hart2014<-unlist(lapply(CFs_sets,function(x){
  100*length(intersect(x,CFs_sets$`Hart 2014`))/length(CFs_sets$`Hart 2014`)
}))

Recall_Hart2017<-unlist(lapply(CFs_sets,function(x){
  100*length(intersect(x,CFs_sets$`Hart 2017`))/length(CFs_sets$`Hart 2017`)
}))

Recall_Behan2019<-unlist(lapply(CFs_sets,function(x){
  100*length(intersect(x,CFs_sets$`ADaM (Behan 2019)`))/length(CFs_sets$`ADaM (Behan 2019)`)
}))

Recall_Sharma2020<-unlist(lapply(CFs_sets,function(x){
  100*length(intersect(x,CFs_sets$`CEN-tools (Sharma 2020)`))/length(CFs_sets$`CEN-tools (Sharma 2020)`)
}))


Recalls<-rbind(Recall_Hart2014[6:11],Recall_Hart2017[6:11],Recall_Behan2019[6:11],Recall_Sharma2020[6:11])
rownames(Recalls)<-c('Hart2014','Hart2017','Behan2019','Sharma2020')
## creating some space for displaying the legend
Recalls<-rbind(Recalls,rep(NA,6))
Recalls<-rbind(Recalls,rep(NA,6))

par(mar=c(6,4,2,1))
barplot(t(Recalls),beside = TRUE,col=col[names(CFs_sets)[6:11]],ylim=c(0,100),ylab='%',
        border=FALSE,main='Recall of prior sets of CFGs',las=2)
legend('right',names(CFs_sets)[6:11],cex=0.9,fill=col[names(CFs_sets)[6:11]],border=NA,bty = 'n')

print(paste('ADaM median Recall across prior known sets:',median(Recalls[,'ADaM'],na.rm = TRUE)))
print(paste('FiPer median Recall across prior known sets avg across variants:',
            mean(apply(Recalls[,2:6],MARGIN = 2,median,na.rm = TRUE))))

## Adding training sets to Sharma2020 and CENtools predictions to explore overall sets similarities
CFs_sets_plus_training<-CFs_sets
CFs_sets_plus_training$`CEN-tools (Sharma 2020)`<-
  union(CFs_sets_plus_training$`CEN-tools (Sharma 2020)`,BAGEL_essential)

CFs_sets_plus_training$`CEN-tools`<-
  union(CFs_sets_plus_training$`CEN-tools`,curated_BAGEL_essential)

##### CF gene set similarity
allEss<-unique(unlist(CFs_sets_plus_training))
membMat<-do.call(rbind,lapply(CFs_sets_plus_training,function(x){is.element(allEss,x)}))
rownames(membMat)<-names(CFs_sets)
colnames(membMat)<-allEss
membMat<-membMat+0

rownames(membMat)[which(rownames(membMat)=="CEN-tools (Sharma 2020)")]<-"CEN-tools (Sharma 2020) + Hart 2017"
rownames(membMat)[which(rownames(membMat)=="CEN-tools")]<-"CEN-tools (Sharma 2020) + curated Hart 2014"

pheatmap(membMat,show_colnames = FALSE,clustering_distance_rows = 'binary',cluster_cols = FALSE,border_color = NA,
         main = 'Similarity across sets',col=c('white','blue'),legend_labels = c('out','in'),legend_breaks = c(0,1))

dmat<-dist(membMat,method = 'binary')

pheatmap(1-as.matrix(dmat), main = 'JS for Pan-cancer core fitness genes across methods',
         legend = FALSE,display_numbers = round(dmat,digits = 2))



## Assembling a set of prior known essential genes that are not included in the
## training sets used by CENtools and/or ADaM

## Loading built sets of essential genes from Iorio et al, 2018
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)
data(EssGenes.ribosomalProteins)
data(EssGenes.HISTONES)

## Adding additional signatures of known CFGs from Pacini et al, 2020
load("data/Kegg.DNArep.Rdata")
load("data/Kegg.Ribosome.Rdata")
load("data/Kegg.Proteasome.Rdata")
load("data/Kegg.Spliceosome.Rdata")
load("data/Kegg.RNApoly.Rdata")
load("data/Histones.Rdata")

positiveControls<-c(EssGenes.DNA_REPLICATION_cons,
  EssGenes.KEGG_rna_polymerase,
  EssGenes.PROTEASOME_cons,
  EssGenes.SPLICEOSOME_cons,
  EssGenes.ribosomalProteins,
  EssGenes.HISTONES,
  Kegg.DNArep,
  Kegg.Ribosome,
  Kegg.Proteasome,
  Kegg.Spliceosome,
  Kegg.RNApoly,
  Histones)

positiveControls_includingTraining<-positiveControls

## removing training sets
positiveControls<-setdiff(positiveControls,TrainingSets)

print(paste(length(positiveControls_includingTraining),'genes in the independent reference sets'))
print(paste(length(positiveControls),'genes not included in any of the training sets will be used as independent positive controls'))

## Assembling a set of genes not expressed in human cancer cell lines (FPKM < 0.1 in more than 1,000 cell lines from
## the Cell Model Passports [10], https://cellmodelpassports.sanger.ac.uk/downloads, version: rnaseq_20191101) or
## whose essentiality is statistically associated with a molecular feature
## (thus very likely to be context specific essential genes) as negative controls

load('data/CMP_RNAseq.Rdata')

## genes not expressed
lowlyExp<-names(which(rowSums(CMP_RNAseq<0.01,na.rm=TRUE)>=ncol(CMP_RNAseq)))

## genes that are differentially essential in the presence of a molecular feature,
## i.e. associated with a biomarker, thus unlikely to be CFGs
wBm <- sort(unique(unlist(read.table('data/dependency_with_biomarkers.txt',stringsAsFactors = FALSE))))

negativeControls<-union(lowlyExp,wBm)

negativeControls_includingTraining<-negativeControls

## removing training sets
negativeControls<-setdiff(negativeControls,TrainingSets)

print(paste(length(negativeControls_includingTraining),'genes lowly expressed or associated with a marker'))
print(paste(length(negativeControls),
            'genes not included in any of the training sets will be used as independent negative controls'))

## Assembling CFGs outputted by a baseline DM predictor
baselineCFGs <- lapply(1:ncol(bdep),function(n){
  names(which(rowSums(bdep)>=n))
})

## Computing sizes (in terms of number of included genes)
## of the baseline CFGs
baselineSizes<-unlist(lapply(baselineCFGs,length))

## Computing baseline TPRs
baselineTPRs<-100*unlist(lapply(baselineCFGs,function(x){
  length(intersect(x,positiveControls))/length(positiveControls)
}))

## Computing baseline FPRs
baselineFPRs<-100*unlist(lapply(baselineCFGs,function(x){
  length(intersect(x,negativeControls))/length(negativeControls)
}))


## Plotting sizes of baseline CFGs, curves of baseline TPRs/FPRs and of FPRs at given level of TPRs.

par(mfrow=c(2,2))
par(mar=c(6,4,1,1))
plot(baselineSizes,xlab='minimal n. cell lines dependent on the CFGs',
     ylab='n. of predicted CFGs',pch=16,log='y')

plot(baselineTPRs,ylab='Recall of positive controls (TPRs)',
     xlab='minimal n. cell lines dependent on the CFGs',pch=16)

plot(baselineFPRs,ylab='Recall of negative controls (FPRs)',
     xlab='minimal n. cell lines dependent on the CFGs',pch=16)

plot(baselineTPRs,baselineFPRs,xlab='Recall of positive controls (TPRs)',
     ylab='Recall of negative controls (FPRs)',pch=16)


# Computing observed TPRs and FPRs
observedTPRs<-100*unlist(lapply(novelCFs_sets[3:length(novelCFs_sets)],
                                function(x){length(intersect(x,positiveControls))/length(positiveControls)}))

observedFPRs<-100*unlist(lapply(novelCFs_sets[3:length(novelCFs_sets)],
                                function(x){length(intersect(x,negativeControls))/length(negativeControls)}))

# Barplotting observed TPRs/FPRs
par(mar=c(15,4,4,1))
barplot(observedTPRs/max(baselineTPRs),col=col[names(observedTPRs)],
        las=2,ylab='TPR / (max baseline TPR)',border=NA,main = 'Observed Relative TPRs')

barplot(observedFPRs/max(baselineFPRs),col=col[names(observedFPRs)],
        las=2,ylab='FPR / (max baseline FPR)', border=NA, main = 'Observed Relative FPRs')


# Plotting observed FPRs at observed level of TPRs, with respect to baseline FPRs at fixed level of TPRs.
par(mfrow=c(1,2))
par(mar=c(4,4,2,1))
plot(spline(baselineTPRs,baselineFPRs),pch=16,
     xlab='% Recall of positive controls (TPRs)',
     ylab='% Recall of negative controls (FPRs)',
     type='l',lwd=5,
     xlim=c(9,32),ylim=c(0.05,0.35))
legend('topleft',legend='baseline predictor',lwd = 5)

points(observedTPRs,observedFPRs,col=col[names(observedTPRs)],pch=16,cex=2)


# Barplotting ratios between observed FPRs and baseline FPRs at observed TPRs
s0fun<-splinefun(baselineTPRs,baselineFPRs)

par(mar=c(4,10,2,1))

barplot(observedFPRs/s0fun(observedTPRs),col=col[names(observedFPRs)],
        las=2,border=NA,xlab=paste('FPRs /\n(baseline FPRs at observed TPRs)'),
        horiz = TRUE,xlim=c(0,2))
abline(v=1,lty=2)


print(sort(observedFPRs/s0fun(observedTPRs)))

print(sort(observedFPRs/s0fun(observedTPRs)))

print(paste('median FPRs vs expectation ratio for unsupervised methods:',
            median((observedFPRs/s0fun(observedTPRs))[5:9])))

## Comparing median numbers of dependant cell lines across sets of predicted CFGs.
screenedGenes<-rownames(bdep)

median_n_dep_cell_lines<-
  unlist(lapply(novelCFs_sets[3:length(novelCFs_sets)],function(x){median(rowSums(bdep[intersect(x,screenedGenes),]))}))

median_perc_dep_cell_lines<-100*median_n_dep_cell_lines/ncol(bdep)

par(mar=c(12,4,2,1))
barplot(median_perc_dep_cell_lines,
        col=col[names(median_n_dep_cell_lines)],
        las=2, border=NA,main = 'CFGs Median % dep cell lines', ylab = '%',ylim=c(0,100))

print(sort(median_perc_dep_cell_lines),decreasing=TRUE)
print(paste('grand median percentaged of dependent cell lines for unsupervised methods:',
            median(median_perc_dep_cell_lines[5:9])))

# Comparing median dependent cell lines for the predicted CFGs and threshold n of the baseline
# DM classifier required to attain the observed TPRs
par(mar=c(6,4,1,1))
plot(100*1:length(baselineTPRs)/length(baselineTPRs),
     baselineTPRs,xlab='threshold n of baseline classifier (as % of screened cell lines)',
     ylab='% Recall of positive controls (TPRs)',pch=16)

s0fun<-splinefun(baselineTPRs,1:length(baselineTPRs))
abline(v=100*s0fun(observedTPRs)/length(baselineTPRs),col=col[names(observedTPRs)],lwd=2)

baseline_n_dep_cl_at_observed_TPR<-s0fun(observedTPRs)

par(mar=c(12,4,2,1))
barplot(median_n_dep_cell_lines/baseline_n_dep_cl_at_observed_TPR,
        col=col[names(median_n_dep_cell_lines)],
        las=2, border=NA,main = 'CFGs Median n.dep cell lines / baseline n at obs TPRs', ylab = 'ratio')
abline(h=1,lty=2)

print(median_n_dep_cell_lines/baseline_n_dep_cl_at_observed_TPR)

print(paste('median ratio for unsupervised methods:',
            median((median_n_dep_cell_lines/baseline_n_dep_cl_at_observed_TPR)[5:9])))

##  Comparing median fitness effects across predicted CFGs and comparison with baseline median
## fitness effect at observed TPRs (excluding training sets)
median_dep <- unlist(lapply(novelCFs_sets[3:11],
                            function(x){median(apply(scaled_depFC[intersect(x,screenedGenes),],1,mean))}))

print(sort(median_dep))

par(mar=c(12,6,4,2))
barplot(median_dep,col=col[names(median_dep)],
        las=2,main = 'Median CFGs fitness effect', border=NA, ylab = 'median fitness effect',ylim=c(-1,0))
abline(h=-1,lty=1)
abline(h=-0.5,lty=2)


baseline_dep <- unlist(lapply(baselineCFGs[round(s0fun(observedTPRs))],
                              function(x){median(apply(scaled_depFC[intersect(setdiff(x,TrainingSets),screenedGenes),],1,mean))}))

par(mar=c(12,6,4,2))
barplot(median_dep/baseline_dep,col=col[names(median_dep)],
        las=2,main = 'Median CFGs fitness effect / baseline', border=NA, ylab = 'ratio',ylim=c(0,1))
abline(h=1,lty=2)

print(sort(median_dep/baseline_dep))

print(median((median_dep/baseline_dep)[5:9]))

## Comparing number of CFGs across predicted set and compared with baseline CFGs at observed TPRs
baseline_size<-unlist(lapply(baselineCFGs[round(s0fun(observedTPRs))],function(x){length(setdiff(x,TrainingSets))}))

observed_sizes<-novelGeneLengths[3:11]
par(mar=c(12,6,4,2))
barplot(observed_sizes/baseline_size,col=col[names(observed_sizes)],
        las=2,main = 'n. CFGs / baseline at observed TPRs',border=NA, ylab = 'ratio')
abline(h=1,lty=2)



### Comparing basal expression of predicted CFGs in Normal tissues
load('data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.Rdata')

bsexp<-lapply(novelCFs_sets[3:11],function(x){
    tmp<-basalExprsNormalTissues[intersect(x,rownames(basalExprsNormalTissues)),]
    unlist(lapply(1:nrow(tmp),function(x){median(tmp[x,])}))
  })

par(mar=c(12,4,4,2))
boxplot(lapply(bsexp,function(x){log(x+1)}),las=2,col=col[names(bsexp)],ylab='grand median across tissues',
main='basal expression in normal tissues')#

##################################################################################################################
## Benchmarking including training sets
## The following code re-executes the benchmark including also the Hart2014 and Hart2017 sets
## and adding to the Sharma2020 and CENtools predicted CFGs the prior known CFGs used in the training phase, respectively Hart2017 and curated Hart2014
## all plots including training

names(CFs_sets_plus_training)[4]<-'CEN-tools (Sharma 2020) + Hart2017'
names(CFs_sets_plus_training)[5]<-'CEN-tools + curated Hart2017'
col_includingTraining<-col
names(col_includingTraining)[4]<-'CEN-tools (Sharma 2020) + Hart2017'
names(col_includingTraining)[5]<-'CEN-tools + curated Hart2017'
col<-c(col,col_includingTraining)

## Computing baseline TPRs Including Training sets
baselineTPRs_includingTraining<-100*unlist(lapply(baselineCFGs,function(x){
  length(intersect(x,positiveControls_includingTraining))/length(positiveControls_includingTraining)
}))

## Computing baseline FPRs Including Training sets
baselineFPRs_includingTraining<-100*unlist(lapply(baselineCFGs,function(x){
  length(intersect(x,negativeControls_includingTraining))/length(negativeControls_includingTraining)
}))

## Plotting baseline DM features Including Training sets
par(mfrow=c(2,2))
par(mar=c(6,4,1,1))
plot(baselineTPRs_includingTraining,ylab='Recall of positive controls (TPRs)',
     xlab='minimal n. cell lines dependent on the CFGs',pch=16)

plot(baselineFPRs_includingTraining,ylab='Recall of negative controls (FPRs)',
     xlab='minimal n. cell lines dependent on the CFGs',pch=16)

plot(baselineTPRs_includingTraining,baselineFPRs_includingTraining,xlab='Recall of positive controls (TPRs)',
     ylab='Recall of negative controls (FPRs)',pch=16)


## computing observed TPRs/FPRs including Training sets
observedTPRs_includingTraining<-100*unlist(lapply(CFs_sets_plus_training,
                                                  function(x){length(intersect(x,positiveControls_includingTraining))/length(positiveControls_includingTraining)}))

observedFPRs_includingTraining<-100*unlist(lapply(CFs_sets_plus_training,
                                                  function(x){length(intersect(x,negativeControls_includingTraining))/length(negativeControls_includingTraining)}))

## Barplotting observed TPRs/FPRs including Training sets
par(mfrow=c(1,1))
par(mar=c(15,4,4,1))
barplot(observedTPRs_includingTraining/max(baselineTPRs_includingTraining),col=col[names(observedTPRs_includingTraining)],
        las=2,ylab='TPR / (max baseline TPR)',border=NA,main = 'Observed Relative TPRs')

print(sort(observedTPRs_includingTraining/max(baselineTPRs_includingTraining),decreasing=TRUE))

print(paste('median relative TPRs for unsupervised methods:',
            median((observedTPRs_includingTraining/max(baselineTPRs_includingTraining)))))

barplot(observedFPRs_includingTraining/max(baselineFPRs_includingTraining),
        col=col[names(observedFPRs_includingTraining)],
        las=2,ylab='FPR / (max baseline FPR)', border=NA, main = 'Observed Relative FPRs')

print(sort(observedFPRs_includingTraining/max(baselineFPRs_includingTraining)))

print(paste('median relative FPRs for unsupervised methods:',
            median((observedFPRs_includingTraining/max(baselineFPRs_includingTraining)))))


## Plotting observed FPRs at observed level of TPRs, with respect to baseline FPRs at fixed level of TPRs.

par(mar=c(4,4,2,1))

plot(spline(c(0,baselineTPRs_includingTraining),c(0,baselineFPRs_includingTraining)),pch=16,
     xlab='% Recall of positive controls (TPRs)',
     ylab='% Recall of negative controls (FPRs)',
     type='l',lwd=5,ylim=c(0.05,0.55),xlim=c(15,37))

points(observedTPRs_includingTraining,observedFPRs_includingTraining,
       col=col[names(observedTPRs_includingTraining)],pch=16,cex=2)


#Barplotting ratios between observed FPRs and baseline FPRs at observed TPRs
s0fun<-splinefun(c(0,baselineTPRs_includingTraining),c(0,baselineFPRs_includingTraining))

par(mar=c(15,4,4,4))
barplot(observedFPRs_includingTraining/s0fun(observedTPRs_includingTraining),
        col=col[names(observedFPRs_includingTraining)],
        las=2,border=NA,ylab='Observed FPRs / baseline FPRs at observed level of TPRs')
abline(h=1,lty=2)

print(sort(observedFPRs_includingTraining/s0fun(observedTPRs_includingTraining)))

print(paste('median FPRs vs expectation ratio for unsupervised methods:',
            median((observedFPRs_includingTraining/s0fun(observedTPRs_includingTraining))[5:9])))

## Comparing numbers of dependant cell lines across sets of predicted CFGs.
screenedGenes<-rownames(bdep)

median_n_dep_cell_lines_includingTraining<-
  unlist(lapply(CFs_sets_plus_training,function(x){median(rowSums(bdep[intersect(x,screenedGenes),]))}))

median_perc_dep_cell_lines_includingTraining<-100*median_n_dep_cell_lines_includingTraining/ncol(bdep)

par(mar=c(15,4,2,1))
barplot(median_perc_dep_cell_lines_includingTraining,
        col=col[names(median_n_dep_cell_lines_includingTraining)],
        las=2, border=NA,main = 'CFGs Median % dep cell lines', ylab = '%',ylim=c(0,100))

print(sort(median_perc_dep_cell_lines_includingTraining),decreasing=TRUE)
print(paste('grand median percentaged of dependent cell lines for unsupervised methods:',
            median(median_perc_dep_cell_lines_includingTraining[5:9])))

## Comparing number of dependent cell lines for the predicted CFGs and threshold  𝑛  of minimal number of dependent cell lines required by the baseline DM classifier required to attain the observed TPRs.par(mar=c(6,4,1,1))
plot(100*1:length(baselineTPRs_includingTraining)/length(baselineTPRs_includingTraining),
     baselineTPRs_includingTraining,xlab='threshold n of baseline classifier (as % of screened cell lines)',
     ylab='% Recall of positive controls (TPRs)',pch=16)

s0fun<-splinefun(c(0,baselineTPRs_includingTraining),0:length(baselineTPRs))
abline(v=100*s0fun(observedTPRs_includingTraining)/length(baselineTPRs_includingTraining),
       col=col[names(observedTPRs_includingTraining)],lwd=2)

baseline_n_dep_cl_at_observed_TPR_includingTraining<-s0fun(observedTPRs_includingTraining)

par(mar=c(15,4,2,1))
barplot(median_n_dep_cell_lines_includingTraining/baseline_n_dep_cl_at_observed_TPR_includingTraining,
        col=col[names(median_n_dep_cell_lines_includingTraining)],
        las=2, border=NA,main = 'CFGs Median n.dep cell lines / baseline n at obs TPRs', ylab = 'ratio')
abline(h=1,lty=2)

print(median_n_dep_cell_lines_includingTraining/baseline_n_dep_cl_at_observed_TPR_includingTraining)

print(paste('median ratio for unsupervised methods:',
            median((median_n_dep_cell_lines_includingTraining/baseline_n_dep_cl_at_observed_TPR_includingTraining)[5:9])))

##  Comparing median fitness effects across predicted CFGs and comparison with baseline median
## fitness effect at observed TPRs (excluding training sets)
median_dep_includingTraining <- unlist(lapply(CFs_sets_plus_training,
                                              function(x){median(apply(scaled_depFC[intersect(x,screenedGenes),],1,mean))}))

par(mar=c(15,6,4,2))
barplot(median_dep_includingTraining,col=col[names(median_dep_includingTraining)],
        las=2,main = 'Median CFGs fitness effect', border=NA, ylab = 'median fitness effect')
abline(h=1,lty=2)

baseline_dep_includingTraining <- unlist(lapply(baselineCFGs[round(s0fun(observedTPRs_includingTraining))],
                                                function(x){median(apply(scaled_depFC[intersect(x,screenedGenes),],1,mean))}))

par(mar=c(15,6,4,2))
barplot(median_dep_includingTraining/baseline_dep_includingTraining,col=col[names(median_dep_includingTraining)],
        las=2,main = 'Median CFGs fitness effect / baseline median at observed TPRs', border=NA, ylab = 'ratio',ylim=c(0,1))
abline(h=1,lty=2)

##Comparing number of CFGs across predicted set and compared with baseline CFGs at observed TPRs

baseline_size_includingTraining<-
  unlist(lapply(baselineCFGs[round(s0fun(observedTPRs_includingTraining))],length))

observed_sizes_includingTraining<-
  unlist(lapply(CFs_sets_plus_training,length))

par(mar=c(15,6,4,2))
barplot(observed_sizes_includingTraining/baseline_size_includingTraining,
        col=col[names(observed_sizes_includingTraining)],
        las=2,main = 'n. CFGs / baseline at observed TPRs',border=NA, ylab = 'ratio')


##Comparing basal expression of predicted CFGs in Normal tissues

load('data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.Rdata')

bsexp_includingTraining<-lapply(CFs_sets_plus_training,function(x){
  tmp<-basalExprsNormalTissues[intersect(x,rownames(basalExprsNormalTissues)),]
  unlist(lapply(1:nrow(tmp),function(x){median(tmp[x,])}))
})

par(mar=c(15,4,4,2))
boxplot(lapply(bsexp_includingTraining,
               function(x){log(x+1)}),las=2,col=col[names(bsexp_includingTraining)],
        ylab='grand median across tissues',
        main='basal expression in normal tissues')

#################################################################################




### Comparison of covered prior known CFGs not included in any of the trainin sets
### across supervised methods
recalledPK_CFGs<-lapply(novelCFs_sets,
                   function(x){intersect(x,positiveControls)})

allNewHits<-unique(unlist(recalledPK_CFGs))

newHitsMemb<-do.call(rbind,lapply(recalledPK_CFGs,function(x){is.element(allNewHits,x)}))+0
colnames(newHitsMemb)<-allNewHits

par(mfrow=c(3,1))
vennDiagram(t(newHitsMemb[c(3,6),]),main='prior known CFGs not included in any of the trainin sets')
vennDiagram(t(newHitsMemb[c(4,6),]))
vennDiagram(t(newHitsMemb[c(5,6),]))






