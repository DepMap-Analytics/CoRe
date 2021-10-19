## Working directory needs to be /CoRe/notebooks

# system('git clone https://github.com/DepMap-Analytics/CoRe')
# setwd('CoRe/notebooks')
# source('CENtools/install.R')
# system('pip3 install -r CENtools/requirements.txt')

options(warn=-1)

library(tidyverse)
library(pheatmap)
library(CoRe)
library(CRISPRcleanR)
library(magrittr)
library(nVennR)
library(limma)
library(stringr)

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


## downloading cell lines' annotations
clannotation<-CoRe.download_AnnotationModel(URL='https://cog.sanger.ac.uk/cmp/download/model_list_20210326.csv')

## removing samples with NAs or not annotated
depFC <- depFC[,-which(is.na(colSums(depFC)))]
cells <- colnames(depFC)

tot_id <- clannotation %>% filter(model_id %in% cells | BROAD_ID %in% cells) %>%
  select(model_name,model_id,BROAD_ID) %>%
  pivot_longer(cols = c("model_id","BROAD_ID"), names_to = "institute", values_to = "model_id") %>%
  select(model_name,model_id) %>% filter(model_id %in% cells) %>%
  arrange(factor(model_id, levels = cells))

depFC <- depFC[,which(cells %in% tot_id$model_id)]

## assigning cell line names as column headers
colnames(depFC) <- tot_id$model_name

## scaling
scaled_depFC <- CoRe.scale_to_essentials(depFC,curated_BAGEL_essential,curated_BAGEL_nonEssential)

## deriving binary essentiality scores
bdep <- apply(scaled_depFC, 2, function(x){
  x[which(x >= -0.5)] <- 0
  x[which(x < -0.5)] <- 1
  x
})

print('Done: quantitative and binary dependency matrices have been created')
print(paste('(accounting for',dim(bdep)[1],'genes and',dim(bdep)[2],'cell lines)'))


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
#                            verbose=FALSE,
#                            TruePositives = curated_BAGEL_essential)

print(paste('Done: ADaM has identified',length(ADaM),'CFGs'))

## The procomputed FiPer CEGs are also available and can be loaded uncommenting
## and executing this line instead of the the following 5:
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

## consensus core-fitness genes across the three less stringent variants
#Perc_Consensus <- Reduce(intersect,list(Perc_fixed,Perc_slope,Perc_AUC))

print('Done')

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



## The procomputed CEN-tools Logistic Regression CFGs are also available and can be loaded uncommenting
## and executing this line instead of the the following ones:
load('data/preComputed/CENtools.RData')

## Identification of core essential genes by CEN-tools on the integrated CERES scaled dataset
#write.csv(scaled_depFC, 'CENtools/data/curated_data/CERES_scaled_depFC.csv', quote = FALSE)

## The following line installs all the python packages required to execute
## the logisting based regression method.
## Uncomment it if necessary.

## system('pip3 install -r CENtools/requirements.txt')

#system('mkdir -p CENtools/data/objects')
#system('python3 CENtools/LR.py', wait = TRUE)

#source('CENtools/clustering.R')

#project_path = 'CENtools/prediction_output/INTEGRATED/'
#CENtools <- ClusterEssentiality(Chosen_project= 'INTEGRATED',
#                                binPath= paste0(project_path, 'INTEGRATED_Histogram_DF_20_BIN.txt'),
#                                resultPath = project_path)

#system('rm -r CENtools/prediction_output/')
#system('rm CENtools/data/curated_data/CERES_scaled_depFC.csv')

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

names(CFs_sets)<-c('Hart2014',
                   'Hart2017',
                   'Behan2019',
                   'Sharma2020',
                   'CEN-tools',
                   'ADaM',
                   'FiPer Average',
                   'FiPer Consensus',
                   'FiPer Slope',
                   'FiPer AUC',
                   'FiPer Fixed')

col=c("#03B2C8","#034DD9","#E6AB02","#1B9E77","#337100","#FC8D62",
      '#F4CAE4','#800080','#CAB2D6','#BEBADA','#777892')
names(col)<-names(CFs_sets)


TrainingSets<-unique(c(BAGEL_essential,
                       BAGEL_nonEssential,
                       curated_BAGEL_essential,
                       curated_BAGEL_nonEssential))

print(paste(length(TrainingSets),'genes used for training by at least one method'))

novelCFs_sets<-lapply(CFs_sets,function(x){
  base::setdiff(x, TrainingSets)
})

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

print(paste('FiPer consensus n.genes (total): ',length(CFs_sets$`FiPer Consensus`)))
print(paste('FiPer consensus n.genes (novel hits):',length(novelCFs_sets$`FiPer Consensus`)))

par(mar=c(4,12,2,2))
barplot(rbind(novelGeneLengths,GeneLengths-novelGeneLengths),
        las=2,border = FALSE,horiz = TRUE,xlim=c(0,2300),
        xlab='n. genes')

legend('bottomright',legend = c('Training sets','(Hart2017 + curated Hart2014)'),
       fill=c('lightgray',NA),border=NA,bty = 'n')


Recall_Hart2014<-unlist(lapply(CFs_sets,function(x){
  100*length(base::intersect(x,CFs_sets$Hart2014))/length(CFs_sets$Hart2014)
}))

Recall_Hart2017<-unlist(lapply(CFs_sets,function(x){
  100*length(base::intersect(x,CFs_sets$Hart2017))/length(CFs_sets$Hart2017)
}))

Recall_Behan2019<-unlist(lapply(CFs_sets,function(x){
  100*length(base::intersect(x,CFs_sets$Behan2019))/length(CFs_sets$Behan2019)
}))

Recall_Sharma2020<-unlist(lapply(CFs_sets,function(x){
  100*length(base::intersect(x,CFs_sets$Sharma2020))/length(CFs_sets$Sharma2020)
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


CFs_sets_plus_training<-CFs_sets
CFs_sets_plus_training$Sharma2020 <-
  union(CFs_sets_plus_training$Sharma2020,BAGEL_essential)

CFs_sets_plus_training$`CEN-tools`<-
  union(CFs_sets_plus_training$`CEN-tools`,curated_BAGEL_essential)

##### CF gene set similarity
allEss<-unique(unlist(CFs_sets_plus_training))
membMat<-do.call(rbind,lapply(CFs_sets_plus_training,function(x){is.element(allEss,x)}))
rownames(membMat)<-names(CFs_sets)
colnames(membMat)<-allEss
membMat<-membMat+0

rownames(membMat)[which(rownames(membMat)=="Sharma2020")]<-"Sharma2020 + Hart 2017"
rownames(membMat)[which(rownames(membMat)=="CEN-tools")]<-"CEN-tools + curated Hart 2014"

pheatmap(membMat,show_colnames = FALSE,clustering_distance_rows = 'binary',cluster_cols = FALSE,border_color = NA,
         main = 'Similarity across sets',col=c('white','blue'),legend_labels = c('out','in'),legend_breaks = c(0,1))

dmat<-dist(membMat,method = 'binary')

pheatmap(1-as.matrix(dmat), main = 'JS for Pan-cancer core fitness genes across methods',
         legend = FALSE,display_numbers = round(dmat,digits = 2))

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
load("data/histones.RData")

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

positiveControls_includingTraining<-unique(positiveControls)

## removing training sets
positiveControls<-setdiff(positiveControls,TrainingSets)

print(paste(length(positiveControls_includingTraining),'genes in the independent reference sets'))
print(paste(length(positiveControls),'genes not included in any of the training sets will be used as independent positive controls'))




load('data/CMP_RNAseq.Rdata')

## genes not expressed
lowlyExp<-names(which(rowSums(CMP_RNAseq<0.01,na.rm=TRUE)>=ncol(CMP_RNAseq)))

## genes that are differentially essential in the presence of a molecular feature,
## i.e. associated with a biomarker, thus unlikely to be CFGs
wBm <- sort(unique(unlist(read.table('data/dependency_with_biomarkers.txt',stringsAsFactors = FALSE))))

negativeControls<-union(lowlyExp,wBm)

negativeControls_includingTraining<-unique(negativeControls)

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

observedTPRs<-100*unlist(lapply(novelCFs_sets[3:length(novelCFs_sets)],
                                function(x){length(intersect(x,positiveControls))/length(positiveControls)}))

observedFPRs<-100*unlist(lapply(novelCFs_sets[3:length(novelCFs_sets)],
                                function(x){length(intersect(x,negativeControls))/length(negativeControls)}))

par(mar=c(15,4,4,1))
barplot(observedTPRs/max(baselineTPRs),col=col[names(observedTPRs)],
        las=2,ylab='TPR / (max baseline TPR)',border=NA,main = 'Observed Relative TPRs')

print(sort(observedTPRs/max(baselineTPRs),decreasing=TRUE))

print(paste('median relative TPRs for unsupervised methods:',
            median((observedTPRs/max(baselineTPRs))[5:9])))

barplot(observedFPRs/max(baselineFPRs),col=col[names(observedFPRs)],
        las=2,ylab='FPR / (max baseline FPR)', border=NA, main = 'Observed Relative FPRs')

print(sort(observedFPRs/max(baselineFPRs)))

print(paste('median relative FPRs for unsupervised methods:',
            median((observedFPRs/max(baselineFPRs))[5:9])))

par(mar=c(4,4,2,1))
plot(spline(baselineTPRs,baselineFPRs),pch=16,
     xlab='% Recall of positive controls (TPRs)',
     ylab='% Recall of negative controls (FPRs)',
     type='l',lwd=5,
     xlim=c(9,32),ylim=c(0.05,0.35))

points(observedTPRs,observedFPRs,col=col[names(observedTPRs)],pch=16,cex=2)


s0fun<-splinefun(baselineTPRs,baselineFPRs)

par(mar=c(12,4,4,4))
barplot(observedFPRs/s0fun(observedTPRs),col=col[names(observedFPRs)],
        las=2,border=NA,ylab='Observed FPRs / baseline FPRs at observed level of TPRs')
abline(h=1,lty=2)

print(sort(observedFPRs/s0fun(observedTPRs)))

print(paste('median FPRs vs expectation ratio for unsupervised methods:',
            median((observedFPRs/s0fun(observedTPRs))[5:9])))

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

par(mar=c(12,6,4,2))
barplot(median_dep,col=col[names(median_dep)],
        las=2,main = 'Median CFGs fitness effect', border=NA, ylab = 'median fitness effect')
abline(h=1,lty=2)

baseline_dep <- unlist(lapply(baselineCFGs[round(s0fun(observedTPRs))],
                              function(x){median(apply(scaled_depFC[intersect(setdiff(x,TrainingSets),screenedGenes),],1,mean))}))

par(mar=c(12,6,4,2))
barplot(median_dep/baseline_dep,col=col[names(median_dep)],
        las=2,main = 'Median CFGs fitness effect / baseline median at observed TPRs', border=NA, ylab = 'ratio',ylim=c(0,1))
abline(h=1,lty=2)

baseline_size<-unlist(lapply(baselineCFGs[round(s0fun(observedTPRs))],function(x){length(setdiff(x,TrainingSets))}))

observed_sizes<-novelGeneLengths[3:11]
par(mar=c(12,6,4,2))
barplot(observed_sizes/baseline_size,col=col[names(observed_sizes)],
        las=2,main = 'n. CFGs / baseline at observed TPRs',border=NA, ylab = 'ratio')
abline(h=1,lty=2)

load('data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.Rdata')

bsexp<-lapply(novelCFs_sets[3:11],function(x){
  tmp<-basalExprsNormalTissues[intersect(x,rownames(basalExprsNormalTissues)),]
  unlist(lapply(1:nrow(tmp),function(x){median(tmp[x,])}))
})

par(mar=c(12,4,4,2))
boxplot(lapply(bsexp,function(x){log(x+1)}),las=2,col=col[names(bsexp)],ylab='grand median across tissues',
        main='basal expression in normal tissues')


names(CFs_sets_plus_training)[4]<-'Sharma2020 + Hart2017'
names(CFs_sets_plus_training)[5]<-'CEN-tools + curated Hart2014'
col_includingTraining<-col
names(col_includingTraining)[4]<-'Sharma2020 + Hart2017'
names(col_includingTraining)[5]<-'CEN-tools + curated Hart2014'
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
par(mar=c(15,4,4,1))
barplot(observedTPRs_includingTraining/max(baselineTPRs_includingTraining),
        col=col[names(observedTPRs_includingTraining)],
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
     type='l',lwd=5,ylim=c(0.05,0.55),xlim=c(21,53))

points(observedTPRs_includingTraining,observedFPRs_includingTraining,
       col=col[names(observedTPRs_includingTraining)],pch=16,cex=2)

abline(v=min(baselineTPRs_includingTraining),lty=2)
abline(h=min(baselineFPRs_includingTraining),lty=2)
legend('topleft','basal DM n* = 1 cell line',lty=2)

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

nc_cfgs<-setdiff(intersect(baselineCFGs[[ncol(bdep)]],
                           positiveControls_includingTraining),
                 CFs_sets_plus_training$`Hart 2014`)

print(paste(length(nc_cfgs),
            'always depleted positive controls not covered by Hart2014:'))
par(mar=c(16,5,4,1))
barplot(100*unlist(lapply(CFs_sets_plus_training,
                          function(x){length(intersect(x,nc_cfgs))/length(nc_cfgs)})),
        col=col_includingTraining[names(CFs_sets_plus_training)],las=2,
        ylab=paste('% of always essential positive\ncontrols not covered by Hart2014'))
abline(h=100,lty=2)

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


## Downloading DEMETER RNAi Cancer Dependency dataset, version DEMETER data v6, 04/20
url <- 'https://ndownloader.figshare.com/files/11489669'
temp <- tempfile()
download.file(url, temp, mode="wb")

## Loading and formatting DEMETER RNAi Cancer Dependency dataset
DEMETER_dep<-as.matrix(read.csv(temp,stringsAsFactors = FALSE,row.names = 1))
gnames<-rownames(DEMETER_dep)
gnames<-unlist(lapply(str_split(gnames,' '),function(x){x[1]}))
rownames(DEMETER_dep)<-gnames

scaled_DEMETER_de<-
  CoRe.scale_to_essentials(DEMETER_dep,ess_genes = curated_BAGEL_essential,noness_genes = curated_BAGEL_nonEssential)


setToconsider<-CFs_sets_plus_training[3:11]

medianRNAi_dep<-lapply(setToconsider,function(x){subM<-apply(scaled_DEMETER_de[intersect(x,rownames(scaled_DEMETER_de)),],
                                                             MARGIN=1,'median',na.rm=TRUE)})
## barplotting median DEMETER dependency scores
par(mar=c(15,5,4,2))
boxplot(medianRNAi_dep,las=2,
        col=col_includingTraining[names(medianRNAi_dep)],
        ylab='Median DEMETER gene fitness score')

grandMedian<-unlist(lapply(medianRNAi_dep,'median'))

print('absolute grand median:')
print(sort(grandMedian))

print('median for FiPer variants:')
print(median(grandMedian[5:9]))


## Computing baseline DM DEMETER depenendency scores at the observed TPRs
setToconsider<-baselineCFGs[s0fun(observedTPRs_includingTraining[3:11])]

baselineDM_medianRNAi_dep<-lapply(setToconsider,function(x){subM<-
  apply(scaled_DEMETER_de[intersect(x,rownames(scaled_DEMETER_de)),],
        MARGIN = 1,
        'median',na.rm=TRUE)})

baselinegrandMedian<-unlist(lapply(baselineDM_medianRNAi_dep,'median'))

print('grand median / baseline:')
print(sort(grandMedian/baselinegrandMedian,decreasing=TRUE))

print('median for FiPer variants:')
print(median(grandMedian[5:9]/baselinegrandMedian[5:9]))
## barplotting median DEMETER dependency scores / baseline
par(mar=c(15,5,4,2))
barplot(grandMedian/baselinegrandMedian,col=col_includingTraining[names(medianRNAi_dep)],
        ylab='Median DEMETER gene fitness score\n/ baseline DM',las=2)
abline(h=1)

my.hypTest<-function(x,k,n,N){


  PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)

  return(PVALS)
}
retrievegeneInfo<-function(GENES, fc){

  commong<-intersect(rownames(fc),GENES)

  fc<-fc[commong,]

  description<-rep('',length(GENES))
  names(description)<-GENES

  hgnc_id<-rep('',length(GENES))
  names(hgnc_id)<-GENES

  entrez_id<-rep('',length(GENES))
  names(entrez_id)<-GENES

  ensemble_id<-rep('',length(GENES))
  names(ensemble_id)<-GENES

  location<-rep('',length(GENES))
  names(location)<-GENES

  family<-rep('',length(GENES))
  names(family)<-GENES

  pubmed_id<-rep('',length(GENES))
  names(pubmed_id)<-GENES

  description[commong]<-fc[,'name']
  family[commong]<-fc[,'gene_family']
  hgnc_id[commong]<-fc[,'hgnc_id']
  entrez_id[commong]<-fc[,'entrez_id']
  pubmed_id[commong]<-fc[,'pubmed_id']

  res<-cbind(description,family,hgnc_id,entrez_id,pubmed_id)

  return(res)
}
enrichedGeneFamilies<-function(geneset,BGgs,fc){
  geneInfos<-retrievegeneInfo(BGgs,fc)[,2]
  tokened<-str_split(geneInfos,'[|]')
  names(tokened)<-names(geneInfos)

  observed<-tokened[geneset]

  observed<-observed[observed!='']

  k<-length(observed)

  tokened<-tokened[tokened!='']

  N<-length(tokened)

  toTest<-summary(as.factor(unlist(observed)))
  toTest<-sort(toTest[which(toTest>1)],decreasing=TRUE)
  toTest<-names(toTest)
  toTest<-setdiff(toTest,'(Other)')


  RES<-do.call(rbind,lapply(toTest,function(x){
    n<-length(which((unlist(lapply(lapply(tokened,'intersect',x),'length')))>0))
    currSet<-names(which((unlist(lapply(lapply(observed,'intersect',x),'length')))>0))
    x<-length(currSet)
    pvals<-my.hypTest(x,k,n,N)
    GENES<-paste(sort(currSet),collapse=', ')
    res<-data.frame(pval=pvals,Recall=x/n,Covered=x,n_in_category=n,n_picked=k,n_background=N,GENES=GENES,stringsAsFactors = FALSE)
  }))

  rownames(RES)<-toTest

  fdr<-100*p.adjust(RES$pval,'fdr')

  RES<-cbind(rownames(RES),fdr,RES)
  colnames(RES)[1]<-'Gene Family'

  RES[,1]<-as.character(RES[,1])
  RES<-RES[order(RES$fdr),]

  return(RES)
}

## loading pre-computed CFG sets. Alternatively comment this line and execute the whole notebook down to this point
load('data/preComputed/CFs_sets_plus_training.RData')

## Loading gene set to be used as background population, i.e. all screened genes in the DepMap dataset
load('data/screenedGenes.RData')

## Deriving gene annotations
fc_annotation<-read.table('data/protein-coding_genes_annot.txt',sep = '\t',header=TRUE,stringsAsFactors = FALSE)
rownames(fc_annotation)<-fc_annotation$symbol

## declaring a vector of colors
col_includingTraining=c("#03B2C8","#034DD9","#E6AB02","#1B9E77","#337100","#FC8D62",
                        '#F4CAE4','#800080','#CAB2D6','#BEBADA','#777892')
names(col_includingTraining)<-names(CFs_sets_plus_training)

## Deriving gene info tables for all the CFG/CEG sets
CFGs_infos<-lapply(CFs_sets_plus_training,retrievegeneInfo,fc_annotation)

## Performing enrichment analysis of gene families
GFs<-lapply(CFs_sets_plus_training,function(x){enrichedGeneFamilies(x,screenedGenes,fc_annotation)})
names(GFs)<-names(CFs_sets_plus_training)

## Extracting significant hits only (FDR < 5%)
NN<-names(CFs_sets_plus_training)
significantOnly<-lapply(NN,function(nn){

  x<-GFs[[nn]]
  x<-x[which(x$fdr<5),]

  return(x$`Gene Family`)
})

## Extracting families enriched in all supervised methods and state-of-the-art CFGs sets
always_enriched_CFGs<-Reduce(intersect,significantOnly[1:6])


## Computing recall of the always enriched in CFG families across methods and sets
RecallOfEnrichedFamilies<-do.call(rbind,lapply(always_enriched_CFGs,function(x){
  unlist(lapply(GFs,function(g){g[x,'Covered']}))
}))

rownames(RecallOfEnrichedFamilies)<-always_enriched_CFGs

## loading 13 distinct colors
load('data/col_13_distinct.RData')
par(mar=c(17,4,4,2), pty = "s")
par(xpd=TRUE)

RecallOfEnrichedFamilies <- RecallOfEnrichedFamilies[,order(colSums(RecallOfEnrichedFamilies), decreasing = TRUE)]

barplot(RecallOfEnrichedFamilies,col=col_13_distinct,xlim=c(0,400),xlab='n. genes',main='Coverage of gene families
        always enriched in CFGs (FDR < 5%)', border = NA, las = 1, horiz = TRUE)
plot(0,0,col=NA,frame.plot=FALSE,xlab='',ylab='',xaxt='n',yaxt='n')
legend('center',rownames(RecallOfEnrichedFamilies),fill=col_13_distinct,title='gene families always enriched in CFGs (FDR < 5%)')

print(sort(colSums(RecallOfEnrichedFamilies),decreasing = TRUE))

## Extracting families enriched in all unsupervised methods
always_enriched_CEGs<-Reduce(intersect,significantOnly[7:11])

allGinfos<-retrievegeneInfo(screenedGenes,fc_annotation)

## computing recalls of always enriched families
RecallOfEnrichedFamilies_CEGs<-do.call(rbind,lapply(always_enriched_CEGs,function(x){
  genes_in_the_family<-rownames(allGinfos)[grep(x,allGinfos[,'family'])]
  unlist(lapply(CFs_sets_plus_training,function(y){length(intersect(y,genes_in_the_family))}))
}))

rownames(RecallOfEnrichedFamilies_CEGs)<-always_enriched_CEGs
colnames(RecallOfEnrichedFamilies_CEGs)<-names(CFs_sets_plus_training)

## Plotting results
## loading 57 distinct colors
load(file='data/col_57_distinct.RData')
par(mar=c(17,4,6,2), pty = "s")
par(xpd=TRUE)

RecallOfEnrichedFamilies_CEGs <- RecallOfEnrichedFamilies_CEGs[,order(colSums(RecallOfEnrichedFamilies_CEGs), decreasing = TRUE)]

barplot(RecallOfEnrichedFamilies_CEGs,col=col_57_distinct,
        xlim=c(0,800),xlab='n. genes',main='Coverage of gene families always enriched in CEGs (FDR < 5%)',
        border=NA, horiz = TRUE, las = 1)
par(mar=c(4,2,3,2))
plot(0,0,col=NA,frame.plot=FALSE,xlab='',ylab='',xaxt='n',yaxt='n')
legend('center',rownames(RecallOfEnrichedFamilies_CEGs),cex=0.55,
       fill=col_57_distinct,title='gene families always enriched in CEGs (FDR < 5%)',border = NA)

## Loading time-point essential genes identified by Tzelepis et al., 2016

load('data/early_ess.RData')
load('data/mid_ess.RData')
load('data/late_ess.RData')

## Determining early/mid/late essential gene families
early_ess_enrich_families<-enrichedGeneFamilies(early_ess,screenedGenes,fc_annotation)
early_ess_enrich_families<-early_ess_enrich_families$`Gene Family`[which(early_ess_enrich_families$fdr<5)]

mid_ess_enrich_families<-enrichedGeneFamilies(mid_ess,screenedGenes,fc_annotation)
mid_ess_enrich_families<-mid_ess_enrich_families$`Gene Family`[which(mid_ess_enrich_families$fdr<5)]

late_ess_enrich_families<-enrichedGeneFamilies(late_ess,screenedGenes,fc_annotation)
late_ess_enrich_families<-late_ess_enrich_families$`Gene Family`[which(late_ess_enrich_families$fdr<5)]

common_CFGs_CEGs_time_comp<-c(length(intersect(intersect(always_enriched_CEGs,always_enriched_CFGs),early_ess_enrich_families)),
                              length(intersect(intersect(always_enriched_CEGs,always_enriched_CFGs),union(mid_ess_enrich_families,late_ess_enrich_families))))
CEGs_exclusive_time_comp<-c(length(intersect(setdiff(always_enriched_CEGs,always_enriched_CFGs),early_ess_enrich_families)),
                            length(intersect(setdiff(always_enriched_CEGs,always_enriched_CFGs),union(mid_ess_enrich_families,late_ess_enrich_families))))

barplot(100*cbind(common_CFGs_CEGs_time_comp/sum(common_CFGs_CEGs_time_comp),
                  CEGs_exclusive_time_comp/sum(CEGs_exclusive_time_comp)),
        names.arg = c('CFG and CEG','CEG only'),main='Enriched Gene Families',ylab='%',
        col=c('darkgray','lightblue'),border=NA)
plot(0,0,col=NA,frame.plot=FALSE,xlab='',ylab='',xaxt='n',yaxt='n')
legend('center',c('early essential gene families','mid/late essential gene families'),
       fill=c('darkgray','lightblue'),border = NA)

# ## Executing BAGEL on CRISPRcleanR corrected FCs
# system('bash data/run_bagel.sh', wait = TRUE)

# ## Assemble BF matrices
# CFGset <- list.dirs("data/BAGEL_output/")[-1]
# CFGlabs <- gsub("data/BAGEL_output//","",CFGset)

# common_g <- scan('data/common_genes.txt',what = "character")
# cellLines <- list.files("data/BAGEL_output/ADaM/")

# BFprofile <- list()

# for (i in 1:length(CFGset)){
#   for (j in 1:length(cellLines)){
#     if (j == 1){
#       BFprofile[[i]] <- read_tsv(paste0(CFGset[i],"/",cellLines[j]))
#       BFprofile[[i]] <- BFprofile[[i]] %>% select("GENE","BF") %>% group_by(GENE) %>% summarise(BF = mean(BF))
#     } else {
#       tmp <- read_tsv(paste0(CFGset[i],"/",cellLines[j]))
#       tmp <- tmp %>% select("GENE","BF") %>% group_by(GENE) %>% summarise(BF = mean(BF))

#       BFprofile[[i]] <- full_join(BFprofile[[i]],tmp,by = "GENE")
#     }
#   }

#   BFprofile[[i]] <- column_to_rownames(BFprofile[[i]], var = "GENE")
#   colnames(BFprofile[[i]]) <- gsub("_corrected_logFCs.tsv","",cellLines)

#   BFprofile[[i]] <- BFprofile[[i]][common_g,]
#   BFprofile[[i]] <- round(BFprofile[[i]], digits = 5)
# }

# names(BFprofile) <- CFGlabs

# ## cell-wise scaling
# data("curated_BAGEL_essential")
# data("curated_BAGEL_nonEssential")

# load('data/preComputed/CFs_sets_plus_training.RData')

# Sets <- CFs_sets_plus_training[1:6]
# names(Sets) <- c("Hart2014","Hart2017","Behan2019","Sharma2020","CENtools","ADaM")
# Sets[[7]] <- curated_BAGEL_essential
# names(Sets)[7] <- "curated_Hart2014"

# for (i in 1:length(BFprofile)){
#   CFGlabs <- names(BFprofile)[i]
#   common_g <- rownames(BFprofile[[i]])

#   for (j in 1:ncol(BFprofile[[i]])){
#     FCprofile <- -BFprofile[[i]][,j]
#     names(FCprofile) <- common_g

#     sigthr <- -ccr.PrRc_Curve(FCprofile,Sets[[CFGlabs]],
#                               curated_BAGEL_nonEssential,FDRth = 0.05, display = FALSE)$sigthreshold

#     BFprofile[[i]][,j] <- BFprofile[[i]][,j] - sigthr
#   }
#   BFprofile[[i]] <- round(BFprofile[[i]], digits = 5)

#   assign(CFGlabs, BFprofile[[i]])
#   do.call(save, list(CFGlabs, file = paste0("data/BF_matrix/Sanger_scaled_",CFGlabs,'.RData')))
# }

# system('rm -r data/BAGEL_output')

load('data/PANCAN_simple_MOBEM.rdata')
load('data/CMP_RNAseq.Rdata')

fullPathCFG <- list.files('data/BF_matrix/',full.names = TRUE)

CFGlabs <- gsub("\\.RData","",fullPathCFG)
CFGlabs <- strsplit(CFGlabs,"/")
CFGlabs <- sapply(CFGlabs,function(x){x[length(x)]})

BFSs <- lapply(1:length(fullPathCFG),function(x){
  load(fullPathCFG[x],.GlobalEnv)
  get(CFGlabs[x])
})

CFGlabs <- gsub("Sanger_scaled_","",CFGlabs)
names(BFSs) <- CFGlabs

BFs <- BFSs$ADaM

clannotation<-CoRe.download_AnnotationModel()

colnames(MoBEM)<-clannotation$model_name[match(colnames(MoBEM),clannotation$COSMIC_ID)]
colnames(CMP_RNAseq)<-clannotation$model_name[match(colnames(CMP_RNAseq),clannotation$model_id)]
cls<-intersect(intersect(colnames(BFs),colnames(MoBEM)),colnames(CMP_RNAseq))

BFs<-BFs[,cls]
MoBEM<-MoBEM[,cls]
CMP_RNAseq<-CMP_RNAseq[,cls]

## typically copy number amplified oncogenes
MoBEM["ERBB2_mut",]<-sign(MoBEM["ERBB2_mut",]+MoBEM["gain:cnaPANCAN301 (CDK12,ERBB2,MED24)",])
MoBEM["MYC_mut",]<-sign(MoBEM["MYC_mut",]+MoBEM["gain:cnaPANCAN91 (MYC)",])
MoBEM["MYCN_mut",]<-sign(MoBEM["MYCN_mut",]+MoBEM["gain:cnaPANCAN344 (MYCN)",])
MoBEM["EGFR_mut",]<-sign(MoBEM["EGFR_mut",]+MoBEM["gain:cnaPANCAN124 (EGFR)",])
MoBEM["KRAS_mut",]<-sign(MoBEM["KRAS_mut",]+MoBEM["gain:cnaPANCAN164 (KRAS)",])

MoBEM<-MoBEM[grep('_mut',rownames(MoBEM)),]
rownames(MoBEM)<-unlist(lapply(str_split(rownames(MoBEM),'_mut'),function(x){x[1]}))

compendium<-read.table('data/Compendium_Cancer_Genes.tsv',sep='\t',header=F,stringsAsFactors = FALSE)
colnames(compendium) <- c("SYMBOL","ROLE")
ACT<-compendium$SYMBOL[which(compendium$ROLE=='Act')]
LoF<-compendium$SYMBOL[which(compendium$ROLE=='LoF')]
ambigous<-compendium$SYMBOL[which(compendium$ROLE=='ambigous')]

OG<-setdiff(unique(setdiff(ACT,LoF)),ambigous)

MoBEM<-MoBEM[intersect(intersect(rownames(MoBEM),OG),rownames(BFs)),]
MoBEM<-MoBEM[which(rowSums(MoBEM)>0),]

CMP_RNAseq<-CMP_RNAseq[intersect(rownames(CMP_RNAseq),rownames(BFs)),]
lowExpMat<-CMP_RNAseq<0.1
CMP_RNAseq<-CMP_RNAseq[which(rowSums(lowExpMat)>0),]
CMP_RNAseq<-CMP_RNAseq[intersect(rownames(CMP_RNAseq),rownames(MoBEM)),]
CMP_RNAseq<-(CMP_RNAseq<0.1)+0

controls<-MoBEM

for (i in 1:nrow(CMP_RNAseq)){
  controls[rownames(CMP_RNAseq)[i],]<-controls[rownames(CMP_RNAseq)[i],]-CMP_RNAseq[i,]
}

## positive instance --> oncogene mutated and expressed in the cell line
## negative instance --> oncogene wild-type and not expressed in the cell line
classRes<-
  do.call(rbind,lapply(names(BFSs),function(x){
    localcontrols<-c(controls)
    localBFs<-c(as.matrix(BFSs[[x]][rownames(controls),colnames(controls)]))

    id <- which(localcontrols!=0)
    localcontrols <- localcontrols[id]
    localBFs <- localBFs[id]

    localcontrols<-localcontrols[order(localBFs,decreasing = TRUE)]
    localBFs<-sort(localBFs,decreasing=TRUE)
    dd<-t.test(localBFs~localcontrols)

    predictions<- -localBFs
    names(predictions)<-as.character(1:length(predictions))
    posc<-names(predictions)[which(localcontrols==1)]
    negc<-names(predictions)[which(localcontrols==-1)]
    AUROC<-ccr.ROC_Curve(FCsprofile = predictions,posc,negc,expName = x,FDRth = 0.05, display = F)
    AUPRC<-ccr.PrRc_Curve(FCsprofile = predictions,posc,negc,expName = x,FDRth = 0.05, display = F)

    return(list(tt.pval=dd$p.value,AUROC=AUROC$AUC,AUPRC=AUPRC$AUC,
                Recall_at_fixedFDR=AUROC$Recall))
  }))

rownames(classRes)<-names(BFSs)

col=c("#F4CAE4","#03B2C8","#034DD9","#E6AB02","#1B9E77","#337100","#FC8D62")
ordering <- c('curated_Hart2014', 'Hart2014', 'Hart2017', 'Behan2019', 'Sharma2020', 'CENtools', 'ADaM')

par(mar=c(12,4,4,4))

auprc <- unlist(classRes[,3])
auprc <- auprc[ordering]
barplot(auprc,las=2,ylim=c(0.72,0.77),col = col,main='AUPRC')

recall_5_fdr <- unlist(classRes[,4])
names(recall_5_fdr) <- gsub('\\.sensitivity','',names(recall_5_fdr))
recall_5_fdr <- recall_5_fdr[ordering]

par(mar=c(12,4,4,4))
barplot(recall_5_fdr,las=2,ylim=c(0.53,0.58),col = col,main='Recall at 5% FDR')
