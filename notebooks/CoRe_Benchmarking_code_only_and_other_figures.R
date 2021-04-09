options(warn=-1)

library(tidyverse)
library(pheatmap)
library(CoRe)
library(magrittr)
library(nVennR)

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

write.csv(scaled_depFC, 'CENtools/data/curated_data/CERES_scaled_depFC.csv', quote = FALSE)

system('pip3 install -r CENtools/requirements.txt')

system('mkdir -p CENtools/data/objects')
system('python3 CENtools/LR.py', wait = TRUE)

source('CENtools/clustering.R')

project_path = 'CENtools/prediction_output/INTEGRATED/'
CENtools <- ClusterEssentiality(Chosen_project= 'INTEGRATED',
                                binPath= paste0(project_path, 'INTEGRATED_Histogram_DF_20_BIN.txt'),
                                resultPath = project_path)

system('rm -r CENtools/prediction_output/')
system('rm CENtools/data/curated_data/CERES_scaled_depFC.csv')

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

pheatmap(membMat,show_colnames = FALSE,clustering_distance_rows = 'binary',cluster_cols = FALSE,border_color = NA,
         main = 'Similarity across sets',col=c('white','blue'))

dmat<-dist(membMat,method = 'binary')

pheatmap(1-as.matrix(dmat), main = 'JS for Pan-cancer core fitness genes across methods',
         legend = FALSE,display_numbers = round(dmat,digits = 2))




