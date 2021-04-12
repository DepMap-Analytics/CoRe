
my.hypTest<-function(x,k,n,N){


  PVALS<-phyper(x-1,n,N-n,k,lower.tail=FALSE)

  return(PVALS)
}

retrievegeneInfo<-function(GENES){

  fc<-read.table('data/protein-coding_genes_annot.txt',sep = '\t',header=TRUE,stringsAsFactors = FALSE)
  rownames(fc)<-fc$symbol

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
  ensemble_id[commong]<-fc[,'ensembl_gene_id']
  location[commong]<-fc[,'location_sortable']
  pubmed_id[commong]<-fc[,'pubmed_id']

  res<-cbind(description,family,hgnc_id,entrez_id,ensemble_id,location,pubmed_id)

  return(res)
}
enrichedGeneFamilies<-function(geneset,BGgs){
  geneInfos<-retrievegeneInfo(BGgs)[,2]
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

load('data/preComputed/CFs_sets_plus_training.RData')
load('data/preComputed/CFs_sets.RData')

load('data/screenedGenes.RData')

col_includingTraining=c("#03B2C8","#034DD9","#E6AB02","#1B9E77","#337100","#FC8D62",
      '#F4CAE4','#800080','#CAB2D6','#BEBADA','#777892')
names(col_includingTraining)<-names(CFs_sets_plus_training)

CFGs_infos<-lapply(CFs_sets_plus_training,retrievegeneInfo)
GFs<-lapply(CFs_sets_plus_training,function(x){enrichedGeneFamilies(x,screenedGenes)})
names(GFs)<-names(CFs_sets_plus_training)

CFGs_infos_no_training<-lapply(CFs_sets,retrievegeneInfo)
names(CFGs_infos_no_training)<-names(CFs_sets)

NN<-names(CFGs_infos_no_training)

lapply(NN,function(x){
  tmp<-CFGs_infos_no_training[[x]]
  tmp<-cbind(rownames(tmp),tmp)
  colnames(tmp)[1]<-'Gene'

  write.table(tmp,row.names = FALSE,
              quote=FALSE,
              sep='\t',file=paste('../../../Other Paper Results/SuppTable1/',x,'_set.tsv',sep=''))
})

significantOnly<-lapply(NN,function(nn){
  print(nn)
  pdf(paste('../../../Other Paper Results/',nn,'_GF_enr.pdf',sep=''),11,15)
  par(mar=c(4,20,4,6))
  par(xpd=TRUE)

  x<-GFs[[nn]]
  x<-x[which(x$fdr<5),]
  tt<-barplot(rev(-log10(x$pval)),horiz = TRUE,xlim=c(1,10),log='x',
              names.arg = rev(x$`Gene Family`),las=2,border = FALSE,
              xlab='-log10 pval',main = paste(nn,'(FDR < 5%)'))

  tex<-rev(x$GENES)
  nchar<-100
  toCut<-which(str_length(tex)>nchar)
  tex[toCut]<-str_sub(tex[toCut],1,nchar)
  tex[toCut]<-paste(tex[toCut],', ...',sep='')
  text(rep(1,length(tt)),tt-0.06,tex,pos = 4,cex=0.7)
  dev.off()

  return(x$`Gene Family`)
})

always_enriched_CFGs<-Reduce(intersect,significantOnly[1:6])

RecallOfEnrichedFamilies<-do.call(rbind,lapply(always_enriched_CFGs,function(x){
  unlist(lapply(GFs,function(g){g[x,'Covered']}))
}))

rownames(RecallOfEnrichedFamilies)<-always_enriched_CFGs

load('data/col_13_distinct.RData')
par(mar=c(15,4,4,2))
par(xpd=TRUE)
barplot(RecallOfEnrichedFamilies,las=2,col=col_13_distinct,ylim=c(0,400),ylab='n. genes',main='Coverage of gene families
        always enriched in CFGs (FDR < 5%)')
plot(0,0,col=NA,frame.plot=FALSE,xlab='',ylab='',xaxt='n',yaxt='n')
legend('center',rownames(RecallOfEnrichedFamilies),fill=col_13_distinct,title='gene families always enriched in CFGs (FDR < 5%)')

always_enriched_CEGs<-Reduce(intersect,significantOnly[7:11])

allGinfos<-retrievegeneInfo(screenedGenes)

RecallOfEnrichedFamilies_CEGs<-do.call(rbind,lapply(always_enriched_CEGs,function(x){
  genes_in_the_family<-rownames(allGinfos)[grep(x,allGinfos[,'family'])]
  unlist(lapply(CFs_sets_plus_training,function(y){length(intersect(y,genes_in_the_family))}))
}))

rownames(RecallOfEnrichedFamilies_CEGs)<-always_enriched_CEGs
colnames(RecallOfEnrichedFamilies_CEGs)<-names(CFs_sets_plus_training)

load(file='data/col_57_distinct.RData')
par(mar=c(15,4,6,2))
par(xpd=TRUE)
barplot(RecallOfEnrichedFamilies_CEGs,las=2,col=col_57_distinct,
  ylim=c(0,800),ylab='n. genes',main='Coverage of gene families always enriched in CEGs (FDR < 5%)',
  border=NA)
par(mar=c(4,2,3,2))
plot(0,0,col=NA,frame.plot=FALSE,xlab='',ylab='',xaxt='n',yaxt='n')
legend('center',rownames(RecallOfEnrichedFamilies_CEGs),cex=0.55,
       fill=col_57_distinct,title='gene families always enriched in CEGs (FDR < 5%)',border = NA)

####################################################################################################








