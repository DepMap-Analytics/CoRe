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

  pvals<-vector()
  GENES<-vector()

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

  RES<-RES[order(RES$fdr),]

  return(RES)
}
enrichedPathways<-function(geneset,BGgs,ming=2){

  k<-length(geneset)
  N<-length(BGgs)

  pvals<-vector()
  GENES<-vector()

  nt<-length(PATHCOM_HUMAN$PATHWAY)
  flag<-1


  ii<-unlist(lapply(lapply(PATHCOM_HUMAN$HGNC_SYMBOL,is.element,geneset),'sum'))
  names(ii)<-NULL

  toTest<-which(ii>=ming)

  nt<-length(toTest)
  testedP<-NULL
  for (i in 1:nt){
    print(i)
    n<-length(intersect(BGgs,PATHCOM_HUMAN$HGNC_SYMBOL[[toTest[i]]]))

    if (n >=ming){

      currSet<-intersect(geneset,PATHCOM_HUMAN$HGNC_SYMBOL[[toTest[i]]])
      x<-length(currSet)

      pvals[flag]<-my.hypTest(x,k,n,N)

      GENES[flag]<-paste(sort(currSet),collapse=', ')
      flag<-flag+1
      testedP<-c(testedP,as.character(PATHCOM_HUMAN$PATHWAY[[toTest[i]]]))
    }
  }


  id<-which(p.adjust(pvals,'fdr')<0.05)

  names(pvals)<-testedP
  pvals<-pvals[id]
  GENES<-GENES[id]

  GENES<-GENES[order(pvals)]
  pvals<-sort(pvals)

  RES<-cbind(names(pvals),pvals,GENES)
  RES<-data.frame(RES,stringsAsFactors = FALSE)
  RES$pvals<-as.numeric(as.character(RES$pvals))
  colnames(RES)[1]<-'pathway'
  rownames(RES)<-NULL
  return(RES)
}

data(curated_BAGEL_essential)


load('data/preComputed/CFs_sets_plus_training.RData')
load('data/screenedGenes.RData')

CFGs_infos<-lapply(CFs_sets_plus_training,retrievegeneInfo)
GFs<-lapply(CFs_sets_plus_training,function(x){enrichedGeneFamilies(x,screenedGenes)})
names(GFs)<-names(CFs_sets_plus_training)

x<-GFs$`Hart 2017`
x<-x[which(x$fdr<0.05),]
tt<-barplot(rev(-log10(x$pval)),horiz = TRUE,xlim=c(1,10),log='x',
            names.arg = rev(x$`Gene Family`),las=2,border = FALSE,
            xlab='-log10 pval')
text(rep(1,length(tt)),tt-0.06,rev(x$GENES),pos = 4,cex=0.7)




