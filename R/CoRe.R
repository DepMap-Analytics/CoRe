## Non Documented
CoRe.scale_to_essentials <- function(ge_fit,ess_genes,noness_genes){
  essential_indices <- which(row.names(ge_fit) %in% ess_genes)
  nonessential_indices <- which(row.names(ge_fit) %in% noness_genes)
  scaled_ge_fit <- ge_fit %>%
    apply(2, function(x){
      (x - median(x[nonessential_indices], na.rm=T)) %>%
        divide_by(median(x[nonessential_indices], na.rm=T) - median(x[essential_indices], na.rm=T))
    })
  return(scaled_ge_fit)
}


## Documented

## Apply CERES scaling on cell line fold change scores using two reference sets of essential and non-essential
## genes
scale_to_essentials <- function(ge_fit,ess_genes,noness_genes){
  essential_indices <- which(row.names(ge_fit) %in% ess_genes)
  nonessential_indices <- which(row.names(ge_fit) %in% noness_genes)

  scaled_ge_fit <- ge_fit %>%

    apply(2, function(x){
      (x - median(x[nonessential_indices], na.rm=T)) %>%
        divide_by(median(x[nonessential_indices], na.rm=T) - median(x[essential_indices], na.rm=T))
    })

  return(scaled_ge_fit)
}

## This function implements an heuristic algrorithm that takes in input a sparse binary matrix and sorts its rows
## and column in a way that the patterns of non null entries have a minimal overlap across rows.
HeuristicMutExSorting<-function(mutPatterns){

  mutPatterns<-sign(mutPatterns)

  ngenes<-nrow(mutPatterns)
  nsamples<-ncol(mutPatterns)

  if (is.null(rownames(mutPatterns))){
    rownames(mutPatterns) <- 1:ngenes
  }

  if (is.null(colnames(mutPatterns))){
    colnames(mutPatterns) <- 1:nsamples
  }

  if (nrow(mutPatterns)>1 & ncol(mutPatterns)>1){

    RowNull<-names(which(rowSums(mutPatterns)==0))
    RowNonNull<-which(rowSums(mutPatterns)>0)

    ColNull<-names(which(colSums(mutPatterns)==0))
    ColNonNull<-which(colSums(mutPatterns)>0)

    mutPatterns<-matrix(c(mutPatterns[RowNonNull,ColNonNull]),
                        length(RowNonNull),length(ColNonNull),
                        dimnames=list(rownames(mutPatterns)[RowNonNull],colnames(mutPatterns)[ColNonNull]))

    if (nrow(mutPatterns)>1 & ncol(mutPatterns)>1){

      coveredGenes<-NA
      uncoveredGenes<-rownames(mutPatterns)

      coveredSamples<-NA
      uncoveredSamples<-colnames(mutPatterns)
      BS<-NA

      while(length(uncoveredGenes)>0 & length(uncoveredSamples)>0){

        patterns<-matrix(c(mutPatterns[uncoveredGenes,uncoveredSamples]),
                         nrow = length(uncoveredGenes),
                         ncol = length(uncoveredSamples),
                         dimnames = list(uncoveredGenes,uncoveredSamples))

        if(length(uncoveredGenes)>1){
          bestInClass<-findBestInClass(patterns)
        }else{
          bestInClass<-uncoveredGenes
        }

        if(is.na(BS[1])){
          BS<-bestInClass
        }else{
          BS<-c(BS,bestInClass)
        }

        if(is.na(coveredGenes[1])){
          coveredGenes<-bestInClass
        }else{
          coveredGenes<-c(coveredGenes,bestInClass)
        }

        uncoveredGenes<-setdiff(uncoveredGenes,coveredGenes)
        toCheck<-matrix(c(patterns[bestInClass,uncoveredSamples]),nrow = 1,ncol=ncol(patterns),dimnames = list(bestInClass,uncoveredSamples))

        if (length(coveredGenes)==1){
          coveredSamples<-names(which(colSums(toCheck)>0))
        }else{
          coveredSamples<-c(coveredSamples,names(which(colSums(toCheck)>0)))
        }

        uncoveredSamples<-setdiff(uncoveredSamples,coveredSamples)

      }

      BS<-c(BS,uncoveredGenes)

      CID<-rearrangeMatrix(mutPatterns,BS)

      FINALMAT<-mutPatterns[BS,CID]

      nullCol<-matrix(0,nrow(FINALMAT),length(ColNull),
                      dimnames = list(rownames(FINALMAT),ColNull))

      FINALMAT<-cbind(FINALMAT,nullCol)

      nullRow<-matrix(0,length(RowNull),ncol(FINALMAT),
                      dimnames = list(RowNull,colnames(FINALMAT)))

      FINALMAT<-rbind(FINALMAT,nullRow)

      return(FINALMAT)

    } else {
      stop('Matrix must have at least 2 non-null rows and 2 non-null columns')
    }
  } else {
    stop('Matrix must have at least 2 rows and 2 columns')
  }
}

## This function finds the gene (i.e. row) with the highest exclusive coverage. The exclusive coverage for a gene g
## is defined as the number of uncovered samples in which this gene is mutated minus the number of samples in which
## at least another uncovered gene is mutated.
findBestInClass<-function(patterns){

  if(nrow(patterns)==1){
    return(rownames(patterns))
  }

  if(ncol(patterns)==1){
    return(rownames(patterns)[1])
  }

  exclCov<-colSums(t(2*patterns)-colSums(patterns))
  names(exclCov)<-rownames(patterns)

  return(names(sort(exclCov,decreasing=TRUE))[1])
}

## Rearrange Binary Matrix columns in order to minimise row-wise entry overlap based on exclusive coverage.
rearrangeMatrix<-function(patterns,
                          GENES){

  remainingSamples<-colnames(patterns)

  toAdd<-NULL

  DD<-t(t(2*patterns)-colSums(patterns))
  colnames(DD) <- remainingSamples

  for (g in GENES){
    cols <- remainingSamples[order(DD[g,remainingSamples],decreasing = TRUE)]
    toAdd<-c(toAdd,names(which(patterns[g,cols]>0)))
    remainingSamples<-setdiff(remainingSamples,toAdd)

    if(length(remainingSamples)==0){
      break
    }
  }

  toAdd<-c(toAdd,remainingSamples)

  return(toAdd)
}

## Exported
## Documentation Revised
CoRe.ADaM<-function(depMat,display=TRUE,
                    main_suffix='fitness genes in at least 1 cell line',
                    xlab='n. dependent cell lines',
                    ntrials=1000,verbose=TRUE,TruePositives){

  if(verbose){
    print('- Profiling of number of fitness genes across fixed numbers of cell lines and its cumulative sums')}
  pprofile<-CoRe.panessprofile(depMat=depMat,display = display,xlab = xlab,main_suffix = main_suffix)
  if(verbose){print('+ Done!')
    print('- Null modeling numbers of fitness genes across numbers of cell lines and their cumulative sums')}
  nullmodel<-CoRe.generateNullModel(depMat=depMat,ntrials = ntrials,verbose = verbose,display = display)
  if(verbose){print('+ Done!')
    print('- Computing empirical odds of numbers of fitness genes per number of cell lines')}
  EO<-CoRe.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM)
  if(verbose){print('+ Done')
    print('- Profiling true positive rates')}
  TPR<-CoRe.truePositiveRate(depMat,TruePositives)
  if(verbose){print('- Done!')
    print('+ Calculating ADaM threshold (min. n. of dependent cell lines for core fitness genes)')}
  crossoverpoint<-CoRe.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'curated BAGEL essential',display = display)
  if(verbose){print(paste('ADaM threshold =',crossoverpoint,'(out of',ncol(depMat),'cell lines)'))
    print('- Done!')
    print('+ Estimating set of core fitness genes')}
  coreFitnessGenes<-CoRe.coreFitnessGenes(depMat,crossoverpoint)
  if(verbose){print('- Done!')}
  return (coreFitnessGenes)
}

#--- Assemble expression based false positives
CoRe.AssembleFPs<-function(URL='https://ndownloader.figshare.com/files/26261476'){
  dir.create(tmp <- tempfile())
  dir.create(file.path(tmp, "mydir"))
  print('Downloading zipped CCLE expression data from DepMap portal')
  download.file(URL,file.path(tmp, "mydir","CCLE_expression.csv"))
  print('...done')

  print('Reading Expression data matrix...')
  X <- read.csv(file.path(tmp,'mydir','CCLE_expression.csv'),
                stringsAsFactors = FALSE,
                header=TRUE,
                row.names = 1)

  gnames<-rownames(X)
  clnames<-colnames(X)

  numdata<-as.matrix(X)

  numdata<-log2(numdata+1)
  numdata<-t(numdata)
  print('Done')
  print('Selecting overall lowly expressed genes...')

  LowlyExpr<-CoRe.FiPer(depMat = numdata,percentile = 0.9,display = FALSE)$cfgenes

  LowlyExpr<-strsplit(LowlyExpr,'[..]')

  LowlyExpr<-sort(unlist(lapply(LowlyExpr,function(x){x[1][1]})))

  print('Done')
  return(LowlyExpr)

}

CoRe.panessprofile<-function(depMat,display=TRUE,
                             main_suffix='fitness genes in at least 1 cell line',
                             xlab='n. dependent cell lines'){

    depMat<-depMat[which(rowSums(depMat)>0),]
    panessprof<-rep(0,ncol(depMat))
    names(panessprof)<-as.character(1:ncol(depMat))
    paness<-summary(as.factor(rowSums(depMat)),maxsum = length(unique(as.factor(rowSums(depMat)))))
    panessprof[as.character(names(paness))]<-paness

    CUMsums<-rev(cumsum(rev(panessprof)))

    names(CUMsums)<-paste('>=',names(CUMsums),sep='')

    if(display){
        par(mfrow=c(2,1))
        par(mar=c(6,4,4,1))

        main=paste(nrow(depMat),main_suffix)
        barplot(panessprof,ylab='n.genes',xlab=xlab,cex.axis = 0.8,cex.names = 0.8,
                las=2,main=main)

        barplot(CUMsums,ylab='n.genes',xlab=xlab,cex.axis = 0.8,cex.names = 0.6,
                las=2,main='Cumulative sums')
    }
    return(list(panessprof=panessprof,CUMsums=CUMsums))
}

CoRe.generateNullModel<-function(depMat,ntrials=1000,display=TRUE,verbose=TRUE){

    set.seed(100812)
    depMat<-depMat[which(rowSums(depMat)>0),]
    nullProf<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,1:ncol(depMat)))
    nullCumSUM<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,paste('≥',1:ncol(depMat),sep='')))
    if(verbose){
      print('Generating null model...')
      pb <- txtProgressBar(min=1,max=ntrials,style=3)
    }
    for (i in 1:ntrials){
      if(verbose){setTxtProgressBar(pb, i)}
        rMat<-
            CoRe.randomisedepMat(depMat)
        Ret<-
            CoRe.panessprofile(rMat,display = FALSE)
        nullProf[i,]<-Ret$panessprof
        nullCumSUM[i,]<-Ret$CUMsums
    }
    if(verbose){
      Sys.sleep(1)
      close(pb)
      print('')
      print('Done')
    }

    if (display){

        par(mfrow=c(2,1))
        main=c(paste(ntrials,' randomised essentiality profiles of\n',nrow(depMat),' genes across ',ncol(depMat),' cell lines',
                     sep=''))
        boxplot(nullProf,las=2,xlab='n. cell lines',ylab='fitness genes in n cell lines',main=main)

        colnames(nullCumSUM)<-paste(">=",1:ncol(nullCumSUM))
        boxplot(log10(nullCumSUM+1),las=2,main='Cumulative sums',xlab='n. cell lines',
                ylab='log10 [number of fitness genes + 1]',
                cex.axis=0.8)
    }

    return(list(nullProf=nullProf,nullCumSUM=nullCumSUM))
}

CoRe.randomisedepMat<-function(depMat){
    rmat<-apply(depMat,2,sample)
}

CoRe.empiricalOdds<-function(observedCumSum,simulatedCumSum){
  nsamples<-length(observedCumSum)
  ntrials<-nrow(simulatedCumSum)
  odds<-rep(NA,1,nsamples)
  names(odds)<-paste('≥',1:nsamples,sep='')
  for (i in 1:nsamples){
    PDF<-density(simulatedCumSum[,i])
    odds[i]<- log10(observedCumSum[i]/mean(simulatedCumSum[,i]))
  }
  return(odds)
}

CoRe.truePositiveRate<-function(depMat,essentialGeneSet){
  nsamples<-ncol(depMat)

  essentialGeneSet<-intersect(essentialGeneSet,rownames(depMat))

  TPR<-rep(NA,1,nsamples)
  names(TPR)<-paste('≥',1:nsamples,sep='')

  ncells<-rowSums(depMat)

  TP<-rep(NA,1,nsamples)
  names(TP)<-paste('≥',1:nsamples,sep='')

  P<-rep(NA,1,nsamples)
  names(P)<-paste('≥',1:nsamples,sep='')

  for (i in nsamples:1){
    positiveSet<-names(which(ncells>=i))
    P[i]<-length(positiveSet)
    truepositives<-intersect(positiveSet,essentialGeneSet)
    TP[i]<-length(truepositives)
    TPR[i]<-TP[i]/length(essentialGeneSet)
  }

  return(list(P=P,TP=TP,TPR=TPR))
}

CoRe.tradeoffEO_TPR<-function(EO,TPR,test_set_name,display=TRUE){
  x<-EO
  x[x==Inf]<-max(x[x<Inf])
  x<-(x-min(x))/(max(x)-min(x))

  y<-TPR
  y<-(y-min(y))/(max(y)-min(y))

  orEO<-EO
  orEO[orEO==Inf]<-max(orEO[orEO<Inf])
  orTPR<-TPR

  EO<-x
  TPR<-y
  point<-min(which(!y>x))

  if(display){
    CCOL<-'red'
    par(mar=c(4,4,4,4))
    par(mfrow=c(1,1))
    MAIN<-c('log10 (obs/Expct) n.genes [red, left]',
            paste('% covered ',test_set_name,' [blue, right]',sep=''))
    plot(EO,type='l',xlab='genes depleted in >= # cell lines',ylab='',axes=FALSE,lwd=4,main=MAIN,col=CCOL,cex.main=0.8,
         xlim=c(0,length(EO)))
    axis(2,at = seq(0,1,0.2),format(seq(min(orEO),max(orEO),(max(orEO)-min(orEO))/5),digits=2))
    axis(1)
    par(new=TRUE)
    plot(TPR,type='l',xlab='',ylab='',axes=FALSE,lwd=4,col='blue',ylim=c(0,1),xlim=c(0,length(EO)))
    axis(4,at = seq(0,1,0.2),format(seq(min(orTPR),max(orTPR),(max(orTPR)-min(orTPR))/5),digits=2))



    abline(v=point)
    abline(h=y[point],lty=2)

    points(point,y[point],pch=16,cex=2)

    legend('top',paste(format(100*orTPR[point],digits=2),'% covered',sep=''),bg = NULL,bty = 'n')
  }

  return(point)
}

CoRe.coreFitnessGenes<-function(depMat,crossoverpoint){
  coreFitnessGenes<-rownames(depMat)[rowSums(depMat)>=crossoverpoint]
  return (coreFitnessGenes)
}

## not Documented

#--- Downloading Binary Dependency Matrix (introduced in Behan 2019) from Project Score
CoRe.download_BinaryDepMatrix<-function(URL='https://cog.sanger.ac.uk/cmp/download/binaryDepScores.tsv.zip'){
  if(url.exists(URL)){
    temp <- tempfile()
    download.file(URL,temp)
    X <- read.table(unz(temp,'binaryDepScores.tsv'),stringsAsFactors = FALSE,sep='\t',header=TRUE,row.names = 1)
    colnames(X)<-str_replace_all(colnames(X),'[.]','-')
    colnames(X)<-unlist(lapply(colnames(X),function(x){
      if(str_sub(x,1,1)=='X'){
        x<-str_replace(x,'X','')
      }else{
        x
      }
    }))
  }else{
    X <- NULL
  }
  return(X)
}

## Non Documented

#--- Downloading Cell Passport models annotation file
CoRe.download_AnnotationModel<-function(URL='https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz'){
  if(url.exists(URL)){
    X <- read_csv(URL)
  }else{
    X <- NULL
  }
  return(X)
}

#--- Downloading and scaling Quantitative Dependency Matrix (introduced in Behan 2019) from Project Score
CoRe.download_DepMatrix<-function(URL='https://cog.sanger.ac.uk/cmp/download/essentiality_matrices.zip',
                                  scaled=FALSE,
                                  ess=NULL,
                                  noness=NULL){

  if(url.exists(URL)){
    dir.create(tmp <- tempfile())
    dir.create(file.path(tmp, "mydir"))
    print('Downloading zipped essentiality matrices from Project Score...')
    download.file(URL,file.path(tmp, "mydir","essentiality_matrices.zip"))
    print('...done')

    dir(file.path(tmp,'mydir'))

    print('Uncompressing zipped essentiality...')
    unzip(file.path(tmp,'mydir','essentiality_matrices.zip'),exdir = file.path(tmp,'mydir'))
    print('...done')

    print('Reading CRIPRcleanR corrected essentiality logFCs...')
    X <- read.table(file.path(tmp,'mydir','EssentialityMatrices','01_corrected_logFCs.tsv'),
                    stringsAsFactors = FALSE,
                    sep='\t',
                    header=TRUE,
                    row.names = 1)

    colnames(X)<-str_replace_all(colnames(X),'[.]','-')
    colnames(X)<-unlist(lapply(colnames(X),function(x){
      if(str_sub(x,1,1)=='X'){
        x<-str_replace(x,'X','')
      }else{
        x
      }
    }))
    print('...done')

    if(scaled){
      X<-scale_to_essentials(X,ess,noness)
    }

  }else{
    X <- NULL
  }
  return(X)
}

#--- Extracting Dependency SubMatrix for a given tissue or cancer type, among those included
#--- in the latest model annotation file on the cell model passports (cite Donny's paper and website URL)
CoRe.extract_tissueType_SubMatrix<-function(fullDepMat,tissue_type="Non-Small Cell Lung Carcinoma"){
  cmp<-read_csv('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz')
  cls<-cmp$model_name[which(cmp$tissue==tissue_type | cmp$cancer_type==tissue_type)]
  cls<-intersect(cls,colnames(fullDepMat))
  return(fullDepMat[,cls])
}

#--- Execute ADaM on tissue or cancer type specifc dependency submatrix
CoRe.CS_ADaM<-function(pancan_depMat,
                       tissue_ctype = 'Non-Small Cell Lung Carcinoma',
                       clannotation = NULL,
                       display=TRUE,
                       main_suffix='fitness genes in at least 1 cell line',
                       xlab='n. dependent cell lines',
                       ntrials=1000,verbose=TRUE,TruePositives){

  cls<-clannotation$model_name[clannotation$tissue==tissue_ctype | clannotation$cancer_type==tissue_ctype]
  cls<-intersect(colnames(pancan_depMat),cls)
  cs_depmat<-pancan_depMat[,cls]

  return(CoRe.ADaM(cs_depmat,display=display,
                   main_suffix = main_suffix,
                   xlab=xlab,
                   ntrials=ntrials,
                   verbose=verbose,TruePositives = TruePositives))
}

#--- Execute ADaM tissue by tissue then at the pancancer level to compute pancancer core fintess genes
CoRe.PanCancer_ADaM<-function(pancan_depMat,
                              tissues_ctypes,
                              clannotation = NULL,
                              display=TRUE,
                              ntrials=1000,verbose=TRUE,TruePositives){


  systematic_CS_ADaM_res<-lapply(tissues_ctypes,function(x){
      if(verbose){
        print(paste('Running ADaM for',x))
      }
      CoRe.CS_ADaM(pancan_depMat,tissue_ctype = x,
                   clannotation,display = display,
                   ntrials = ntrials,
                   verbose=verbose,TruePositives = TruePositives)
    })

  all_TS_CF_genes<-sort(unique(unlist(systematic_CS_ADaM_res)))

  TS_CF_matrix<-
    do.call(cbind,lapply(systematic_CS_ADaM_res,function(x){is.element(all_TS_CF_genes,x)+0}))

  rownames(TS_CF_matrix)<-all_TS_CF_genes
  colnames(TS_CF_matrix)<-tissues_ctypes

  CoRe.ADaM(depMat = TS_CF_matrix, display = display,
            main_suffix = 'genes predicted to be core fitness for at least 1 tissue/cancer-type',
            xlab = 'n. tissue/cancer-type',verbose = FALSE,TruePositives = TruePositives)


}

#--- Computes recall and other ROC indicators for identified core fitness genes
#--- with respect to pre-defined signatures of essential genes
CoRe.CF_Benchmark<-function(testedGenes,background,priorKnownSignatures,falsePositives,displayBar=TRUE){
  priorKnownSignatures<-lapply(priorKnownSignatures,function(x){intersect(x,background)})
  falsePositives<-intersect(falsePositives,background)

  memb<-do.call(rbind,lapply(priorKnownSignatures,function(x){
    is.element(testedGenes,x)+0
  }))

  colnames(memb)<-testedGenes

  totals<-unlist(lapply(priorKnownSignatures,function(x){
      length(intersect(x,background))
  }))

  memb<-HeuristicMutExSorting(memb)

  TPRs<-rowSums(memb)/totals[rownames(memb)]

  x<-rowSums(memb)
  N<-length(background)
  n<-length(testedGenes)
  k<-totals

  ps<-phyper(x-1,n,N-n,k,lower.tail=FALSE)[rownames(memb)]

  if(displayBar){
    pheatmap(memb,cluster_rows = FALSE,cluster_cols = FALSE,col=c('white','blue'),show_colnames = FALSE,
             main=paste(length(testedGenes),'core fitness genes'),legend = FALSE,
             width = 5,height = 3)

    par(mfrow=c(1,2))
    par(mar=c(5,1,0,1))
    barplot(100*rev(TPRs),horiz = TRUE,names.arg = NA,xlab='% covered genes',border=FALSE)
    abline(v=0)
    abline(v=c(20,40,60,80,100),lty=2,col='gray',lwd=2)
    ps[ps==0]<-min(ps[ps>0]/10)
    barplot(rev(-log10(ps)),horiz = TRUE,names.arg = NA,xlab='-log10 pval',border=FALSE,xlim=c(1,200),log = 'x')
    abline(v=1)
    abline(v=seq(10,200,20),lty=2,col='gray',lwd=2)
  }


  FPR<-length(intersect(falsePositives,testedGenes))/length(falsePositives)
  PPV<-length(intersect(testedGenes,unlist(priorKnownSignatures)))/length(testedGenes)

  TPRs<-TPRs[order(names(TPRs))]
  ps<-ps[order(names(ps))]

  return(list(TPRs=data.frame(Recall=TPRs,EnrichPval=ps),PPV=PPV,FPR=FPR))
}

#--- Calculate the Core Fitness genes using the  90th-percentile least dependent cell line from
#--- Quantative knockout screen dependency matrix.
CoRe.FiPer<-function(depMat,display=TRUE,percentile=0.9,method='AUC'){

  depMat<-as.matrix(depMat)

  rankCL<-t(apply(depMat,1,function(x){sx<-order(x)
                                       x<-match(1:length(x),sx)}))

  rownames(rankCL)<- rownames(depMat)
  colnames(rankCL)<- colnames(depMat)

  rankG<-apply(depMat,2,function(x){sx<-order(x)
                                    x<-match(1:length(x),sx)})

  rownames(rankG)<- rownames(depMat)
  colnames(rankG)<- colnames(depMat)

  CLnumber<-ncol(depMat)
  threshold = as.integer(CLnumber*percentile)

  nG<-nrow(rankG)
  nCL<-ncol(rankG)

  if(method=='fixed'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){rankG[x,match(threshold,rankCL[x,])]}))
    Label = "Gene rank in 90th perc. least dep cell line"
  }

  if(method=='average'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){mean(rankG[x,names(which(rankCL[x,]>=threshold))])}))
    Label = "Gene average rank in ≥ 90th perc. of least dep cell lines"
  }

  if(method=='slope'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      b<-as.data.frame(a)
      p<-lm(a ~ seq(1:nCL) , data=b)
      coef(p)[2]
    }))
    Label = "Slope of gene ranks across ranked dep cell lines"
  }

  if(method=='AUC'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      sum(a)
    }))
    Label = "AUC of gene ranks across ranked dep cell lines"
  }

  rownames(LeastDependentdf)<-rownames(rankG)
  doR <- density(LeastDependentdf, bw = "nrd0")

  if (display){
    par(mfrow=c(2,1))
    hist(LeastDependentdf,breaks = 100,xlab="Rank",main=Label)
    plot(doR,main="Gaussian Kernel Density Estimation",xlab="Rank")

  }

  localmin <- which(diff(-1*sign(diff(doR$y)))==-2)[1]+1
  myranks<- doR$x
  rankthreshold <- as.integer(myranks[localmin])+1

  if(display){
    abline(v=rankthreshold,col='red')
    legend('topleft',legend='Discriminative Threshold',cex=0.7,col='red',lty=1)
  }

  cfgenes <- rownames(LeastDependentdf)[which(LeastDependentdf < rankthreshold[1])]

  return(list(cfgenes=cfgenes,geneRanks=LeastDependentdf,LocalMinRank=rankthreshold[1]))
}

CoRe.VisCFness<-function(depMat,gene,percentile=0.9,posControl='RPL8',negControl='MAP2K1',method='fixed'){
  gg<-gene
  depMat<-as.matrix(depMat)

  rankCL<-t(apply(depMat,1,function(x){
      sx<-order(x)
      x<-match(1:length(x),sx)
      }))

  rownames(rankCL)<- rownames(depMat)
  colnames(rankCL)<- colnames(depMat)

  rankG<-apply(depMat,2,function(x){
    sx<-order(x)
    x<-match(1:length(x),sx)})

  rownames(rankG)<- rownames(depMat)
  colnames(rankG)<- colnames(depMat)

  nG<-nrow(rankG)
  nCL<-ncol(rankG)

  par(mfrow=c(1,2))

  plot(rankG[negControl,names(sort(rankCL[negControl,]))],
       col=rgb(150,0,0,alpha = 180,maxColorValue = 255),
       pch=16,ylim=c(0,nrow(depMat)),
       xlab='cell line dependency ranks',ylab='gene dependency rank in x dependant cell line')
  points(rankG[posControl,names(sort(rankCL[posControl,]))],pch=16,
         col=rgb(0,200,100,alpha = 180,maxColorValue = 255))

  points(rankG[gg,names(sort(rankCL[gg,]))],pch=16,
         col=rgb(0,0,100,alpha = 180,maxColorValue = 255))

  threshold = as.integer(nCL*percentile)

  if(method=='fixed'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){rankG[x,match(threshold,rankCL[x,])]}))
    Label = "Gene rank in 90th perc. least dep cell line"
  }

  if(method=='average'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){mean(rankG[x,names(which(rankCL[x,]>=threshold))])}))
    Label = "Gene average rank in ≥ 90th perc. of least dep cell lines"
  }

  if(method=='slope'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      b<-as.data.frame(a)
      p<-lm(a ~ seq(1:nCL) , data=b)
      coef(p)[2]
    }))
    Label = "Slope of gene ranks across ranked dep cell lines"
  }

  if(method=='AUC'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      sum(a)
    }))
    Label = "AUC of gene ranks across ranked dep cell lines"
  }

  cc<-c(rgb(0,200,100,alpha = 180,maxColorValue = 255),
        rgb(150,0,0,alpha = 180,maxColorValue = 255),
        rgb(0,0,100,alpha = 180,maxColorValue = 255))
  legend('topleft',legend=c(posControl,negControl,gg),col=cc,pch=16)

  names(LeastDependentdf)<-rownames(rankG)

  hist(LeastDependentdf,main=Label)
  abline(v = LeastDependentdf[negControl],lwd=4,col=rgb(150,0,0,alpha = 180,maxColorValue = 255))
  abline(v = LeastDependentdf[posControl],lwd=4,col=rgb(0,200,100,alpha = 180,maxColorValue = 255))
  abline(v = LeastDependentdf[gg],lwd=4,col=rgb(0,0,100,alpha = 180,maxColorValue = 255))

  doR <- density(LeastDependentdf, bw = "nrd0")

  localmin <- which(diff(-1*sign(diff(doR$y)))==-2)[1]+1
  myranks<- doR$x
  rankthreshold <- as.integer(myranks[localmin])+1

  abline(v = rankthreshold,
         lwd=3,
         col='darkgrey',
         lty=2)

}

