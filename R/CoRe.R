

## Exported
## Non Documented

## Documented
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
CoRe.AdAM<-function(depMat,display=TRUE,
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
     print('+ Calculating AdAM threshold (min. n. of dependent cell lines for core fitness genes)')}
     crossoverpoint<-CoRe.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'curated BAGEL essential',display = display)
     if(verbose){print(paste('AdAM threshold =',crossoverpoint,'(out of',ncol(depMat),'cell lines)'))
     print('- Done!')
     print('+ Estimating set of core fitness genes')}
     coreFitnessGenes<-CoRe.coreFitnessGenes(depMat,crossoverpoint)
     if(verbose){print('- Done!')}
     return (coreFitnessGenes)
     }


## Non Documented
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

#--- Extracting Binary Dependency SubMatrix for a given tissue or cancer type, among those included
#--- in the latest model annotation file on the cell model passports (cite Donny's paper and website URL)
CoRe.extract_tissueType_BinDepMatrix<-function(fullBinDepMat,tissue_type="Non-Small Cell Lung Carcinoma"){
  cmp<-read_csv('https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz')
  cls<-cmp$model_name[which(cmp$tissue==tissue_type | cmp$cancer_type==tissue_type)]
  cls<-intersect(cls,colnames(fullBinDepMat))
  return(fullBinDepMat[,cls])
}

#--- Execute AdAM on tissue or cancer type specifc dependency submatrix
CoRe.CS_AdAM<-function(pancan_depMat,
                       tissue_ctype = 'Non-Small Cell Lung Carcinoma',
                       clannotation = NULL,
                       display=TRUE,
                       main_suffix='fitness genes in at least 1 cell line',
                       xlab='n. dependent cell lines',
                       ntrials=1000,verbose=TRUE,TruePositives){

  cls<-clannotation$model_name[clannotation$tissue==tissue_ctype | clannotation$cancer_type==tissue_ctype]
  cls<-intersect(colnames(pancan_depMat),cls)
  cs_depmat<-pancan_depMat[,cls]

  return(CoRe.AdAM(cs_depmat,display=display,
                   main_suffix = main_suffix,
                   xlab=xlab,
                   ntrials=ntrials,
                   verbose=verbose,TruePositives = TruePositives))
}

#--- Execute AdAM tissue by tissue then at the pancancer level to compute pancancer core fintess genes
CoRe.PanCancer_AdAM<-function(pancan_depMat,
                              tissues_ctypes,
                              clannotation = NULL,
                              display=TRUE,
                              ntrials=1000,verbose=TRUE,TruePositives){


  systematic_CS_AdAM_res<-lapply(tissues_ctypes,function(x){
      if(verbose){
        print(paste('Running AdAM for',x))
      }
      CoRe.CS_AdAM(pancan_depMat,tissue_ctype = x,
                   clannotation,display = display,
                   ntrials = ntrials,
                   verbose=verbose,TruePositives = TruePositives)
    })

  all_TS_CF_genes<-sort(unique(unlist(systematic_CS_AdAM_res)))

  TS_CF_matrix<-
    do.call(cbind,lapply(systematic_CS_AdAM_res,function(x){is.element(all_TS_CF_genes,x)+0}))

  rownames(TS_CF_matrix)<-all_TS_CF_genes
  colnames(TS_CF_matrix)<-tissues_ctypes

  CoRe.AdAM(depMat = TS_CF_matrix,
            main_suffix = 'genes predicted to be core fitness for at least 1 tissue/cancer-type',
            xlab = 'n. tissue/cancer-type',verbose = FALSE,TruePositives = TruePositives)


}

#--- Computes recall and other ROC indicators for identified core fitness genes
#--- with respect to pre-defined signatures of essential genes
CoRe.CF_Benchmark<-function(testedGenes,background,priorKnownSignatures){

  memb<-do.call(rbind,lapply(priorKnownSignatures,function(x){
    is.element(testedGenes,x)
  }))

  colnames(memb)<-testedGenes

  totals<-unlist(lapply(priorKnownSignatures,function(x){
      length(intersect(x,background))
  }))

  memb<-MExMaS.HeuristicMutExSorting(memb)

  pheatmap(memb,cluster_rows = FALSE,cluster_cols = FALSE,col=c('white','blue'),show_colnames = FALSE,
           main=paste(length(testedGenes),'core fitness genes'),legend = FALSE,
           width = 5,height = 3)


  TPRs<-rowSums(memb)/totals[rownames(memb)]

  x<-rowSums(memb)
  N<-length(background)
  n<-length(testedGenes)
  k<-totals

  ps<-phyper(x-1,n,N-n,k,lower.tail=FALSE)[rownames(memb)]

  par(mfrow=c(1,2))
  par(mar=c(5,1,0,1))
  barplot(100*rev(TPRs),horiz = TRUE,names.arg = NA,xlab='% covered genes',border=FALSE)
  abline(v=0)
  abline(v=c(20,40,60,80,100),lty=2,col='gray',lwd=2)
  ps[ps==0]<-min(ps[ps>0]/10)
  barplot(rev(-log10(ps)),horiz = TRUE,names.arg = NA,xlab='-log10 pval',border=FALSE,xlim=c(1,200),log = 'x')
  abline(v=1)
  abline(v=seq(10,200,20),lty=2,col='gray',lwd=2)

  return(data.frame(Recall=TPRs,EnrichPval=ps))
}

rearrangeMatrix<-function(patterns,GENES){

  remainingSamples<-colnames(patterns)

  toAdd<-NULL

  for (g in GENES){
    remainingGenes<-setdiff(GENES,g)

    P1<-matrix(c(patterns[g,remainingSamples]),length(g),length(remainingSamples),dimnames = list(g,remainingSamples))
    P2<-matrix(c(patterns[remainingGenes,remainingSamples]),length(remainingGenes),length(remainingSamples),
               dimnames=list(remainingGenes,remainingSamples))

    if(length(remainingGenes)>1){
      DD<-colnames(P1)[order(P1-colSums(P2),decreasing=TRUE)]
    }else{
      DD<-colnames(P1)[order(P1-P2,decreasing=TRUE)]
    }

    toAdd<-c(toAdd,names(which(patterns[g,DD]>0)))
    remainingSamples<-setdiff(remainingSamples,toAdd)
    if(length(remainingSamples)==0){
      break
    }
  }

  toAdd<-c(toAdd,remainingSamples)

  return(toAdd)
}
#'
#' #' Calculate the Core Fitness genes using the  90th-percentile least dependent cell line from Quantative knockout screen dependency matrix.
#' #'
#' #' @description This function identifies the Core Fitness genes from a given Quantative knockout screen dependency matrix where each row is gene and each column the cell line. The function uses all the cell lines and identifies the genes that are essential in majority of the cell lines.
#' #' @usage ADAM2.PercentileCF(depMat,
#' #'              display=TRUE,
#' #'              percentile=0.9,
#' #'              prefix='PercentileMethod')
#' #' @param depMat Quantative knockout screen dependency matrix where rows are genes and columns are samples. A real number in position \emph{[i,j]} represents the strength of dependency which indicates the amaount of loss of fitness in the \emph{j}-th sample in case of the inactivation of the \emph{i}-th gene. Higher strength of dependency indicates higher probability of beign a core fitness gene. These values are used for ranking the genes in terms of their dependecy strength.
#' #' @param display Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
#' #' @param percentile percentage of the cell lines where the given gene should show depletion. The default value is 0.9 indicating 90-th percentile least dependent cell line.
#' #' @param prefix if the display is false the plots are generated in the working directory using the prefix.
#' #' @details This function implements the idea that if a gene is essential then it should fall in the top Z most depleted genes in at least 90% of cell lines. Here the threshold Z is calculated in a data driven way.
#' #' For a given gene, we can rank its gene effect score in each cell line, then arrange cell lines in order of increasing gene effect score for that gene. This creates a bimodal distribution of gene ranks in their 90th-percentile least depleted lines.
#' #' Z is choosen as the minimum density between the two normal distributions that are estimated from these ranks. All genes with rank less than this threshold in their 90th percentile cell lines are reported
#' #' @return A list of the following vectors:
#' #' \item{cfgenes}{Vector of number of genes that are core fitness genes}
#' #' \item{LeastDependent}{A dataframe where each row corresponds to a gene.There are two columns: \emph{Value} stores the rank of the gene at the \emph{N}-th percentile least dependent cell line and the \emph{Gene} stores the gene name}
#' #' \item{threshold}{The rank threshold for core fitness genes}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleSBFData)
#' #' results <- ADAM2.PercentileCF(depMat=exampleSBFData,display=TRUE)
#' #' cfgenes <- results$cfgenes
#' #' @keywords functions
#' #' @export
#' ADAM2.PercentileCF<-function(depMat,display=TRUE,percentile=0.9,prefix='PercentileMethod'){
#'
#'   mydata_transpose<-t(depMat)
#'   mydata_transpose_df <- as.data.frame(mydata_transpose)
#'   CLnames <-rownames(mydata_transpose_df)
#'   Genenames <- colnames(mydata_transpose_df)
#'   CLnumber <- length(CLnames)
#'   Genenumber <- length(Genenames)
#'
#'   rankCL <- data.frame(matrix(nrow=Genenumber,ncol=CLnumber))
#'   rankCL<-t(apply(-depMat,1,order))
#'   rownames(rankCL)<-Genenames
#'   colnames(rankCL)<-seq(1,CLnumber)
#'
#'   rankGene<- data.frame(matrix(nrow=CLnumber,ncol=Genenumber))
#'   rankGene<-t(apply(-mydata_transpose_df,1,rank))
#'   rownames(rankGene)<-CLnames
#'   colnames(rankGene)<-Genenames
#'
#'   threshold = as.integer(CLnumber*percentile)
#'   LeastDependent90df <- data.frame(matrix(nrow=Genenumber,ncol=2))
#'   count=1
#'   for (i in Genenames){
#'     LeastDependent90df[count,1] <- rankGene[rankCL[i,threshold],i]
#'     LeastDependent90df[count,2] <- i
#'     count=count+1
#'   }
#'   rownames(LeastDependent90df)<-Genenames
#'   colnames(LeastDependent90df)<-c("Value","Gene")
#'
#'   doR <- density(LeastDependent90df$Value, bw = "nrd0")
#'
#'   if (display){
#'     par(mfrow=c(2,1))
#'     hist(LeastDependent90df$Value,breaks = 100,xlab="Rank",main="Distribution of gene ranks in their 90th percentile least depleting lines")
#'     plot(doR,main="Gaussian Kernel Density Estimation",xlab="Rank")
#'   }
#'   else{
#'     filename<-paste0("./",prefix,".jpeg")
#'     jpeg(filename, width = 8, height = 8, units = 'in', res = 600)
#'     par(mfrow=c(2,1))
#'     hist(LeastDependent90df$Value,breaks = 100,xlab="Rank",main="Distribution of gene ranks in their 90th percentile least depleting lines")
#'     plot(doR,main="Gaussian Kernel Density Estimation",xlab="Rank")
#'     dev.off()
#'   }
#'
#'   localmin <- which(diff(-1*sign(diff(doR$y)))==-2)+1
#'   myranks<- doR$x
#'   rankthreshold <- as.integer(myranks[localmin])+1
#'
#'   cfgenes <- LeastDependent90df[which(LeastDependent90df$Value < rankthreshold[1]),]$Gene
#'   return(list(cfgenes=cfgenes,LeastDependent=LeastDependent90df,threshold=rankthreshold[1]))
#' }
#'
#'
#' #' Calculate the Core Fitness genes using the Average 90th-percentile least dependent cell line from Quantative knockout screen dependency matrix.
#' #'
#' #' @description This function identifies the Core Fitness genes from a given Quantative knockout screen dependency matrix where each row is gene and each column the cell line. The function uses all the cell lines and identifies the genes that are essential in majority of the cell lines.
#' #' @usage ADAM2.PercentileAverageCF(depMat,
#' #'                     display=TRUE,
#' #'                     percentile=0.9,
#' #'                     prefix='PercentileAverageMethod')
#' #' @param depMat Quantative knockout screen dependency matrix where rows are genes and columns are samples. A real number in position \emph{[i,j]} represents the strength of dependency which indicates the amaount of loss of fitness in the \emph{j}-th sample in case of the inactivation of the \emph{i}-th gene. Higher strength of dependency indicates higher probability of beign a core fitness gene. These values are used for ranking the genes in terms of their dependecy strength.
#' #' @param display Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
#' #' @param percentile percentage of the cell lines where the given gene should show depletion. The default value is 0.9 indicating least dependent 90th percentile cell line.
#' #' @param prefix if the display is false the plots are generated in the working directory using the prefix.
#' #' @details This function implements the idea that if a gene is essential then it should fall in the top Z most depleted genes in at least 90% of cell lines. Here the threshold Z is calculated in a data driven way.
#' #' For a given gene, we can rank its gene effect score in each cell line, then arrange cell lines in order of increasing gene effect score for that gene. The average ranks of the genes in the 90th percentile of least depleted genes are calculated.
#' #' Z is choosen as the minimum density between the two gaussian distributions that are estimated from these average rankings. All genes with average rank less than this threshold in their 90th percentile least depleted cell lines are reported.
#' #' @return A list of the following vectors:
#' #' \item{cfgenes}{Vector of number of genes that are core fitness genes}
#' #' \item{LeastDependent}{A dataframe where each row corresponds to a gene.There are two columns: \emph{Value} stores the average rank of the gene at the \emph{N}-th percentile least dependent cell lines and the \emph{Gene} stores the gene name}
#' #' \item{threshold}{The rank threshold for core fitness genes}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleSBFData)
#' #' results <- ADAM2.PercentileAverageCF(depMat=exampleSBFData)
#' #' cfgenes <- results$cfgenes
#' #' @keywords functions
#' #' @export
#' ADAM2.PercentileAverageCF<-function(depMat,display=TRUE,percentile=0.9,prefix='PercentileAverageMethod'){
#'
#'   mydata_transpose<-t(depMat)
#'   mydata_transpose_df <- as.data.frame(mydata_transpose)
#'   CLnames <-rownames(mydata_transpose_df)
#'   Genenames <- colnames(mydata_transpose_df)
#'   CLnumber <- length(CLnames)
#'   Genenumber <- length(Genenames)
#'
#'   rankCL <- data.frame(matrix(nrow=Genenumber,ncol=CLnumber))
#'   rankCL<-t(apply(-depMat,1,order))
#'   rownames(rankCL)<-Genenames
#'   colnames(rankCL)<-seq(1,CLnumber)
#'
#'   rankGene<- data.frame(matrix(nrow=CLnumber,ncol=Genenumber))
#'   rankGene<-t(apply(-mydata_transpose_df,1,rank))
#'   rownames(rankGene)<-CLnames
#'   colnames(rankGene)<-Genenames
#'
#'   threshold = as.integer(CLnumber*percentile)
#'   Dependentdf <- data.frame(matrix(nrow=Genenumber,ncol=2))
#'   count=1
#'   for (i in Genenames){
#'     Dependentdf[count,1] <- mean(rankGene[,i][unlist(rankCL[i,threshold:CLnumber])])
#'     Dependentdf[count,2] <- i
#'     count=count+1
#'   }
#'   rownames(Dependentdf)<-Genenames
#'   colnames(Dependentdf)<-c("Value","Gene")
#'
#'   doR <- density(Dependentdf$Value, bw = "nrd0")
#'
#'   if (display){
#'     par(mfrow=c(2,1))
#'     hist(Dependentdf$Value,breaks = 100,xlab="Rank",main="Distribution of gene ranks in their 90th percentile least depleting lines")
#'     plot(doR,main="Gaussian Kernel Density Estimation",xlab="Rank")
#'   }
#'   else{
#'     filename<-paste0("./",prefix,".jpeg")
#'     jpeg(filename, width = 8, height = 8, units = 'in', res = 600)
#'     par(mfrow=c(2,1))
#'     hist(Dependentdf$Value,breaks = 100,xlab="Rank",main="Distribution of gene ranks in their 90th percentile least depleting lines")
#'     plot(doR,main="Gaussian Kernel Density Estimation",xlab="Rank")
#'     dev.off()
#'   }
#'
#'   localmin <- which(diff(-1*sign(diff(doR$y)))==-2)+1
#'   myranks<- doR$x
#'   rankthreshold <- as.integer(myranks[localmin])+1
#'
#'   cfgenes <- Dependentdf[which(Dependentdf$Value < rankthreshold[1]),]$Gene
#'   return(list(cfgenes=cfgenes,LeastDependent=Dependentdf,threshold=rankthreshold[1]))
#' }
#'
#'
#' #' Calculate the Core Fitness genes using linear modeling of the ranks of genes in all cell lines ordered wrt to their gene score effects
#' #'
#' #' @description This function identifies the Core Fitness genes from a given Quantative knockout screen dependency matrix where each row is gene and each column the cell line. The function uses all the cell lines and identifies the genes that are essential in majority of the cell lines.
#' #' @usage ADAM2.SlopeCF(depMat,
#' #'                 display=TRUE,
#' #'                 prefix='SlopeMethod')
#' #' @param depMat Quantative knockout screen dependency matrix where rows are genes and columns are samples. A real number in position \emph{[i,j]} represents the strength of dependency which indicates the amount of loss of fitness in the \emph{j}-th sample in case of the inactivation of the \emph{i}-th gene. Higher strength of dependency indicates higher probability of beign a core fitness gene. These values are used for ranking the genes wrt to gene effect scores.
#' #' @param display Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
#' #' @param prefix if the display is false the plots are generated in the working directory using the prefix.
#' #' @details This function implements the idea that if a gene is essential then it should have ranked better in all cell lines including the least dependent cell lines. Instead of calculating rank threshold the ranks are modeled as a linear relation.
#' #' For a given gene, we can rank its gene effect score in each cell line, then arrange cell lines in order of increasing gene effect score for that gene. The ranks of the genes in these cell lines are fitted to a linear model where smaller slope values indicates higher dependency in all the cell lines.
#' #' A slope threshold is choosen as the minimum density between the two gaussian distributions that are estimated from the distribution of slopes. All genes with a slope less than this threshold is reported. Notice that we do not need to put a constraint such as 90th percentile least depleated cell lines.
#' #' @return A list of the following vectors:
#' #' \item{cfgenes}{Vector of number of genes that are core fitness genes}
#' #' \item{LeastDependent}{A dataframe where each row corresponds to a gene.There are two columns: \emph{Value} stores the slope of linear model that fits the rank of the gene in allcell lines and the \emph{Gene} stores the gene name}
#' #' \item{threshold}{The slope threshold for core fitness genes}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleSBFData)
#' #' results <- ADAM2.SlopeCF(depMat=exampleSBFData)
#' #' cfgenes <- results$cfgenes
#' #' @keywords functions
#' #' @export
#' ADAM2.SlopeCF<-function(depMat,display=TRUE,prefix='SlopeMethod'){
#'
#'   mydata_transpose<-t(depMat)
#'   mydata_transpose_df <- as.data.frame(mydata_transpose)
#'   CLnames <-rownames(mydata_transpose_df)
#'   Genenames <- colnames(mydata_transpose_df)
#'   CLnumber <- length(CLnames)
#'   Genenumber <- length(Genenames)
#'
#'   rankCL <- data.frame(matrix(nrow=Genenumber,ncol=CLnumber))
#'   rankCL<-t(apply(-depMat,1,order))
#'   rownames(rankCL)<-Genenames
#'   colnames(rankCL)<-seq(1,CLnumber)
#'
#'   rankGene<- data.frame(matrix(nrow=CLnumber,ncol=Genenumber))
#'   rankGene<-t(apply(-mydata_transpose_df,1,rank))
#'   rownames(rankGene)<-CLnames
#'   colnames(rankGene)<-Genenames
#'
#'   Dependentdf <- data.frame(matrix(nrow=Genenumber,ncol=2))
#'   count=1
#'
#'   for (i in Genenames){
#'     a <- rankGene[,i][unlist(rankCL[i,])]
#'     b<-as.data.frame(a)
#'     p<-lm(a ~ seq(1:CLnumber) , data=b)
#'     Dependentdf[count,1] <- coef(p)[2]
#'     Dependentdf[count,2] <- i
#'     count=count+1
#'   }
#'   rownames(Dependentdf)<-Genenames
#'   colnames(Dependentdf)<-c("Value","Gene")
#'
#'   doR <- density(Dependentdf$Value, bw = "nrd0")
#'
#'   if (display){
#'     par(mfrow=c(2,1))
#'     hist(Dependentdf$Value,breaks = 100,xlab="Slope",main="Distribution of slopes of the ranks of genes in all cell lines")
#'     plot(doR,main="Gaussian Kernel Density Estimation",xlab="Slope")
#'   }
#'   else{
#'     filename<-paste0("./",prefix,".jpeg")
#'     jpeg(filename, width = 8, height = 8, units = 'in', res = 600)
#'     par(mfrow=c(2,1))
#'     hist(Dependentdf$Value,breaks = 100,xlab="Slope",main="Distribution of slopes of the ranks of genes in all cell lines")
#'     plot(doR,main="Gaussian Kernel Density Estimation",xlab="Slope")
#'     dev.off()
#'   }
#'
#'   localmin <- which(diff(-1*sign(diff(doR$y)))==-2)+1
#'   myranks<- doR$x
#'   rankthreshold <- myranks[localmin+1]
#'
#'   cfgenes <- Dependentdf[which(Dependentdf$Value < rankthreshold[1]),]$Gene
#'   return(list(cfgenes=cfgenes,LeastDependent=Dependentdf,threshold=rankthreshold[1]))
#' }
#'
#' #' Calculate the Log transformed Bayesian Factor given a bimodal distribution of gene ranks or slopes of the linear fit of ranks.
#' #'
#' #' @description This function calculates the log transformed bayesian factor between essential and non-essential genes, by fitting two normal distributins to the bimodal distribution of the ranks of genes in their \emph{N}-th percentile least dependent cell lines ordered wrt gene effect score/estimated slopes from linear fitting of the ordered ranks of genes in cell lines. The function fits a two normal distributions on the bimodal distribution and estimates the mean the standard deviation of these distributions. Based on the estimations the bayesian factor is calculated for each gene using the probabilities of genes coming from essential or non-essential distributions.
#' #' @usage ADAM2.CalculateBayesianfactor(RankDistribution,
#' #'                 display=TRUE,
#' #'                 prefix='BayesianFactor')
#' #' @param RankDistribution Quantative knockout screen dependency matrix where rows are genes and columns are samples. A real number in position \emph{[i,j]} represents the strength of dependency which indicates the amaount of loss of fitness in the \emph{j}-th sample in case of the inactivation of the \emph{i}-th gene. Higher strength of dependency indicates higher probability of beign a core fitness gene. These values are used for ranking the genes in terms of their dependecy strength.
#' #' @param display Boolean, default is TRUE. Should plots of the bimodal normal distribution fitting
#' #' @param prefix if the display is false the plots are generated in the working directory using the prefix.
#' #' @details This function is using the mixdist R library to fit the given data into a distribution that can be represented as a mixture of two normal distributions. The mean and the standard deviations for each normal distribution is estimated
#' #' Each normal distribution corresponds to the essential and non-essential genes. For each gene a bayesian factor is calculated which can be defined as the Prob(gene|essential)/Prob(gene|non-essential). The log transformed bayesian factors are reported by this function.
#' #' @return A data frame with the following columns:
#' #' \item{Gene}{Gene name}
#' #' \item{logBF}{Log of the Bayesian Factor}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleSBFData)
#' #' results <- ADAM2.SlopeCF(depMat=exampleSBFData,display=TRUE)
#' #' bfresults <- ADAM2.CalculateBayesianfactor(RankDistribution=results$LeastDependent)
#' #' @keywords functions
#' #' @import mixdist
#' #' @export
#' ADAM2.CalculateBayesianfactor<-function(RankDistribution,display=TRUE,prefix='BayesianFactor'){
#'   his<-hist(RankDistribution$Value,breaks=100)
#'   df<-data.frame(mid=his$mids, cou=his$counts)
#'
#'   doR <- density(RankDistribution$Value, bw = "nrd0")
#'   localmax <- which(diff(sign(diff(doR$y)))==-2)+1
#'   myranks<- doR$x
#'   rankthreshold <- as.integer(myranks[localmax])+1
#'
#'   estimate_m1 <- rankthreshold[1]
#'   estimate_m2 <- tail(rankthreshold,n=1)
#'
#'   estimate_sd = estimate_m2-estimate_m1 / 5;
#'
#'
#'   fitpro <- mix(as.mixdata(df), mixparam(mu=c(estimate_m1,estimate_m2), sigma=c(estimate_sd,3*estimate_sd)), dist="norm")
#'
#'   if (display){
#'     par(mfrow=c(2,1))
#'     hist(RankDistribution$Value,breaks = 100,xlab="Value",main="Bimodal Distribution")
#'     plot(fitpro)
#'   }
#'   else{
#'     filename<-paste0("./",prefix,".jpeg")
#'     jpeg(filename, width = 8, height = 8, units = 'in', res = 600)
#'     par(mfrow=c(2,1))
#'     hist(RankDistribution$Value,breaks = 100,xlab="Value",main="Bimodal Distribution")
#'     plot(fitpro)
#'     dev.off()
#'   }
#'
#'   m1 <- fitpro$parameters$mu[1]
#'   m2 <- fitpro$parameters$mu[2]
#'
#'   s1 <- fitpro$parameters$sigma[1]
#'   s2 <- fitpro$parameters$sigma[2]
#'
#'   m<-matrix(data=RankDistribution$Value)
#'
#'   bak<-cbind(RankDistribution$Gene,apply(m,1, function(x) log(dnorm(x,mean=m1,sd=s1)/dnorm(x,mean=m2,sd=s2))))
#'   colnames(bak) <- c('Gene','logBF')
#'   bak_df <- as.data.frame(bak)
#'
#'   return(bak_df)
#' }


