

#' Calculate profile of number of fitness genes across fixed numbers of cell lines and cumulative sums.
#'
#' @description This function calculates the numbers (and cumulative numbers) of genes whose inactivation exerts a fitness effect in \emph{n} cell lines, varying \emph{n} from 1 to the number of cell lines in the dataset in input.
#' @usage ADAM2.panessprofile(depMat,
#'                    display=TRUE,
#'                    main_suffix='fitness genes in at least 1 cell line',
#'                    xlab='n. dependent cell lines')
#' @param depMat Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.
#' @param display Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
#' @param main_suffix If display=TRUE, title suffix to give to plot of number of genes depleted in a give number of cell lines, default is 'genes depleted in at least 1 cell line'
#' @param xlab If display=TRUE, label to give to x-axis of the plots, default is 'n. cell lines'
#' @return A list with the following two named vectors:
#' \item{panessprof}{Number of genes that are depleted for a number of cell lines}
#' \item{CUMsums}{Cumulative number of genes depleted in at least x cell lines}
#' @author C. Pacini, E. Karakoc & F. Iorio
#' @examples
#' data(exampleDepMat)
#' pprofile <- ADAM2.panessprofile(depMat = exampleDepMat)
#' @keywords functions
#' @export
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
#'
#' #' Generate null profile of number of fitness genes across fixed numbers of cell lines and cumulative sums.
#' #'
#' #' @description This function randomly perturbs the binary dependency matrix to generate a null distribution of profiles of fitness genes across fixed number of cell lines, and corresponding null distribution of cumulative sums.
#' #' @usage ADAM2.generateNullModel(depMat,
#' #'                        ntrials=1000,
#' #'                        display=TRUE)
#' #' @param depMat Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.
#' #' @param ntrials Integer, default =100. How many times to randomly perturb dependency matrix to generate the null distributions.
#' #' @param display Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
#' #' @details For a number of trials specified in (\code{ntrials}) the inputted binary dependency matrix is randomised, keeping its column marginal sums. The profiles of fitness genes across fixed number of cell lines, and corresponding cumulative sums, are returned for each random perturbation.
#' #' @return A list with the following two named vectors:
#' #' \item{nullProf}{Matrix of number of fitness genes for fixed number of cell lines from. Each rows of matrix corresponds to a random trial.}
#' #' \item{nullCumSum}{Matrix of profile of cumulative number of fitness genes in fixed number of cell lines. Each row of matrix is one random trial.}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleDepMat)
#' #' pprofile <- ADAM2.generateNullModel(depMat = exampleDepMat,ntrials=1000)
#' #' @keywords functions
#' #' @export
#' ADAM2.generateNullModel<-function(depMat,ntrials=1000,display=TRUE){
#'
#'     set.seed(100812)
#'     depMat<-depMat[which(rowSums(depMat)>0),]
#'     nullProf<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,1:ncol(depMat)))
#'     nullCumSUM<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,paste('≥',1:ncol(depMat),sep='')))
#'     print('Generating null model...')
#'     pb <- txtProgressBar(min=1,max=ntrials,style=3)
#'
#'     for (i in 1:ntrials){
#'         setTxtProgressBar(pb, i)
#'         rMat<-
#'             ADAM2.randomisedepMat(depMat)
#'         Ret<-
#'             ADAM2.panessprofile(rMat,display = FALSE)
#'         nullProf[i,]<-Ret$panessprof
#'         nullCumSUM[i,]<-Ret$CUMsums
#'     }
#'     Sys.sleep(1)
#'     close(pb)
#'     print('')
#'     print('Done')
#'
#'     if (display){
#'
#'         par(mfrow=c(2,1))
#'         main=c(paste(ntrials,' randomised essentiality profiles of\n',nrow(depMat),' genes across ',ncol(depMat),' cell lines',
#'                      sep=''))
#'         boxplot(nullProf,las=2,xlab='n. cell lines',ylab='genes depleted in n cell lines',main=main)
#'
#'         colnames(nullCumSUM)<-paste(">=",1:ncol(nullCumSUM))
#'         boxplot(log10(nullCumSUM+1),las=2,main='Cumulative sums',xlab='n. cell lines',
#'                 ylab='log10 [number of genes + 1]',
#'                 cex.axis=0.8)
#'     }
#'
#'     return(list(nullProf=nullProf,nullCumSUM=nullCumSUM))
#' }
#'
#' #' Binary matrix randomisation preserving column totals
#' #'
#' #' @description This function takes in input a matrix and shuffles its entries column wisely. If the matrix is binary then then matrix resulting from this shuffling will have the same column marginal totals of the inpputted one.
#' #' @usage ADAM2.randomisedepMat(depMat)
#' #' @param depMat A numerical matrix.
#' #' @return The matrix given in input with entries shuffled column wisely.
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleDepMat)
#' #' rnd_exampleDepMat<-ADAM2.randomisedepMat(exampleDepMat)
#' #' @keywords functions
#' #' @export
#' ADAM2.randomisedepMat<-function(depMat){
#'     rmat<-apply(depMat,2,sample)
#' }
#'
#' #' Set reference set of predefined essential genes
#' #'
#' #' @description This function takes in input a filename that contains the predefined essential sets that are used for calculating true positive rates
#' #' @usage ADAM2.setEssentialGenes(reffile)
#' #' @param reffile A text file that contains a gene name per line. These genes are predefined essential genes.
#' #' @return vector of predefined reference genes that can be used in ADAM.truePositiveRate function
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @keywords functions
#' #' @export
#' ADAM2.setEssentialGenes<-function(reffile){
#'   essentialGenes <- read.table(file=reffile,header=FALSE)
#' }
#'
#'
#' #' Empirical odds of number of fitness genes per number of cell lines
#' #'
#' #' @description This function calculates log10 odd ratios of observed vs. expected profiles of cumulative number of fitness genes in fixed number of cell lines. Expected values are the mean of those observed across randomised version of the observed binary matrix.
#' #' @usage ADAM2.empiricalOdds(observedCumSum,
#' #'                    simulatedCumSum)
#' #' @param observedCumSum Observed profile of cumulative sum of numbers of fitness genes in fixed number of cell lines. This is generated by the \code{ADAM.panessprofile} function.
#' #' @param simulatedCumSum Random profiles of cumulative sum of fitness genes in fixed number of cell lines. This is generated by the function \code{ADAM.generateNullModel}.
#' #' @return A named vector:
#' #' \item{odds}{log base 10 odd ratios of observed versus expected cumulative sums of number of fitness genes across fixed numbers of cell lines.}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleDepMat)
#' #' observed<-ADAM2.panessprofile(depMat=exampleDepMat)
#' #' null_m<-ADAM2.generateNullModel(depMat=exampleDepMat)
#' #' logOdds <- ADAM2.empiricalOdds(observedCumSum=observed$CUMsums,simulatedCumSum=null_m$nullCumSUM)
#' #' logOdds
#' #' @seealso ADAM2.panessprofile, ADAM2.generateNullModel
#' #' @keywords functions
#' #' @export
#' ADAM2.empiricalOdds<-function(observedCumSum,simulatedCumSum){
#'
#'     nsamples<-length(observedCumSum)
#'     ntrials<-nrow(simulatedCumSum)
#'
#'     odds<-rep(NA,1,nsamples)
#'     names(odds)<-paste('≥',1:nsamples,sep='')
#'     for (i in 1:nsamples){
#'
#'         PDF<-density(simulatedCumSum[,i])
#'
#'
#'         odds[i]<- log10(observedCumSum[i]/mean(simulatedCumSum[,i]))
#'
#'     }
#'     return(odds)
#' }
#'
#'
#' #' Profile of True Positive Rates
#' #'
#' #' @description This function calculates a profile of True Positive Rates for fitness genes in at least n cell lines, with positive cases from a reference set of essential genes.
#' #' @usage ADAM2.truePositiveRate(depMat,
#' #'                       essentialGeneSet)
#' #' @param depMat Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.
#' #' @param essentialGeneSet Reference set of predefined essential genes. This is used to define positive cases.
#' #' @details This function calculates true positive rates for fitness genes in at least n cell lines (for each n). First, this function calculates the number of cell lines each gene is a fitness gene. Second, for a given number of cell lines, the set of genes that are fitness genes in at least that number of cell lines is determined. Finally, this set of genes is then compared to the reference set of essential genes to calculate a true positive rate.
#' #' @return A list of the following vectors:
#' #' \item{P}{Vector of number of genes that are depleted for a number of cell lines.}
#' #' \item{TP}{Vector of number of genes in sets of P are true positives, i.e. in the essentialGeneSet.}
#' #' \item{TPR}{TP divided by number of genes in set essentialGeneSet to give the true positive rate.}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleDepMat)
#' #' pprofile<-ADAM2.panessprofile(depMat=exampleDepMat)
#' #' nullmodel<-ADAM2.generateNullModel(depMat=exampleDepMat,ntrials = 1000)
#' #' data(curated_BAGEL_essential)
#' #' EO<-ADAM2.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )
#' #' TPR<-ADAM2.truePositiveRate(exampleDepMat,curated_BAGEL_essential)
#' #' @keywords functions
#' #' @export
#' ADAM2.truePositiveRate<-function(depMat,essentialGeneSet){
#'     nsamples<-ncol(depMat)
#'
#'     essentialGeneSet<-intersect(essentialGeneSet,rownames(depMat))
#'
#'     TPR<-rep(NA,1,nsamples)
#'     names(TPR)<-paste('≥',1:nsamples,sep='')
#'
#'     ncells<-rowSums(depMat)
#'
#'     TP<-rep(NA,1,nsamples)
#'     names(TP)<-paste('≥',1:nsamples,sep='')
#'
#'     P<-rep(NA,1,nsamples)
#'     names(P)<-paste('≥',1:nsamples,sep='')
#'
#'     for (i in nsamples:1){
#'         positiveSet<-names(which(ncells>=i))
#'         P[i]<-length(positiveSet)
#'         truepositives<-intersect(positiveSet,essentialGeneSet)
#'         TP[i]<-length(truepositives)
#'         TPR[i]<-TP[i]/length(essentialGeneSet)
#'     }
#'
#'     return(list(P=P,TP=TP,TPR=TPR))
#' }
#'
#'
#' #' Calculate ADAM model threshold
#' #'
#' #' @description This function finds the minimum number of cell lines in which a gene needs to be fitness in order to be called core-fitness. This is defined as the \emph{n} providing the best trade-off between i) coverage of priori-known essential genes in the resulting set of predicted core-fitness genes, i.e. fitness in at least \emph{n} cell lines, and ii) deviance from expectation of the number of fitness genes in \emph{n} cell lines.
#' #' @usage ADAM2.tradeoffEO_TPR(EO,
#' #'                     TPR,
#' #'                     test_set_name)
#' #' @param EO Profile of empirical odds values. Computed with the \code{ADAM2.empiricalOdds} function.
#' #' @param TPR Profile of True positive rates for across number of cell line. Computed with the \code{ADAM2.truePositiveRate} function.
#' #' @param test_set_name Name to give to the analysis, used for plotting titles.
#' #' @details Compare and plot the log10 odds ratios with the true positive rates to find the cross over point where the true positive rate falls below the odds ratio.
#' #' @return ADAM model threshold:
#' #' \item{point}{Number of cell lines for which a gene needs to be a fitness gene in order to be predicted as core-fitness gene.}
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @seealso \code{\link{ADAM2.empiricalOdds}},
#' #' \code{\link{ADAM2.truePositiveRate}}
#' #' @examples
#' #' #load in example binary depletion matrix
#' #' data(exampleDepMat)
#' #'
#' #' # Generate the profiles of number of fitness genes across number of cell lines from
#' #' # observed data and corresponding comulative sums.
#' #' pprofile<-ADAM2.panessprofile(depMat=exampleDepMat)
#' #'
#' #' # Generate a set of random profiles of number of genes depleted for a number of cell lines
#' #' # and corresponding cumulative sums by perturbing observed data.
#' #' nullmodel<-ADAM2.generateNullModel(depMat=exampleDepMat,ntrials = 1000)
#' #'
#' #' #load a reference set of essential genes
#' #' data(curated_BAGEL_essential)
#' #'
#' #' # Calculate log10 odd ratios of observed/expected profiles of cumulative number of fitness
#' #' # genes in fixed number of cell lines.
#' #' # Observed values are from the ADAM.panessprofile function and expected are the average of
#' #' # random set from ADAM2.generateNullModel
#' #' EO<-ADAM2.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )
#' #'
#' #' # Calculate True positive rates for fitness genes in at least n cell lines in the observed
#' #' # dependency matrix, with positive cases from a reference set of essential genes
#' #' TPR<-ADAM2.truePositiveRate(exampleDepMat,curated_BAGEL_essential)
#' #'
#' #' # Calculate minimum number of cell lines a gene needs to be a fitness gene in order to
#' #' # be considered as a core-fitness gene
#' #' crossoverpoint<-ADAM2.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'curated BAGEL essential')
#' #' @keywords functions
#' #' @export
#' ADAM2.tradeoffEO_TPR<-function(EO,TPR,test_set_name){
#'
#'     CCOL<-'red'
#'
#'     x<-EO
#'     x[x==Inf]<-max(x[x<Inf])
#'     x<-(x-min(x))/(max(x)-min(x))
#'
#'     y<-TPR
#'     y<-(y-min(y))/(max(y)-min(y))
#'
#'     orEO<-EO
#'     orEO[orEO==Inf]<-max(orEO[orEO<Inf])
#'     orTPR<-TPR
#'
#'     EO<-x
#'     TPR<-y
#'     par(mar=c(4,4,4,4))
#'     MAIN<-c('log10 (obs/Expct) n.genes [red, left]',
#'             paste('% covered ',test_set_name,' [blue, right]',sep=''))
#'     plot(EO,type='l',xlab='genes depleted in >= # cell lines',ylab='',axes=FALSE,lwd=4,main=MAIN,col=CCOL,cex.main=0.8,
#'          xlim=c(0,length(EO)))
#'     axis(2,at = seq(0,1,0.2),format(seq(min(orEO),max(orEO),(max(orEO)-min(orEO))/5),digits=2))
#'     axis(1)
#'     par(new=TRUE)
#'     plot(TPR,type='l',xlab='',ylab='',axes=FALSE,lwd=4,col='blue',ylim=c(0,1),xlim=c(0,length(EO)))
#'     axis(4,at = seq(0,1,0.2),format(seq(min(orTPR),max(orTPR),(max(orTPR)-min(orTPR))/5),digits=2))
#'
#'
#'     point<-min(which(!y>x))
#'
#'     abline(v=point)
#'     abline(h=y[point],lty=2)
#'
#'     points(point,y[point],pch=16,cex=2)
#'
#'     legend('top',paste(format(100*orTPR[point],digits=2),'% covered',sep=''),bg = NULL,bty = 'n')
#'
#'     return(point)
#' }
#'
#'
#' #' Calculate the Core Fitness genes, starting from the binary dependency matrix.
#' #'
#' #' @description This function identifies the Core Fitness genes using the Adaptive Daisy Model (implemented ADAM) starting from a binary dependency matrix.
#' #' @usage ADAM2.coreFitnessGenesWrapper(depMat,
#' #'                       display=TRUE,
#' #'                       main_suffix='fitness genes in at least 1 cell line',
#' #'                       xlab='n. dependent cell lines',
#' #'                       ntrials=1000)
#' #' @param depMat Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.
#' #' @param display Boolean, default is TRUE. Should bar plots of the dependency profiles be plotted
#' #' @param main_suffix If display=TRUE, title suffix to give to plot of number of genes depleted in a give number of cell lines, default is 'genes depleted in at least 1 cell line'
#' #' @param xlab label to give to x-axis of the plots, default is 'n. cell lines'
#' #' @param ntrials Integer, default =1000. How many times to randomly perturb dependency matrix to generate the null distributions.
#' #' @details This function calculates the Core Fitness essential genes based on the calculated minimum number of cell lines that optimizes the True positive rates with log10 odds ratios. log10 odd ratios are calculated of observed vs. expected profiles of cumulative number of fitness genes in fixed number of cell lines. Expected values are the mean of those observed across randomised version of the observed binary matrix.
#' #' @return A vector that containing the Core Fitness Genes:
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleDepMat)
#' #' cfgenes <- ADAM2.coreFitnessGenesWrapper(depMat=exampleDepMat,ntrials=1000)
#' #' @keywords functions
#' #' @export
#' ADAM2.coreFitnessGenesWrapper<-function(depMat,display=TRUE,
#'                                 main_suffix='fitness genes in at least 1 cell line',
#'                                 xlab='n. dependent cell lines',
#'                                 ntrials=1000){
#'     pprofile<-ADAM2.panessprofile(depMat=depMat)
#'     nullmodel<-ADAM2.generateNullModel(depMat=depMat,ntrials = ntrials)
#'     EO<-ADAM2.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )
#'     TPR<-ADAM2.truePositiveRate(depMat,curated_BAGEL_essential)
#'     crossoverpoint<-ADAM2.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'curated BAGEL essential')
#'     coreFitnessGenes<-rownames(exampleDepMat)[rowSums(exampleDepMat)>=crossoverpoint]
#'     return (coreFitnessGenes)
#' }
#'
#' #' Calculate the Core Fitness genes given the binary dependency matrix and the minimal number of cell line threshold.
#' #'
#' #' @description This function identifies the Core Fitness genes using the Adaptive Daisy Model (implemented ADAM) starting from a binary dependency matrix.
#' #' @usage ADAM2.coreFitnessGenes(depMat,
#' #'                               crossoverpoint)
#' #' @param depMat Binary dependency matrix, rows are genes and columns are samples. 1 in position \emph{[i,j]} indicates that inactivation of the \emph{i}-th gene exerts a significant loss of fitness in the \emph{j}-th sample, 0 otherwise.
#' #' @param crossoverpoint minimum number of cell lines in which a gene needs to be fitness in order to be called core-fitness
#' #' @details This function calculates the Core Fitness essential genes based on the calculated minimum number of cell lines that optimizes the True positive rates with log10 odds ratios. log10 odd ratios are calculated of observed vs. expected profiles of cumulative number of fitness genes in fixed number of cell lines. Expected values are the mean of those observed across randomised version of the observed binary matrix.
#' #' @return A vector that containing the Core Fitness Genes:
#' #' @author C. Pacini, E. Karakoc & F. Iorio
#' #' @examples
#' #' data(exampleDepMat)
#' #' cfgenes <- ADAM2.coreFitnessGenes(depMat=exampleDepMat,crossoverpoint=3800)
#' #' @keywords functions
#' #' @export
#' ADAM2.coreFitnessGenes<-function(depMat,crossoverpoint){
#'   coreFitnessGenes<-rownames(depMat)[rowSums(depMat)>=crossoverpoint]
#'   return (coreFitnessGenes)
#' }
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


