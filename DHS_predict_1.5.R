library("snowfall")
library("randomForest")
library("e1071")
library("pROC")
library("ROCR")
library("ggplot2")
###0.require function for cross validation by RAP
###function CrossValidation() for randomfroest
###function CrossValidation1() for support vector machine
.cvSampleIndex <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}


.one_cross_validation <- function( cv, featureMat, positives, negatives, posSample_cv, negSample_cv) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]
  
  posLen <- length(train_genes_p)
  negLen <- length(train_genes_n)
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(train_genes_p, train_genes_n), ] )
  obj <- randomForest(x = fmat, y = factor(label))
  
  positives.train.score <- predict(obj, featureMat[train_genes_p,], type = "vote")[,"1"]
  negatives.train.score <- predict(obj, featureMat[train_genes_n,], type = "vote")[,"1"]
  positives.test.score <- predict(obj, featureMat[test_genes_p,], type = "vote")[,"1"]
  negatives.test.score <- predict(obj, featureMat[test_genes_n,], type = "vote")[,"1"]
  
  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                 positives.test = test_genes_p, negatives.test = test_genes_n,
                 positives.train.score = positives.train.score,
                 negatives.train.score = negatives.train.score,
                 positives.test.score = positives.test.score,
                 negatives.test.score = negatives.test.score,
                 train.AUC = train.AUC,
                 test.AUC = test.AUC) )
  
  res
}


CrossValidation <- function( seed = 1, featureMat, positives, negatives, cross = 10, cpus = 1){
  
  call <- match.call()
  
  #sample index for cv
  posSample_cv <- .cvSampleIndex(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".one_cross_validation", namespace = "RAP")
    sfLibrary( "pROC", character.only=TRUE)
    sfLibrary( "randomForest", character.only=TRUE )
    
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation( cv = j, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv)
    }
  }
  cvRes
}

.cvSampleIndex1 <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}


.one_cross_validation1 <- function( cv, featureMat, positives, negatives, posSample_cv, negSample_cv) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]
  
  posLen <- length(train_genes_p)
  negLen <- length(train_genes_n)
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(train_genes_p, train_genes_n), ] )
  obj <- svm(x = fmat, y = factor(label), probability = TRUE)
  
  positives.train.score <- attr(predict(obj, featureMat[train_genes_p,], probability = TRUE), "probabilities")[,"1"]
  negatives.train.score <- attr(predict(obj, featureMat[train_genes_n,], probability = TRUE), "probabilities")[,"1"]
  positives.test.score <- attr(predict(obj, featureMat[test_genes_p,],  probability = TRUE), "probabilities")[,"1"]
  negatives.test.score <- attr(predict(obj, featureMat[test_genes_n,], probability = TRUE), "probabilities")[,"1"]
  
  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                 positives.test = test_genes_p, negatives.test = test_genes_n,
                 positives.train.score = positives.train.score,
                 negatives.train.score = negatives.train.score,
                 positives.test.score = positives.test.score,
                 negatives.test.score = negatives.test.score,
                 train.AUC = train.AUC,
                 test.AUC = test.AUC) )
  
  res
}

CrossValidation1 <- function( seed = 1, featureMat, positives, negatives, cross = 10, cpus = 1){
  
  call <- match.call()
  
  #sample index for cv
  posSample_cv <- .cvSampleIndex1(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex1(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".one_cross_validation1", namespace = "RAP")
    sfLibrary( "pROC", character.only=TRUE)
    sfLibrary( "randomForest", character.only=TRUE )
    sfLibrary( "e1071", character.only=TRUE)
    
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation1, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation1( cv = j, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv)
    }
  }
  cvRes
}

plotROC <- function(cvRes) {
  
  
  cvListPredictions <- list()
  cvListLabels <- list()
  AUCVec <- rep(0, length(cvRes) )
  for( i in 1:length(cvRes) ) {
    curCV <- cvRes[[i]]
    cvListPredictions[[i]] <- c( curCV$positives.test.score, curCV$negatives.test.score )
    cvListLabels[[i]] <- c( rep(1, length(curCV$positives.test.score)), rep(0, length(curCV$negatives.test.score) ) )
    AUCVec[i] <- curCV$test.AUC
  }
  mAUC <- format( mean(AUCVec), digits= 3)
  
  #if( !require(ROCR) ) {
  #   install.packages("ROCR")
  #   library(ROCR)
  #}
  pred <- prediction(cvListPredictions, cvListLabels)
  perf <- performance(pred,"tpr","fpr")
  
  
  par(mar=c(5,6,4,2))   
  plot(perf, col= "gray", lty=3, main = paste( "AUC = ", mAUC, sep = ""), cex.lab = 2.5, cex.axis = 2, cex.main = 3, mgp = c(4,1.8,0) )
  plot(perf, col = "black",  lwd= 3, avg="vertical", spread.estimate="none", add=TRUE)  
  
}

##note: the current file need contain Pse-in-one scripts and independencies
##
##1.preprocess
##multiintersect
setwd("/home/malab1/wlt/03_grad_paper/02_workspace/randomforest/run/data/tissue")
system("bedtools multiinter -header -i  DHS_Ath_inflorescence_normal.bed DHS_Ath_leaf_normal.bed DHS_Ath_open_flower_normal.bed DHS_Ath_root_normal.bed DHS_Ath_seedling_normal.bed -empty -g AT10_chr.genome > DHSsmultiintersect.txt")
DHSs.multi <- read.delim(file="DHSsmultiintersect.txt",sep = "\t",header = T)
##preprocess of motif scan result
##run motifscan using fimo --thresh 1e-5 --o TAIR10 Ath_TF_binding_motifs.meme TAIR10_Chr.all.fasta
AT_motif_posi <- read.delim(file = "/home/malab1/wlt/03_grad_paper/01_data/619_motif_scan/result/fimo.gff",sep = "\t",header = F) 
AT_motif_posi <- AT_motif_posi[-1,]
AT_motif_posi <- AT_motif_posi[,c(1,4,5,7,9)]
colnames(AT_motif_posi) <- c("chr","start","end","strand","anno")
AT_motif_posi <- AT_motif_posi[order(AT_motif_posi$chr,AT_motif_posi$start),]
write.table(AT_motif_posi,file = "AT_motif_posi.bed",sep = "\t",quote = F,row.names = F,col.names = F)

##2.get positive sets
DHSpos_pri <-  DHSs.multi[DHSs.multi$num != 0,][,c(1:3)]
DHSpos_pri$len <- DHSpos_pri$end - DHSpos_pri$start
DHSpos_hub <- DHSpos_pri[DHSpos_pri$len==150,]
write.table(DHSpos_hub[,1:3],file = "DHSpos.bed",sep = "\t",quote = F,row.names = F,col.names = F)
system("bedtools getfasta -fi TAIR10_Chr.all.fasta -bed DHSpos.bed  -fo DHSpos.txt")
###remove non ATGC character manually, there is'nt non ATGC character in DHSpos.txt
###calculate positive sets feture
system("python kmer.py DHSpos.txt DHSpos_2mer.txt DNA -r 0  -k 2 -f tab")
system("python kmer.py DHSpos.txt DHSpos_RC2mer.txt DNA -r 1  -k 2 -f tab")
system("python pse.py DHSpos.txt DHSpos_PseDNC.txt DNA PC-PseDNC-General -lamada 3 -w 0.2 -e user_indices.txt -f tab")
system("bedtools coverage -a DHSpos.bed -b AT_motif_posi.bed  > pos_motif_cov.txt")

###get positive feature matrix
pos2mer <- as.matrix(read.table(file = "DHSpos_2mer.txt",sep = "\t",header = F, quote = ""))
posRC2mer <- as.matrix(read.table(file = "DHSpos_RC2mer.txt",sep = "\t",header = F, quote = ""))
posPseDNC <- as.matrix(read.table(file = "DHSpos_PseDNC.txt",sep = "\t",header = F))
pos_motif_cov <- as.matrix(read.delim(file = "pos_motif_cov.txt",sep = "\t",header = F)) 
aaa <- gsub(pattern = " ",replacement = "",x = pos_motif_cov[,4])  
class(aaa) <- "numeric"  
posfmat <- cbind(pos2mer,posRC2mer)
posfmat <- cbind(posfmat,posPseDNC)
posfmat <- cbind(posfmat,aaa)
colnames(posfmat) <- paste0("feature_", 1:ncol(posfmat))
DHSpos <- paste0("DHSs",c(1:nrow(posfmat)))
rownames(posfmat) <- DHSpos

##3.get negative sets
DHSneg_pri <-  DHSs.multi[DHSs.multi$num == 0,][,c(1:3)]
DHSneg_pri$length <- DHSneg_pri$end-DHSneg_pri$start
DHSneg_hub <- DHSneg_pri[DHSneg_pri$length>=2000,][,c(1:3)]
DHSneg_hub[DHSneg_hub[,2]==0,][,2] <- 1
## find midpoint
DHSneg_hub$res <- floor(apply(DHSneg_hub[,c(2:3)],1,mean))
DHSneg_hub$start1 <- DHSneg_hub$res-75
DHSneg_hub$end1   <- DHSneg_hub$res+75
DHSneg <- DHSneg_hub[,c(1,5,6)]
write.table(DHSneg,file = "DHSneg.bed",sep = "\t",quote = F,row.names = F,col.names = F)
system("bedtools getfasta -fi TAIR10_Chr.all.fasta -bed DHSneg.bed  -fo DHSneg.txt")
###remove non ATGC character manually, there is'nt non ATGC character in DHSneg.txt
###calculate negative sets feture
system("python kmer.py DHSneg.txt DHSneg_2mer.txt DNA -r 0  -k 2 -f tab")
system("python kmer.py DHSneg.txt DHSneg_RC2mer.txt DNA -r 1  -k 2 -f tab")
system("python pse.py DHSneg.txt DHSneg_PseDNC.txt DNA PC-PseDNC-General -lamada 3 -w 0.2 -e user_indices.txt -f tab")
system("bedtools coverage -a DHSneg.bed -b AT_motif_posi.bed  > neg_motif_cov.txt")
###get negative feature matrix
neg2mer <- as.matrix(read.table(file = "DHSneg_2mer.txt",sep = "\t",header = F, quote = ""))
negRC2mer <- as.matrix(read.table(file = "DHSneg_RC2mer.txt",sep = "\t",header = F, quote = ""))
negPseDNC <- as.matrix(read.table(file = "DHSneg_PseDNC.txt",sep = "\t",header = F))
neg_motif_cov <- read.delim(file = "neg_motif_cov.txt",sep = "\t",header = F) 
aaa1 <- gsub(pattern = " ",replacement = "",x = neg_motif_cov[,4])  
class(aaa1) <- "numeric"  
negfmat <- cbind(neg2mer,negRC2mer)
negfmat <- cbind(negfmat, negPseDNC)
negfmat <- cbind(negfmat, aaa1)

colnames(negfmat) <- paste0("feature_", 1:ncol(negfmat))
DHSneg <- paste0("notDHSs",c(1:nrow(negfmat)))
rownames(negfmat) <- DHSneg

fMat <-  rbind(posfmat,negfmat)
fMat  <- as.matrix(fMat)
class(fMat) <- "numeric"

##4.anlalysis features
## calculate positvie DHSs sets distribution
DHSpos_hub <- DHSpos_hub[,1:3]
DHSpos_hub$mid <- floor(apply(DHSpos_hub[,c(2:3)],1,mean))
DHSpos_hub.distance.lower <- DHSpos_hub[1:(nrow(DHSpos_hub)-1),c(1,4)]
DHSpos_hub.distance.higher <- DHSpos_hub[2:nrow(DHSpos_hub),c(1,4)]
DHSpos_hub.distance.lower$distance <- DHSpos_hub.distance.higher$mid-DHSpos_hub.distance.lower$mid
DHSpos_hub.distance <- DHSpos_hub.distance.lower[,c(1,3)]
colnames(DHSpos_hub.distance) <- c("chr","dist")
DHSpos_hub.distance1 <- DHSpos_hub.distance[-grep(pattern = "-",x = DHSpos_hub.distance[,2]),]
View(DHSpos_hub.distance1)
summary(DHSpos_hub.distance1$dist)
ggplot(data = DHSpos_hub.distance1, aes(x=chr,y=dist,color = chr,legend=F))+
  geom_boxplot()+
  ylim(100,3000)+
  labs(y="DHSs distance")
ggsave(filename = "~/wlt/03_grad_paper/DHSs_distribution.jpg")
###Evalutate different types of feature 

ttfMat <- as.data.frame(fMat)
ttfMat$type <- "positive"
ttfMat[48358:nrow(ttfMat),47] <- "negative"

aafMat1 <- data.frame(ttfMat[,c(1,47)])
aafMat1$feature <- "2mer_AA"
colnames(aafMat1) <- c("value","type","feature")

aafMat2 <- data.frame(ttfMat[,c(2,47)])
aafMat2$feature <- "2mer_AC"
colnames(aafMat2) <- c("value","type","feature")

aafMat3 <- data.frame(ttfMat[,c(3,47)])
aafMat3$feature <- "2mer_AG"
colnames(aafMat3) <- c("value","type","feature")

aafMat4 <- data.frame(ttfMat[,c(4,47)])
aafMat4$feature <- "2mer_AT"
colnames(aafMat4) <- c("value","type","feature")

aafMat5 <- data.frame(ttfMat[,c(5,47)])
aafMat5$feature <- "2mer_CA"
colnames(aafMat5) <- c("value","type","feature")

aafMat6 <- data.frame(ttfMat[,c(6,47)])
aafMat6$feature <- "2mer_CC"
colnames(aafMat6) <- c("value","type","feature")

aafMat7 <- data.frame(ttfMat[,c(7,47)])
aafMat7$feature <- "2mer_CG"
colnames(aafMat7) <- c("value","type","feature")

aafMat8 <- data.frame(ttfMat[,c(8,47)])
aafMat8$feature <- "2mer_CT"
colnames(aafMat8) <- c("value","type","feature")

aafMat9 <- data.frame(ttfMat[,c(9,47)])
aafMat9$feature <- "2mer_GA"
colnames(aafMat9) <- c("value","type","feature")

aafMat10 <- data.frame(ttfMat[,c(10,47)])
aafMat10$feature <- "2mer_GC"
colnames(aafMat10) <- c("value","type","feature")

aafMat11 <- data.frame(ttfMat[,c(11,47)])
aafMat11$feature <- "2mer_GG"
colnames(aafMat11) <- c("value","type","feature")

aafMat12 <- data.frame(ttfMat[,c(12,47)])
aafMat12$feature <- "2mer_GT"
colnames(aafMat12) <- c("value","type","feature")

aafMat13 <- data.frame(ttfMat[,c(13,47)])
aafMat13$feature <- "2mer_TA"
colnames(aafMat13) <- c("value","type","feature")

aafMat14 <- data.frame(ttfMat[,c(14,47)])
aafMat14$feature <- "2mer_TC"
colnames(aafMat14) <- c("value","type","feature")

aafMat15 <- data.frame(ttfMat[,c(15,47)])
aafMat15$feature <- "2mer_TG"
colnames(aafMat15) <- c("value","type","feature")

aafMat16 <- data.frame(ttfMat[,c(16,47)])
aafMat16$feature <- "2mer_TT"
colnames(aafMat16) <- c("value","type","feature")

aafMat17 <- data.frame(ttfMat[,c(17,47)])
aafMat17$feature <- "RC2mer_AA"
colnames(aafMat17) <- c("value","type","feature")

aafMat18 <- data.frame(ttfMat[,c(18,47)])
aafMat18$feature <- "RC2mer_AC"
colnames(aafMat18) <- c("value","type","feature")

aafMat19 <- data.frame(ttfMat[,c(19,47)])
aafMat19$feature <- "RC2mer_AG"
colnames(aafMat19) <- c("value","type","feature")

aafMat20 <- data.frame(ttfMat[,c(20,47)])
aafMat20$feature <- "RC2mer_AT"
colnames(aafMat20) <- c("value","type","feature")

aafMat21 <- data.frame(ttfMat[,c(21,47)])
aafMat21$feature <- "RC2mer_CA"
colnames(aafMat21) <- c("value","type","feature")

aafMat22 <- data.frame(ttfMat[,c(22,47)])
aafMat22$feature <- "RC2mer_CC"
colnames(aafMat22) <- c("value","type","feature")

aafMat23 <- data.frame(ttfMat[,c(23,47)])
aafMat23$feature <- "RC2mer_CG"
colnames(aafMat23) <- c("value","type","feature")

aafMat24 <- data.frame(ttfMat[,c(24,47)])
aafMat24$feature <- "RC2mer_GA"
colnames(aafMat24) <- c("value","type","feature")

aafMat25 <- data.frame(ttfMat[,c(25,47)])
aafMat25$feature <- "RC2mer_GC"
colnames(aafMat25) <- c("value","type","feature")

aafMat26 <- data.frame(ttfMat[,c(26,47)])
aafMat26$feature <- "RC2mer_TA"
colnames(aafMat26) <- c("value","type","feature")

aafMat27 <- data.frame(ttfMat[,c(27,47)])
aafMat27$feature <- "PseDNC_1"
colnames(aafMat27) <- c("value","type","feature")

aafMat28 <- data.frame(ttfMat[,c(28,47)])
aafMat28$feature <- "PseDNC_2"
colnames(aafMat28) <- c("value","type","feature")

aafMat29 <- data.frame(ttfMat[,c(29,47)])
aafMat29$feature <- "PseDNC_3"
colnames(aafMat29) <- c("value","type","feature")

aafMat30 <- data.frame(ttfMat[,c(30,47)])
aafMat30$feature <- "PseDNC_4"
colnames(aafMat30) <- c("value","type","feature")

aafMat31 <- data.frame(ttfMat[,c(31,47)])
aafMat31$feature <- "PseDNC_5"
colnames(aafMat31) <- c("value","type","feature")

aafMat32 <- data.frame(ttfMat[,c(32,47)])
aafMat32$feature <- "PseDNC_6"
colnames(aafMat32) <- c("value","type","feature")

aafMat33 <- data.frame(ttfMat[,c(33,47)])
aafMat33$feature <- "PseDNC_7"
colnames(aafMat33) <- c("value","type","feature")

aafMat34 <- data.frame(ttfMat[,c(34,47)])
aafMat34$feature <- "PseDNC_8"
colnames(aafMat34) <- c("value","type","feature")

aafMat35 <- data.frame(ttfMat[,c(35,47)])
aafMat35$feature <- "PseDNC_9"
colnames(aafMat35) <- c("value","type","feature")

aafMat36 <- data.frame(ttfMat[,c(36,47)])
aafMat36$feature <- "PseDNC_10"
colnames(aafMat36) <- c("value","type","feature")

aafMat37 <- data.frame(ttfMat[,c(37,47)])
aafMat37$feature <- "PseDNC_11"
colnames(aafMat37) <- c("value","type","feature")

aafMat38 <- data.frame(ttfMat[,c(38,47)])
aafMat38$feature <- "PseDNC_12"
colnames(aafMat38) <- c("value","type","feature")

aafMat39 <- data.frame(ttfMat[,c(39,47)])
aafMat39$feature <- "PseDNC_13"
colnames(aafMat39) <- c("value","type","feature")

aafMat40 <- data.frame(ttfMat[,c(40,47)])
aafMat40$feature <- "PseDNC_14"
colnames(aafMat40) <- c("value","type","feature")

aafMat41 <- data.frame(ttfMat[,c(41,47)])
aafMat41$feature <- "PseDNC_15"
colnames(aafMat41) <- c("value","type","feature")

aafMat42 <- data.frame(ttfMat[,c(42,47)])
aafMat42$feature <- "PseDNC_16"
colnames(aafMat42) <- c("value","type","feature")

aafMat43 <- data.frame(ttfMat[,c(43,47)])
aafMat43$feature <- "PseDNC_17"
colnames(aafMat43) <- c("value","type","feature")

aafMat44 <- data.frame(ttfMat[,c(44,47)])
aafMat44$feature <- "PseDNC_18"
colnames(aafMat44) <- c("value","type","feature")

aafMat45 <- data.frame(ttfMat[,c(45,47)])
aafMat45$feature <- "PseDNC_19"
colnames(aafMat45) <- c("value","type","feature")

aafMat46 <- data.frame(ttfMat[,c(46,47)])
aafMat46$feature <- "motif_count"
colnames(aafMat46) <- c("value","type","feature")
tt <- paste0("aafMat",1:46)
fMat_plot <- rbind(aafMat1,aafMat2,aafMat3,aafMat4,
                   aafMat5,aafMat6,aafMat7,aafMat8,
                   aafMat9,aafMat10,aafMat11,aafMat12,
                   aafMat13,aafMat14,aafMat15,aafMat16,
                   aafMat17,aafMat18,aafMat19,aafMat20,
                   aafMat21,aafMat22,aafMat23,aafMat24,
                   aafMat25,aafMat26,aafMat27,aafMat28,
                   aafMat29,aafMat30,aafMat31,aafMat32,
                   aafMat33,aafMat34,aafMat35,aafMat36,
                   aafMat37,aafMat38,aafMat39,aafMat40,
                   aafMat41,aafMat42,aafMat43,aafMat44,
                   aafMat45)
write.table(fMat_plot, file = "fMat_plot.txt" ,sep = "\t",quote = F,row.names = F,col.names = F)
ggplot(data=fMat_plot,aes(x=feature,y=value,colour=type))+
  geom_boxplot()+
  scale_x_discrete(labels=paste0("f",1:45))+
  labs(y="feture_value")
ggsave("~/wlt/03_grad_paper/seqfeaure_distribution.jpg",width = 20,height=8,dpi=600)

# ggplot(data = ttfMat, aes(x=feature_46, colour=type) )+
#   geom_density()+
#   xlim(0,50)
# ggsave("~/wlt/03_grad_paper/motiffeaure_distribution.jpg")

summary(posfmat[,46])
summary(negfmat[,46])


##5.Build and Evaluate DHSs predction model


Res_randomforest <- CrossValidation( seed = 3, featureMat=fMat, positives=DHSpos, negatives=DHSneg, cross = 5, cpus = 3)
Res_svm <- CrossValidation1( seed = 3, featureMat=fMat, positives=DHSpos, negatives=DHSneg, cross = 5, cpus = 3)
par()
# plotROC(Res_randomforest)
# plotROC(Res_svm)

cvRes1 <- Res_randomforest
cvRes2 <- Res_svm
# randomForest
cvListPredictions1 <- list()
cvListLabels1 <- list()
AUCVec1 <- rep(0, length(cvRes1) )
for( i in 1:length(cvRes1) ) {
  curCV1 <- cvRes1[[i]]
  cvListPredictions1[[i]] <- c( curCV1$positives.test.score, curCV1$negatives.test.score )
  cvListLabels1[[i]] <- c( rep(1, length(curCV1$positives.test.score)), rep(0, length(curCV1$negatives.test.score) ) )
  AUCVec1[i] <- curCV1$test.AUC
}
mAUC1 <- format( mean(AUCVec1), digits= 3)

#if( !require(ROCR) ) {
#   install.packages("ROCR")
#   library(ROCR)
#}
pred1 <- prediction(cvListPredictions1, cvListLabels1)
perf1 <- performance(pred1,"tpr","fpr")
###SVM
cvListPredictions2 <- list()
cvListLabels2 <- list()
AUCVec2 <- rep(0, length(cvRes2) )
for( i in 1:length(cvRes2) ) {
  curCV2 <- cvRes2[[i]]
  cvListPredictions2[[i]] <- c( curCV2$positives.test.score, curCV2$negatives.test.score )
  cvListLabels2[[i]] <- c( rep(1, length(curCV2$positives.test.score)), rep(0, length(curCV2$negatives.test.score) ) )
  AUCVec2[i] <- curCV2$test.AUC
}
mAUC2 <- format( mean(AUCVec2), digits= 3)

#if( !require(ROCR) ) {
#   install.packages("ROCR")
#   library(ROCR)
#}
pred2 <- prediction(cvListPredictions1, cvListLabels1)
perf2 <- performance(pred2,"tpr","fpr")

pdf(file = "/home/malab1/wlt/03_grad_paper/ROC.pdf", height = 5, width = 15)
par(mfrow = c(1,2))
plot(perf1, col= "gray", lty=3, type = "n", main = paste( "RFAUC = ", mAUC1, sep = ""), cex.lab = 2.5, cex.axis = 2, cex.main = 3, mgp = c(4,1.8,0) )
plot(perf1, col = "black",  lwd= 3, avg="vertical", spread.estimate="none", add=TRUE)

plot(perf2, col= "gray", lty=3, main = paste( "SVMAUC = ", mAUC2, sep = ""), cex.lab = 2.5, cex.axis = 2, cex.main = 3, mgp = c(4,1.8,0))
plot(perf2, col = "black",  lwd= 3, avg="vertical", spread.estimate="none", add=TRUE)
dev.off()
