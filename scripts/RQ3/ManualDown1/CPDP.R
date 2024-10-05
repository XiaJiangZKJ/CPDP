packages <- c("beanplot", "car", "class", "doSNOW", "e1071", "effsize", "foreach", "FSelector", "hmeasure", "igraph", "kernlab", "kknn", "lmtest", "markovchain", "moments", "monmlp", "neuralnet", "nnet", "perturb", "plyr", "PMCMR", "psych", "randomForest", "rms", "rpart", "RSNNS", "RWeka", "tree", "usdm")
packages <- packages[!(packages %in% installed.packages())]
if (length(packages) != 0) { install.packages(packages, repos="http://cran.rstudio.com/") }

ComputeAUC <- function(data, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE, allpercentcutoff=TRUE)
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  auc <- NULL
  
  cuts <- nrow(data)
  
  sumTRUE  <- sum(data$REL)  ### sum of true instances
  sumFALSE <- cuts - sumTRUE ### sum of false instances
  
  if (sumTRUE==0) { cat("sum of REL is zero"); auc <- -1 }
  
  TPR.old <- cumsum(data$REL) / sumTRUE
  FPR.old <- cumsum(1 - data$REL) / sumFALSE ### FPR.old <- (c(1:cuts) - cumsum(data$REL)) / sumFALSE
  
  TPR <- TPR.old
  FPR <- FPR.old
  if (!allpercentcutoff) {
    pos.equal <- which(data[-cuts, "PRE"]==data[-1, "PRE"])
    if (length(pos.equal)!=0) {
      TPR <- TPR.old[-pos.equal]
      FPR <- FPR.old[-pos.equal]
    }
  }
  
  cuts <- length(TPR)
  subareas <- rep(0, cuts)
  if (cuts==1) {
    subareas[1] <- 0.5
  } else {
    subareas[1]  <- TPR[1] * FPR[1] * 0.5
    subareas[-1] <- (TPR[-1]+TPR[-cuts]) * (abs(FPR[-1]-FPR[-cuts])) * 0.5
  }
  
  auc <- sum(subareas)
  
  if (is.na(auc)) { cat("AUC is NA"); auc <- -1 }
  
  names(auc) <- c("AUC")
  return(auc)
}

SortData <- function(data, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
  if (!all(c("NUM", "REL", "LOC", "PRE") %in% colnames(data))) { stop("ERROR: NUM, REL, LOC, or PRE is not colnames of data frame") }
  
  key.2nd <- "REL"
  if (effortaware) {
    key.2nd <- "density"
    if (!(key.2nd %in% colnames(data))) {
      data[[key.2nd]] <- data$NUM/(data$LOC+1)
    }
    sorted <- FALSE
  }
  
  if (!sorted) {
    if (worstcase) {
      data <- data[order(-data$PRE, +data[[key.2nd]], -data$LOC), ]
    } else if (bestcase) {
      data <- data[order(-data$PRE, -data[[key.2nd]], +data$LOC), ]
    } else if (LOCUP) {
      data <- data[order(-data$PRE, +data$LOC), ]
    } else {
      data <- data[order(-data$PRE), ]
    }
    sorted <- TRUE
  }
  
  return(data)
}

ComputeCPs <- ComputeAllClassificationPerformances <- function(train.dt, valid.dt, thresholds=NULL, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE, allpercentcutoff=TRUE)
{
  Recalls <- Precisions <- PDs <- PFs <- Correctness <- Completeness <- F1s <- F2s<- G1s <- G2s <- G3s <- Balances <- EDs <- MCCs <- Type1s <- Type2s <- Overalls <- NECM_15s <- NECM_20s <- NECM_25s <- NULL
  
  ### FUN <- match.fun(FUNCNAME)
  FUN <- ComputeCP
  if (!sorted) {
    train.dt <- SortData(data=train.dt, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    valid.dt <- SortData(data=valid.dt, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  if (is.null(thresholds)) {
    thresholds <- ComputeThreshold(data=train.dt, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)
  }
  
  thresh.list.names <- c("V.BPP", "V.BCE", "V.MFM")
  for (thresh.list.name in thresh.list.names) {
    temp <- FUN(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, threshold=thresholds[[thresh.list.name]])
    Recalls[[thresh.list.name]] <- temp$RECALL
    Precisions[[thresh.list.name]] <- temp$PRECISION
    PDs[[thresh.list.name]] <- temp$PD
    PFs[[thresh.list.name]] <- temp$PF
    Correctness[[thresh.list.name]] <- temp$Correctness
    Completeness[[thresh.list.name]] <- temp$Completeness
    F1s[[thresh.list.name]] <- temp$F1
    F2s[[thresh.list.name]] <- temp$F2
    G1s[[thresh.list.name]] <- temp$G1
    G2s[[thresh.list.name]] <- temp$G2
    G3s[[thresh.list.name]] <- temp$G3
    Balances[[thresh.list.name]] <- temp$Balance
    EDs[[thresh.list.name]] <- temp$ED
    MCCs[[thresh.list.name]] <- temp$MCC
    Type1s[[thresh.list.name]] <- temp$Type1
    Type2s[[thresh.list.name]] <- temp$Type2
    Overalls[[thresh.list.name]] <- temp$Overall
    NECM_15s[[thresh.list.name]] <- temp$NECM_15
    NECM_20s[[thresh.list.name]] <- temp$NECM_20
    NECM_25s[[thresh.list.name]] <- temp$NECM_25
  }
  
  ### cutoff value of 0.5
  temp <- FUN(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, threshold=0.5)
  Recalls[["V.0.5"]] <- temp$RECALL
  Precisions[["V.0.5"]] <- temp$PRECISION
  PDs[["V.0.5"]] <- temp$PD
  PFs[["V.0.5"]] <- temp$PF
  Correctness[["V.0.5"]] <- temp$Correctness
  Completeness[["V.0.5"]] <- temp$Completeness
  F1s[["V.0.5"]] <- temp$F1
  F2s[["V.0.5"]] <- temp$F2
  G1s[["V.0.5"]] <- temp$G1
  G2s[["V.0.5"]] <- temp$G2
  G3s[["V.0.5"]] <- temp$G3
  Balances[["V.0.5"]] <- temp$Balance
  EDs[["V.0.5"]] <- temp$ED
  MCCs[["V.0.5"]] <- temp$MCC
  Type1s[["V.0.5"]] <- temp$Type1
  Type2s[["V.0.5"]] <- temp$Type2
  Overalls[["V.0.5"]] <- temp$Overall
  NECM_15s[["V.0.5"]] <- temp$NECM_15
  NECM_20s[["V.0.5"]] <- temp$NECM_20
  NECM_25s[["V.0.5"]] <- temp$NECM_25
  
  cutoff.list.names <- c("P.BPP", "P.BCE", "P.MFM")
  for (cutoff.list.name in cutoff.list.names) {
    temp <- FUN(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, cutoff=thresholds[[cutoff.list.name]])
    Recalls[[thresh.list.name]] <- temp$RECALL
    Precisions[[thresh.list.name]] <- temp$PRECISION
    PDs[[thresh.list.name]] <- temp$PD
    PFs[[thresh.list.name]] <- temp$PF
    Correctness[[thresh.list.name]] <- temp$Correctness
    Completeness[[thresh.list.name]] <- temp$Completeness
    F1s[[thresh.list.name]] <- temp$F1
    F2s[[thresh.list.name]] <- temp$F2
    G1s[[thresh.list.name]] <- temp$G1
    G2s[[thresh.list.name]] <- temp$G2
    G3s[[thresh.list.name]] <- temp$G3
    Balances[[thresh.list.name]] <- temp$Balance
    EDs[[thresh.list.name]] <- temp$ED
    MCCs[[thresh.list.name]] <- temp$MCC
    Type1s[[thresh.list.name]] <- temp$Type1
    Type2s[[thresh.list.name]] <- temp$Type2
    Overalls[[thresh.list.name]] <- temp$Overall
    NECM_15s[[thresh.list.name]] <- temp$NECM_15
    NECM_20s[[thresh.list.name]] <- temp$NECM_20
    NECM_25s[[thresh.list.name]] <- temp$NECM_25
  }
  
  cutoff.list.names <- c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.05))
  for (cutoff.list.name in cutoff.list.names) {
    temp <- FUN(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, cutoff=cutoff.list.name)
    cutoff.list.name <- paste("P", cutoff.list.name, sep=".")
    Recalls[[cutoff.list.name]] <- temp$RECALL
    Precisions[[cutoff.list.name]] <- temp$PRECISION
    PDs[[cutoff.list.name]] <- temp$PD
    PFs[[cutoff.list.name]] <- temp$PF
    Correctness[[cutoff.list.name]] <- temp$Correctness
    Completeness[[cutoff.list.name]] <- temp$Completeness
    F1s[[cutoff.list.name]] <- temp$F1
    F2s[[cutoff.list.name]] <- temp$F2
    G1s[[cutoff.list.name]] <- temp$G1
    G2s[[cutoff.list.name]] <- temp$G2
    G3s[[cutoff.list.name]] <- temp$G3
    Balances[[cutoff.list.name]] <- temp$Balance
    EDs[[cutoff.list.name]] <- temp$ED
    MCCs[[cutoff.list.name]] <- temp$MCC
    Type1s[[cutoff.list.name]] <- temp$Type1
    Type2s[[cutoff.list.name]] <- temp$Type2
    Overalls[[cutoff.list.name]] <- temp$Overall
    NECM_15s[[cutoff.list.name]] <- temp$NECM_15
    NECM_20s[[cutoff.list.name]] <- temp$NECM_20
    NECM_25s[[cutoff.list.name]] <- temp$NECM_25
  }
  
  #### average performance
  FUN <- ComputeCPAVG
  temp <- FUN(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)
  Recalls[["AVG"]] <- temp$RECALL
  Precisions[["AVG"]] <- temp$PRECISION
  F1s[["AVG"]] <- temp$F1
  
  Recalls <- unlist(Recalls)
  Precisions <- unlist(Precisions)
  PDs <- unlist(PDs)
  PFs <- unlist(PFs)
  Correctness <- unlist(Correctness)
  Completeness<- unlist(Completeness)
  F1s <- unlist(F1s)
  F2s <- unlist(F2s)
  G1s <- unlist(G1s)
  G2s <- unlist(G2s)
  G3s <- unlist(G3s)
  Balances <- unlist(Balances)
  EDs <- unlist(EDs)
  MCCs <- unlist(MCCs)
  Type1s <- unlist(Type1s)
  Type2s <- unlist(Type2s)
  Overalls <- unlist(Overalls)
  NECM_15s <- unlist(NECM_15s)
  NECM_20s <- unlist(NECM_20s)
  NECM_25s <- unlist(NECM_25s)
  
  names(Recalls) <- paste("RECALL", names(Recalls), sep=".")
  names(Precisions) <- paste("PRECISION", names(Precisions), sep=".")
  names(PDs) <- paste("PD", names(PDs), sep=".")
  names(PFs) <- paste("PF", names(PFs), sep=".")
  names(Correctness) <- paste("Correctness", names(Correctness), sep=".")
  names(Completeness) <- paste("Completeness", names(Completeness), sep=".")
  names(G1s) <- paste("G1", names(G1s), sep=".")
  names(G2s) <- paste("G2", names(G2s), sep=".")
  names(G3s) <- paste("G3", names(G3s), sep=".")
  names(F1s) <- paste("F1", names(F1s), sep=".")
  names(F2s) <- paste("F2", names(F2s), sep=".")
  names(Balances) <- paste("Balance", names(Balances), sep=".")
  names(EDs) <- paste("ED", names(EDs), sep=".")
  names(MCCs) <- paste("MCC", names(MCCs), sep=".")
  names(Type1s) <- paste("Type1", names(Type1s), sep=".")
  names(Type2s) <- paste("Type2", names(Type2s), sep=".")
  names(Overalls) <- paste("Overall", names(Overalls), sep=".")
  names(NECM_15s) <- paste("NECM_15", names(NECM_15s), sep=".")
  names(NECM_20s) <- paste("NECM_20", names(NECM_20s), sep=".")
  names(NECM_25s) <- paste("NECM_25", names(NECM_25s), sep=".")
  
  return(list(RECALL=Recalls, PRECISION=Precisions, PD=PDs, PF=PFs, Correctness=Correctness, Completeness=Completeness, F1=F1s, F2=F2s, G1=G1s, G2=G2s, G3=G3s, Balance=Balances, ED=EDs, MCC=MCCs, Type1=Type1s, Type2=Type2s, Overall=Overalls, NECM_15=NECM_15s, NECM_20=NECM_20s, NECM_25=NECM_25s))
}

ComputeCP <- ComputeClassificationPerformance <- function(data, threshold=NULL, cutoff=NULL, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  if (is.null(cutoff)) {
    if (is.null(threshold)) { stop("threshold and cutoff are all NULL") }
    data$PREBIN <- 1 * (data$PRE >= threshold)
  } else {
    data$PREBIN <- 0
    len <- nrow(data) * cutoff
    data$PREBIN[1:len] <- 1
  }
  
  TP <- sum(data$PREBIN*data$REL)
  TN <- sum((1-data$PREBIN)*(1-data$REL))
  FP <- sum(data$PREBIN*(1-data$REL))
  FN <- sum((1-data$PREBIN)*data$REL)
  
  Recall <- 0; if (sum(data$REL)>0) { Recall <- TP/(TP+FN) }
  Precision <- 0; if (sum(data$PREBIN)>0) { Precision <- TP/(TP+FP) }
  
  PD <- 0; if (TP+FN>0) { PD <- TP/(TP+FN) }
  PF <- 0; if (FP+TN>0) { PF <- FP/(FP+TN) }
  
  Correctness <- 0; if (TP+FP>0) { Correctness <- TP/(TP+FP) }
  Completeness<- 0; if (sum(data$NUM)>0) { Completeness <- sum(data$PREBIN*data$REL*data$NUM)/sum(data$NUM) }
  
  F1 <- 0; if (Recall+Precision>0) { F1 <- 2*Recall*Precision/(Recall+1*Precision) }
  F2 <- 0; if (Recall+Precision>0) { F2 <- 5*Recall*Precision/(Recall+4*Precision) }
  G1 <- 2*PD*(1-PF)/(PD+1-PF)
  G2 <- sqrt(Recall*Precision)
  G3 <- sqrt(Recall*(1-PF))
  Balance <- 1 - sqrt(((1-PD)^2+(0-PF)^2)/2)
  ED <- sqrt(0.6*(1-PD)^2+(1-0.6)*PF^2)
  MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  Type1 <- FP/(FP+TN)
  Type2 <- FN/(TP+FN)
  Overall <- (FP+FN)/(TP+TN+FP+FN)
  NECM_15 <- (FP + 15*FN)/(TP+TN+FP+FN)
  NECM_20 <- (FP + 20*FN)/(TP+TN+FP+FN)
  NECM_25 <- (FP + 25*FN)/(TP+TN+FP+FN)
  
  # f1 <- 0; if ((precision+recall)>0) { f1 <- 2*precision*recall/(precision+recall) }
  
  # accuracy <- (sum(data$PREBIN*data$REL) + sum((1-data$PREBIN)*(1-data$REL))) / nrow(data)
  
  # TP <- A <- sum(data$REL*data$PREBIN)
  # FN <- B <- sum(data$REL*(1-data$PREBIN))
  # FP <- C <- sum((1-data$REL)*data$PREBIN)
  # TN <- D <- sum((1-data$REL)*(1-data$PREBIN))
  
  # MCC <- (A*D-B*C)/sqrt((A+B)*(A+C)*(B+D)*(C+D))
  # gmean <- 2*A/(A+B)*D/(C+D)/(A/(A+B) + D/(C+D))
  
  # PD <- TP/(TP+FN)
  # PF <- FP/(FP+TN)
  # Type1 <- FP/(FP+TN)
  # Type2 <- FN/(FN+TP)
  # overall <- (FP+FN)/(TP+TN+FP+FN)
  # NECM_15 <- (FP + 15*FN)/(TP+TN+FP+FN)
  # NECM_20 <- (FP + 20*FN)/(TP+TN+FP+FN)
  # NECM_25 <- (FP + 25*FN)/(TP+TN+FP+FN)
  ## NECM <- FP/(FP+TN)*(FP+TN)/(TP+TN+FP+FN)+CII_CI*FN/(FN+TP)*(FN+TP)/(TP+TN+FP+FN)
  # bal <- 1 - sqrt(((1-PD)^2+PF^2)/2)
  
  return(list(RECALL=Recall, PRECISION=Precision, PD=PD, PF=PF, Correctness=Correctness, Completeness=Completeness, F1=F1, F2=F2, G1=G1, G2=G2, G3=G3, Balance=Balance, ED=ED, MCC=MCC, Type1=Type1, Type2=Type2, Overall=Overall, NECM_15=NECM_15, NECM_20=NECM_20, NECM_25=NECM_25))
  
  #ACCURACY=accuracy, MCC=MCC, GMEAN=gmean, BAL=bal, PF=PF, TYPE1=Type1, TYPE2=Type2, OVERALL=overall, NECM_15=NECM_15, NECM_20=NECM_20, NECM_25=NECM_25))
}

ComputeThreshold <- function(data, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE, allpercentcutoff=TRUE) # lower than is faulty
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  ##############################################################
  #### all possible thresholds.
  ##############################################################
  thresholds <- data$PRE
  
  ##### new code
  len <- nrow(data)
  TPs <- cumsum(data$REL)     ### true positive
  FPs <- cumsum((1-data$REL)) ### false positive
  
  FNs <- TPs[len] - TPs  ### false negative + true positive = all positive
  TNs <- (len - TPs[len]) - FPs ### true negative + false positive = all negative
  
  TPR <- TPs / (TPs + FNs)
  TPR[is.na(TPR)] <- 0
  
  Recall <- TPR
  
  FPR <- FPs / (FPs + TNs)
  FPR[is.na(FPR)] <- 0
  
  Precision <- TPs / (TPs + FPs)
  Precision[is.na(Precision)] <- 0
  
  BPP <- 1 - sqrt( ( (1-TPR)^2 + (0-FPR)^2 ) / 2 )
  BCE <- abs(1-TPR-FPR)
  MFM <- 2 * Recall * Precision / (Recall+Precision)
  MFM[is.na(MFM)] <- 0
  
  df <- data.frame(index=1:len, threshold=thresholds, BPP=BPP, BCE=BCE, MFM=MFM)
  
  if (!allpercentcutoff) {
    pos.equal <- which(data[-len, "PRE"]==data[-1, "PRE"])
    if (length(pos.equal)!=0) {
      df <- df[-pos.equal, ]
    }
  }
  
  V.BPP <- df[which.max(df$BPP), "threshold"]
  V.BCE <- df[which.min(df$BCE), "threshold"]
  V.MFM <- df[which.max(df$MFM), "threshold"]
  
  P.BPP <- df[which.max(df$BPP), "index"] / len
  P.BCE <- df[which.min(df$BCE), "index"] / len
  P.MFM <- df[which.max(df$MFM), "index"] / len
  
  return(list(V.BPP=V.BPP, V.BCE=V.BCE, V.MFM=V.MFM, P.BPP=P.BPP, P.BCE=P.BCE, P.MFM=P.MFM))
}

ComputeCPAVG <- ComputeClassificationPerformanceAVG <- function(data, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE, allpercentcutoff=TRUE)
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  len <- nrow(data)
  
  ### recall
  recalls <- rep(0, len)
  if (sum(data$REL)!=0) {
    recalls <- cumsum(data$REL)/sum(data$REL)
  }
  
  ### precision
  precisions <- cumsum(data$REL) / c(1:nrow(data))
  
  ### f1
  f1s <- 2*precisions*recalls/(precisions+recalls)
  f1s[is.na(f1s)] <- 0
  
  ### accuracy
  accuracys <- (cumsum(data$REL) + c(rev(cumsum(rev(1-data$REL))[-len]), 0)) / nrow(data)
  
  ### filter duplicated items
  if (!allpercentcutoff) {
    pos.equal  <- which(data[-len, "PRE"]==data[-1, "PRE"])
    if (length(pos.equal)!=0) {
      recalls <- recalls[-pos.equal]
      precisions <- precisions[-pos.equal]
      f1s <- f1s[-pos.equal]
      accuracys <- accuracys[-pos.equal]
    }
  }
  
  recall.avg <- mean(recalls)
  precision.avg <- mean(precisions)
  f1.avg <- mean(f1s)
  accuracy.avg <- mean(accuracys)
  
  return(list(RECALL.AVG=recall.avg, PRECISION.AVG=precision.avg, F1.AVG=f1.avg, ACCURACY.AVG=accuracy.avg))
}

ComputeCEs <- function(data, cutoffs=c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1)), effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  res.E <- ComputeCEwithCutoffs(data=data, cutoffs=cutoffs, effort.percent.flag=TRUE,  effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
  CE <- res.E$CE
  AR <- res.E$AR
  NofB <- res.E$NofB
  PofB <- res.E$PofB
  E    <- res.E$E
  S    <- res.E$S
  K    <- res.E$K
  RECALL    <- res.E$RECALL
  PRECISION <- res.E$PRECISION
  ACCURACY  <- res.E$ACCURACY
  F1 <- res.E$F1
  ER <- res.E$ER
  
  return(list(CE=CE, AR=AR, RECALL=RECALL, PRECISION=PRECISION, ACCURACY=ACCURACY, F1=F1, ER=ER, NofB=NofB, PofB=PofB, E=E, S=S, K=K))
}

ComputeCEwithCutoffs <- function(data, cutoffs=c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1)), effort.percent.flag=TRUE, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  data.mdl <- data
  data.opt <- data[order(-data$density, +data$LOC), ]
  data.min <- data[order(+data$density, -data$LOC), ]
  
  areas.opt <- ComputeCEAreawithCutoffs(data=data.opt, cutoffs=cutoffs, effort.percent.flag=effort.percent.flag, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
  areas.mdl <- ComputeCEAreawithCutoffs(data=data.mdl, cutoffs=cutoffs, effort.percent.flag=effort.percent.flag, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
  
  areas.random <- 0.5 * cutoffs * cutoffs
  
  ARs <- areas.mdl$AREA / areas.opt$AREA
  CEs <- (areas.mdl$AREA - areas.random)/(areas.opt$AREA - areas.random)
  
  E <- min(which(cumsum(data.mdl$NUM)/sum(data.mdl$NUM)>0.8))/nrow(data.mdl)
  
  S <- cor(x=data.mdl$PRE, y=data.mdl$NUM, method="spearman")
  K <- cor(x=data.mdl$PRE, y=data.mdl$NUM, method="kendall")
  
  names(E) <- "E"
  names(S) <- "S"
  names(K) <- "K"
  
  ARs <- pmin(1, ARs)
  CEs <- pmin(1, CEs)
  ARs[is.na(ARs)] <- 0
  CEs[is.na(CEs)] <- 0
  
  if (any(is.na(ARs)) | any(is.na(CEs))) { cat("one element of ARs or CEs is NA") }
  
  ### ERs.REL.RAND <- (areas.mdl$RECALL-cutoffs)/(areas.opt$RECALL-cutoffs)
  ERs        <- areas.mdl$ER
  F1s        <- areas.mdl$F1
  NofBs      <- areas.mdl$NofB
  PofBs      <- areas.mdl$PofB
  RECALLs    <- areas.mdl$RECALL
  ACCURACYs  <- areas.mdl$ACCURACY
  PRECISIONs <- areas.mdl$PRECISION
  
  if (effort.percent.flag) {
    names(ARs) <- paste("AR.E", cutoffs, sep=".")
    names(CEs) <- paste("CE.E", cutoffs, sep=".")
  } else {
    names(ARs) <- paste("AR.P", cutoffs, sep=".")
    names(CEs) <- paste("CE.P", cutoffs, sep=".")
  }
  
  return(list(CE=CEs, AR=ARs, RECALL=RECALLs, PRECISION=PRECISIONs, ACCURACY=ACCURACYs, F1=F1s, ER=ERs, NofB=NofBs, PofB=PofBs, E=E, S=S, K=K))
}

ComputeCEAreawithCutoffs <- function(data, cutoffs=c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1)), effort.percent.flag=TRUE, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  len <- nrow(data)
  ### cumulative sums of LOC or NUM
  cumXs <- cumsum(data$LOC) # x: LOC%
  cumYs <- cumsum(data$NUM) # y: Bug%
  
  Xs <- cumXs/cumXs[len]
  Ys <- cumYs/cumYs[len]
  
  poses <- vector(length=length(cutoffs))
  cutoffs.new <- NULL
  if (effort.percent.flag) {
    for (i in seq(cutoffs)) {
      ### find the lowest index which is larger or equal than cutoff
      ### that's to say cutoff is between pos-1 and pos
      cutoff   <- cutoffs[i]
      poses[i] <- min(which(Xs >= cutoff))
      if (cutoff==1) { poses[i] <- len }
    }
    cutoffs.new <- cutoffs ### x-axis is the corresponding cutoff (effort)
  } else {
    poses <- len*cutoffs ### which
    cutoffs.new <- Xs[poses] ### x-axis
  }
  
  Yposes <- vector(length=length(cutoffs))
  if (effort.percent.flag) {
    for (i in seq(cutoffs)) {
      cutoff <- cutoffs[i]
      pos <- poses[i]
      
      if (pos==1) {
        Yposes[i] <- cutoff * (Ys[1]/Xs[1])
      } else {
        Yposes[i] <- (cutoff-Xs[pos-1]) * ((Ys[pos]-Ys[pos-1])/(Xs[pos]-Xs[pos-1])) + Ys[pos-1]
        if (cutoff==1) { Yposes[i] <- 1 }
      }
    }
  } else {
    Yposes <- Ys[poses]
  }
  
  fix_subareas <- vector(length=len)
  fix_subareas[1] <- 0.5 * Ys[1] * Xs[1]
  fix_subareas[2:len] <- 0.5 * (Ys[1:(len-1)] + Ys[2:len]) * abs(Xs[1:(len-1)] - Xs[2:len])
  
  areas      <- vector(length=length(cutoffs))
  RECALLs    <- PofBs <- NofBs <- vector(length=length(cutoffs))
  PRECISIONs <- vector(length=length(cutoffs))
  ACCURACYs  <- vector(length=length(cutoffs))
  for (i in seq(poses)) {
    pos <- poses[i]
    Ypos <- Yposes[i]
    cutoff <- cutoffs.new[i]
    
    subareas <- vector(length=pos)
    if (pos==1) {
      subareas[1] <- 0.5 * cutoff * Ypos
    } else if (pos==2) {
      subareas[1] <- fix_subareas[1]
      subareas[2] <- (Ypos+Ys[1])*(abs(cutoff-Xs[1]))*0.5
    } else {
      subareas[1:(pos-1)] <- fix_subareas[1:(pos-1)]
      subareas[pos]  <- (Ypos+Ys[pos-1])*(abs(cutoff-Xs[pos-1]))*0.5
    }
    
    subareas[is.na(subareas)] <- 0
    
    areas[i]      <- sum(subareas)
    NofBs[i]      <- cumYs[pos]
    PofBs[i]      <- cumYs[pos] / cumYs[len]
    RECALLs[i]    <- cumYs[pos] / cumYs[len]
    PRECISIONs[i] <- cumYs[pos] / pos
    ACCURACYs[i]  <- (cumYs[pos] + ((len-pos)-(cumYs[len]-cumYs[pos]))) / len
  }
  
  F1s <- 2*RECALLs*PRECISIONs/(RECALLs+PRECISIONs)
  F1s[is.na(F1s)] <- 0
  if (any(is.na(areas)) | any(is.na(RECALLs)) | any(is.na(PRECISIONs)) | any(is.na(ACCURACYs))) {
    if (sum(data$NUM) == 0) {
      areas[is.na(areas)] <- 0
      RECALLs[is.na(RECALLs)] <- 0
      PRECISIONs[is.na(PRECISIONs)] <- 0
      ACCURACYs[is.na(ACCURACYs)] <- 0
    } else {
      stop("one element of areas, recalls, precisions, or accuracys is NA")
    }
  }
  
  if (effort.percent.flag) {
    names(F1s)        <- paste("F1.E", cutoffs, sep=".")
    names(RECALLs)    <- paste("RECALL.E", cutoffs, sep=".")
    names(PofBs)      <- paste("PofB.E", cutoffs, sep=".")
    names(NofBs)      <- paste("NofB.E", cutoffs, sep=".")
    names(ACCURACYs)  <- paste("ACCURACY.E", cutoffs, sep=".")
    names(PRECISIONs) <- paste("PRECISION.E", cutoffs, sep=".")
  } else {
    names(F1s)        <- paste("F1.P", cutoffs, sep=".")
    names(RECALLs)    <- paste("RECALL.P", cutoffs, sep=".")
    names(PofBs)      <- paste("PofB.P", cutoffs, sep=".")
    names(NofBs)      <- paste("NofB.P", cutoffs, sep=".")
    names(ACCURACYs)  <- paste("ACCURACY.P", cutoffs, sep=".")
    names(PRECISIONs) <- paste("PRECISION.P", cutoffs, sep=".")
  }
  
  ERs <- NULL
  for (i in seq(cutoffs)) {
    pos <- poses[i]
    data$PREBIN <- 0
    data$PREBIN[1:pos] <- 1
    
    Effort.M <- cumXs[pos] / cumXs[len] ### the effort required to inspect the top pos modules of model m
    Effort.R <- RECALLs[i] ### the effort required by random model to achieve same recall
    
    ER.TNE <- sum((1-data$PREBIN)*(1-data$REL)*data$LOC)/sum(data$LOC) ### TP/FP肯定要审查, FN部分在运行期间引发失效后也会被审查. 与全部审查模型相比, 模型m节省的工作量约简为 ER' = [Effort(all) - Effort(m)] / Effort(all) = SLOC(TN) / SLOC(TP+FP+FN+TN)
    ER.ABS.REL <- Effort.R - Effort.M ### 绝对的ER 相对于整个系统代码审查的百分比
    ER.ABS.LOC <- ER.ABS.REL*sum(data$LOC) ### 绝对的ER 相对于整个系统代码审查的百分比
    ER.REL.RAND <- (Effort.R - Effort.M) / Effort.R ### 相对于随机模型的ER
    ER.REL.SELF <- (Effort.R - Effort.M) / Effort.M ### 相对于自身模型的ER
    
    if (Effort.R==0) { ER.REL.RAND <- NA }
    if (Effort.M==0) { ER.REL.SELF <- NA }
    
    aER <- c(ER.TNE, ER.ABS.REL, ER.ABS.LOC, ER.REL.RAND, ER.REL.SELF)
    if (effort.percent.flag) {
      names(aER) <- paste(c("ER.TNE", "ER.ABS.REL", "ER.ABS.LOC", "ER.REL.RAND", "ER.REL.SELF"), "E", cutoffs[i], sep=".")
    } else {
      names(aER) <- paste(c("ER.TNE", "ER.ABS.REL", "ER.ABS.LOC", "ER.REL.RAND", "ER.REL.SELF"), "P", cutoffs[i], sep=".")
    }
    
    ERs <- c(ERs, aER)
  }
  
  if (sum(data$NUM) == 0) {
    ERs[is.na(ERs)] <- 0
  }
  
  return(list(AREA=areas, NofB=NofBs, PofB=PofBs, RECALL=RECALLs, PRECISION=PRECISIONs, ACCURACY=ACCURACYs, F1=F1s, ER=ERs))
}

ComputeAUCEC <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{	
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  data.mdl <- data
  data.opt <- data[order(-data$density, +data$LOC), ]
  data.wst <- data[order(+data$density, -data$LOC), ]
  
  area.mdl <- ComputeArea(data=data.mdl, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
  area.opt <- ComputeArea(data=data.opt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
  area.wst <- ComputeArea(data=data.wst, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
  
  AUCEC <- area.mdl/area.opt
  names(AUCEC) <- "AUCEC"
  return(AUCEC)
}

ComputeArea <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
{
  if (!sorted) {
    data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
    sorted <- TRUE
  }
  
  len <- nrow(data)
  ### cumulative sums of LOC or NUM
  cumXs <- cumsum(data$LOC) # x: LOC%
  cumYs <- cumsum(data$NUM) # y: Bug%
  
  Xs <- cumXs/1000 ### Xs <- cumXs/cumXs[len]
  Ys <- cumYs/cumYs[len]
  
  fix_subareas <- vector(length=len)
  fix_subareas[1] <- 0.5 * Ys[1] * Xs[1]
  fix_subareas[2:len] <- 0.5 * (Ys[1:(len-1)] + Ys[2:len]) * abs(Xs[1:(len-1)] - Xs[2:len])
  
  area <- sum(fix_subareas)
  
  return(area)
}








require(RWeka)

#source("SoftEngDataAna.R")

root <- "(final)data"
paths <- list.files(path=root, full.names=FALSE)

paths <- c("")
for (method in c("UP", "DOWN")) {
	for (path in paths) {
		cat(path, "\n")
		# data_path <- file.path(path, "data")
		# files <- list.files(path=file.path(root, data_path), full.names=TRUE)

		files <- list.files(path=file.path(root, path), full.names=TRUE)
		files <- files[!R.utils::isDirectory(files)]

		all.Recall <- all.Precision <- all.PD <- all.PF <- all.Correctness <- all.Completeness <- all.F1 <- all.F2 <- all.G1 <- all.G2 <- all.G3 <- all.Balance <- all.ED <- all.MCC <- all.AUC <- all.Type1 <- all.Type2 <- all.Overall <- all.NECM <- NULL ### classification
		
		all.AUCEC <- all.POPT <- all.NofB <- all.PofB <- all.E <- all.S <- all.K <- NULL ### ranking
		
		df.stat <- NULL
		for (file in files) {
			cat(file, "\n")
			
#			if (grepl("Xia TSE 2016", file)) {
#				data <- RWeka::read.arff(file)
#				data$LOC <- data$loc
#				data$BUG <- 1*(data$bug==1)
#				data$NUM <- data$BUG
#			} else {
#				data <- read.csv(file, sep=",", stringsAsFactors=FALSE)
#				if (grepl("Kamei EMSEJ 2016", file)) {
#					data$BUG <- 1*(data$bug>=1)
#					data$NUM <- as.numeric(data$bug)
#					data$LOC <- (data$la+data$ld)*(data$lt*data$nf)
#				} else {
#					if ("bugs" %in% colnames(data)) {
#						data$BUG <- 1*(data$bugs>=1)
#						data$NUM <- as.numeric(data$bugs)
#					} else if ("bug" %in% colnames(data)) {
#						if (file=="(ok)comparison/Liu TSE 2010/data/kc2.csv") {
#							data$BUG <- 1*(data$bug=="yes")
#							data$NUM <- data$BUG
#						} else {
#							data$BUG <- 1*(data$bug>=1)
#							data$NUM <- as.numeric(data$bug)
#						}
#					} else if ("ERROR_COUNT" %in% colnames(data)) {
#						data$BUG <- 1*(data$ERROR_COUNT>=1)
#						data$NUM <- as.numeric(data$ERROR_COUNT)
#					} else if ("defects" %in% colnames(data)) {
#						data$BUG <- 1*(data$defects=="true")
#						data$NUM <- data$BUG
#					} else if ("isDefective" %in% colnames(data)) {
#						data$BUG <- 1*(data$isDefective==TRUE)
#						data$NUM <- data$BUG
#					} else if ("Defective" %in% colnames(data)) {
#						data$BUG <- 1*(data$Defective=="Y")
#						data$NUM <- data$BUG
#					} else if ("problems" %in% colnames(data)) {
#						data$BUG <- 1*(data$problems=="yes")
#						data$NUM <- data$BUG
#					} else if ("label" %in% colnames(data)) {
#						data$BUG <- 1*(data$label=="Y")
#						data$NUM <- data$BUG
#					} else if ("hasb" %in% colnames(data)) {
#						data$BUG <- 1*(data$hasb=="yes")
#						data$NUM <- data$BUG
#					} else if ("c" %in% colnames(data)) {
#						data$BUG <- 1*(data$c=="TRUE")
#						data$NUM <- data$BUG
#					} else {
#						browser()
#					}
#
#					if ("CountLineCode" %in% colnames(data)) {
#						data$LOC <- data$CountLineCode
#					} else if ("numberOfLinesOfCode" %in% colnames(data)) {
#						data$LOC <- data$numberOfLinesOfCode
#					} else if ("loc" %in% colnames(data)) {
#						data$LOC <- data$loc
#					} else if ("total_loc" %in% colnames(data)) {
#						data$LOC <- data$total_loc
#					} else if ("LOC_TOTAL" %in% colnames(data)) {
#						data$LOC <- data$LOC_TOTAL
#					} else if ("lOCode" %in% colnames(data)) {
#						data$LOC <- data$lOCode
#					} else if ("TLOC" %in% colnames(data)) {
#						data$LOC <- data$TLOC
#					} else if ("size" %in% colnames(data)) {
#						data$LOC <- data$size				
#					} else if ("CountLine" %in% colnames(data)) {
#						data$LOC <- data$CountLine
#					} else if ("loc_total" %in% colnames(data)) {
#						data$LOC <- data$loc_total
#					} else {
#						browser()
#					}
#				}
#			}

      #?????ļ??ĵ?һ????SLOC??????һ????bug??Ϣ
      data <- read.csv(file, sep=",", stringsAsFactors=FALSE)
      
      data$NUM <- as.numeric(data[,ncol(data)])
      data$LOC <- as.numeric(data[,1])
      data$BUG <- 1*(as.numeric(data$NUM>=1))

			data[is.na(data)] <- 0
			

			
			df.stat <- rbind(df.stat, c(nrow(data), sum(data$BUG), sum(data$BUG)/nrow(data), sum(data$NUM)))

			cname <- "LOC"
			
			# valid.dt.up <- data.frame(NUM=data$BUG, REL=data$BUG, LOC=data$LOC, PRE=1000/(data[, cname]+1))
			# valid.dt.dw <- data.frame(NUM=data$BUG, REL=data$BUG, LOC=data$LOC, PRE=data[, cname])

			valid.dt.up <- valid.dt.up.num <- data.frame(NUM=data$NUM, REL=data$BUG, LOC=data$LOC, PRE=1000/(data[, cname]+1))
			valid.dt.dw <- valid.dt.dw.num <- data.frame(NUM=data$NUM, REL=data$BUG, LOC=data$LOC, PRE=data[, cname])
			
			valid.dt <- NULL
			if (method=="UP") {
				valid.dt <- valid.dt.up
			} else {
				valid.dt <- valid.dt.dw
			}
			
			sorted           <- FALSE
			worstcase        <- TRUE
			bestcase         <- FALSE
			LOCUP            <- FALSE
			allpercentcutoff <- FALSE


			AUC <- ComputeAUC(sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff, data=valid.dt)
			cp <- ComputeCPs(sorted=sorted,  worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff, train.dt=valid.dt, valid.dt=valid.dt)
			ce <- ComputeCEs(sorted=sorted, data=valid.dt, cutoffs=seq(0.1, 1, 0.1), worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
			aucec <- ComputeAUCEC(sorted=sorted, data=valid.dt, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)

			
			aRecall <- cp[["RECALL"]][paste("RECALL.P.", seq(0.1, 1, 0.1), sep="")]
			aPrecision <- cp[["PRECISION"]][paste("PRECISION.P.", seq(0.1, 1, 0.1), sep="")]

			aPD <- cp[["PD"]][paste("PD.P.", seq(0.1, 1, 0.1), sep="")]			
			aPF <- cp[["PF"]][paste("PF.P.", seq(0.1, 1, 0.1), sep="")]

			aCorrectness <- cp[["Correctness"]][paste("Correctness.P.", seq(0.1, 1, 0.1), sep="")]
			aCompleteness <- cp[["Completeness"]][paste("Completeness.P.", seq(0.1, 1, 0.1), sep="")]
			
			aF1 <- cp[["F1"]][paste("F1.P.", seq(0.1, 1, 0.1), sep="")]
			aF2 <- cp[["F2"]][paste("F2.P.", seq(0.1, 1, 0.1), sep="")]
			
			aG1 <- cp[["G1"]][paste("G1.P.", seq(0.1, 1, 0.1), sep="")]
			aG2 <- cp[["G2"]][paste("G2.P.", seq(0.1, 1, 0.1), sep="")]
			aG3 <- cp[["G3"]][paste("G3.P.", seq(0.1, 1, 0.1), sep="")]

			aBalance <- cp[["Balance"]][paste("Balance.P.", seq(0.1, 1, 0.1), sep="")]
			aED <- cp[["ED"]][paste("ED.P.", seq(0.1, 1, 0.1), sep="")]
			aMCC <- cp[["MCC"]][paste("MCC.P.", seq(0.1, 1, 0.1), sep="")]
			aAUC <- AUC
			
			aType1 <- cp[["Type1"]][paste("Type1.P.", seq(0.1, 1, 0.1), sep="")]
			aType2 <- cp[["Type2"]][paste("Type2.P.", seq(0.1, 1, 0.1), sep="")]
			aOverall <- cp[["Overall"]][paste("Overall.P.", seq(0.1, 1, 0.1), sep="")]
			
			aNECM_15 <- cp[["NECM_15"]][paste("NECM_15.P.", seq(0.1, 1, 0.1), sep="")]
			aNECM_20 <- cp[["NECM_20"]][paste("NECM_20.P.", seq(0.1, 1, 0.1), sep="")]
			aNECM_25 <- cp[["NECM_25"]][paste("NECM_25.P.", seq(0.1, 1, 0.1), sep="")]
			
			
			###Ranking
			aAUCEC <- aucec
			aPOPT <- 0.5*ce$CE[paste("CE.E.", seq(0.1, 1, 0.1), sep="")]+0.5
			aNofB <- ce$NofB[paste("NofB.E.", seq(0.1, 1, 0.1), sep="")]
			aPofB <- ce$PofB[paste("PofB.E.", seq(0.1, 1, 0.1), sep="")]
			aE <- ce$E
			aS <- ce$S
			aK <- ce$K

			### classification
			all.Recall    <- rbind(all.Recall, aRecall)
			all.Precision <- rbind(all.Precision, aPrecision)
			all.PD <- rbind(all.PD, aPD)
			all.PF <- rbind(all.PF, aPF)
			all.Correctness  <- rbind(all.Correctness, aCorrectness)
			all.Completeness <- rbind(all.Completeness, aCompleteness)
			all.F1 <- rbind(all.F1, aF1)
			all.F2 <- rbind(all.F2, aF2)
			all.G1 <- rbind(all.G1, aG1)
			all.G2 <- rbind(all.G2, aG2)
			all.G3 <- rbind(all.G3, aG3)
			all.Balance <- rbind(all.Balance, aBalance)
			all.ED  <- rbind(all.ED, aED)
			all.MCC <- rbind(all.MCC, aMCC)
			all.AUC <- rbind(all.AUC, aAUC)
			all.Type1 <- rbind(all.Type1, aType1)
			all.Type2 <- rbind(all.Type2, aType2)
			all.Overall <- rbind(all.Overall, aOverall)
			all.NECM <- rbind(all.NECM, c(aNECM_15, aNECM_20, aNECM_25))

			### ranking
			all.AUCEC <- rbind(all.AUCEC, aAUCEC)
			all.POPT  <- rbind(all.POPT, aPOPT)
			all.NofB <- rbind(all.NofB, aNofB)
			all.PofB <- rbind(all.PofB, aPofB)
			all.E <- rbind(all.E, aE)
			all.S <- rbind(all.S, aS)
			all.K <- rbind(all.K, aK)
		}

		path_res <- file.path(root, path, "result")
		if (!dir.exists(path_res)) { dir.create(path_res) }
		
		colnames(df.stat) <- c("N", "#BUGGY", "#buggyrate", "#bugs")
		rownames(df.stat) <- basename(files)
		write.csv(df.stat, file=file.path(path_res, "summary-all.csv"))
		
		rownames(all.Recall) <- rownames(all.Precision) <- rownames(all.PD) <- rownames(all.PF) <- rownames(all.Correctness) <- rownames(all.Completeness) <- rownames(all.F1) <- rownames(all.F2) <- rownames(all.G1) <- rownames(all.G2) <- rownames(all.G3) <- rownames(all.Balance) <- rownames(all.ED) <- rownames(all.MCC) <- rownames(all.AUC) <- rownames(all.Type1) <- rownames(all.Type2) <- rownames(all.Overall) <- rownames(all.NECM) <- basename(files)
		rownames(all.AUCEC) <- rownames(all.POPT) <- rownames(all.NofB) <- rownames(all.PofB) <- rownames(all.E) <- rownames(all.S) <- rownames(all.K) <- basename(files)


		subOtherMean <- function(data) {
			data.revise <- data
			colnames(data.revise) <- paste(colnames(data), ".revise", sep="")
			for (i in 1:nrow(data)) {
				for (j in 1:ncol(data)) {
					data.revise[i, j] <- mean(data[-1*i, j])
				}
			}
			return(cbind(data, data.revise))
		}

		# all.Recall <- subOtherMean(all.Recall)
		# all.Precision <- subOtherMean(all.Precision)
		# all.PD <- subOtherMean(all.PD)
		# all.PF <- subOtherMean(all.PF)
		# all.Correctness <- subOtherMean(all.Correctness)
		# all.Completeness <- subOtherMean(all.Completeness)
		# all.F1 <- subOtherMean(all.F1)
		# all.F2 <- subOtherMean(all.F2)
		# all.G1 <- subOtherMean(all.G1)
		# all.G2 <- subOtherMean(all.G2)
		# all.G3 <- subOtherMean(all.G3)
		# all.Balance <- subOtherMean(all.Balance)
		# all.ED <- subOtherMean(all.ED)
		# all.MCC <- subOtherMean(all.MCC)
		# all.AUC <- subOtherMean(all.AUC)
		# all.Type1 <- subOtherMean(all.Type1)
		# all.Type2 <- subOtherMean(all.Type2)
		# all.Overall <- subOtherMean(all.Overall)
		# all.NECM  <- subOtherMean(all.NECM)
		
		# ### Ranking
		# all.AUCEC <- subOtherMean(all.AUCEC)
		# all.POPT <- subOtherMean(all.POPT)
		# all.NofB <- subOtherMean(all.NofB)
		# all.PofB <- subOtherMean(all.PofB)
		# all.E <- subOtherMean(all.E)
		# all.S <- subOtherMean(all.S)
		# all.K <- subOtherMean(all.K)
		
		### Classification
		write.csv(round(all.Recall, 3), file=file.path(path_res, paste("all-Recall-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.Precision, 3), file=file.path(path_res, paste("all-Precision-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.PD, 3), file=file.path(path_res, paste("all-PD-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.PF, 3), file=file.path(path_res, paste("all-PF-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.Correctness, 3), file=file.path(path_res, paste("all-Correctness-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.Completeness, 3), file=file.path(path_res, paste("all-Completeness-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.F1, 3), file=file.path(path_res, paste("all-F1-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.F2, 3), file=file.path(path_res, paste("all-F2-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.G1, 3), file=file.path(path_res, paste("all-G1-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.G2, 3), file=file.path(path_res, paste("all-G2-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.G3, 3), file=file.path(path_res, paste("all-G3-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.Balance, 3), file=file.path(path_res, paste("all-Balance-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.ED, 3), file=file.path(path_res, paste("all-ED-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.MCC, 3), file=file.path(path_res, paste("all-MCC-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.AUC, 3), file=file.path(path_res, paste("all-AUC-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.Type1, 3), file=file.path(path_res, paste("all-Type1-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.Type2, 3), file=file.path(path_res, paste("all-Type2-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.Overall, 3), file=file.path(path_res, paste("all-Overall-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.NECM, 3), file=file.path(path_res, paste("all-NECM-", method, ".csv", sep="")), row.names=TRUE)
		
		### Ranking
		write.csv(round(all.AUCEC, 3), file=file.path(path_res, paste("all-AUCEC-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.POPT, 3), file=file.path(path_res, paste("all-POPT-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.NofB, 3), file=file.path(path_res, paste("all-NofB-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.PofB, 3), file=file.path(path_res, paste("all-PofB-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.E, 3), file=file.path(path_res, paste("all-E-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.S, 3), file=file.path(path_res, paste("all-S-", method, ".csv", sep="")), row.names=TRUE)
		write.csv(round(all.K, 3), file=file.path(path_res, paste("all-K-", method, ".csv", sep="")), row.names=TRUE)
		
		

	}
}
