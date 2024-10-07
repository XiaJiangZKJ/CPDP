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


ComputeCPs <- function(train.dt, valid.dt, thresholds=NULL, effortaware=FALSE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE, allpercentcutoff=TRUE)
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
