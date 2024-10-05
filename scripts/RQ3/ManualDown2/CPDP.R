
require(RWeka)

source("SoftEngDataAna.R")

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

      #输入文件的第一列是SLOC，最后一列是bug信息
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
