### 将R界面设置成英文
Sys.setenv(LANGUAGE="en")
##############################
### install packages
##############################
packages <- c("beanplot", "car", "class", "doSNOW", "e1071", "effsize", "foreach", "FSelector", "hmeasure", "igraph", "kernlab", "kknn", "lmtest", "markovchain", "moments", "monmlp", "neuralnet", "nnet", "perturb", "plyr", "PMCMR", "psych", "randomForest", "rms", "rpart", "RSNNS", "RWeka", "tree", "usdm")
packages <- packages[!(packages %in% installed.packages())]
if (length(packages) != 0) { install.packages(packages, repos="http://cran.rstudio.com/") }


### 绘制LOC%-defect% Alberg 图
### plot Alberg graph
AlbergCurve <- function(height=1.2)
{
    Xs   <- seq(0, 1, 0.0001)
    Ys.m <- sqrt(1-(1-Xs)^2)

    per <- c(seq(0.8/Ys.m[2000]*1.08, 0.8/Ys.m[2000], length=length(Xs)*0.2), seq(0.8/Ys.m[2000], 1.15, length=length(Xs)*0.3), seq(1.15, 1.0, length=length(Xs)*0.5))
    per <- per[1:length(Xs)]
    Ys.o <- Ys.m*per
    Ys.o[Ys.o>1] <- 1

    per <- c(seq(0.85, 0.95, length=length(Xs)*0.5), seq(0.95, 1.0, length=length(Xs)*0.5))
    per <- per[1:length(Xs)]
    Ys.m <- Ys.m*per

    Ys.r <- Xs

    Xs <- Xs*100
    Ys.m <- Ys.m*100
    Ys.o <- Ys.o*100
    Ys.w <- 100-Ys.o
    Ys.r <- Ys.r*100

    golden.ratio <- 2/(sqrt(5) + 1)
    width <- height/golden.ratio
    png(file="Alberg.png", width=width, height=height, units="in", res=600, pointsize=10)
    par(mai=c(0, 0, 0, 0), omi=c(0.35, 0.35, 0.1, 0.1), mex=0.4, cex=1.0)
    plot(Ys.o~Xs, ylab="", xlab="", xlim=c(0, 100), ylim=c(0, 100), axes=FALSE, type="l", lwd=2)
    lines(Xs, Ys.r, col="gray", lty='solid',  lwd=2)
    lines(Xs, Ys.m, col="blue", lty='dashed', lwd=2)
    # lines(rev(Xs), Ys.w, col="gray", lty='solid',  lwd=2)

    abline(v=20, lty='dashed')

    box(bty='l'); ### 连接x轴与y轴
    axis(1, mgp=c(0.5, 0.5, 0), tck=-0.02, cex.axis=0.6)
    axis(2, mgp=c(0.5, 0.5, 0), tck=-0.02, cex.axis=0.6);
    mtext("% faults",  side=2.5, line=2.5, adj=0.5, cex=0.6)
    # mtext("% code churn", side=1, line=2.5, adj=0.5, cex=0.6)
    mtext("% SLOC", side=1, line=2, adj=0.5, cex=0.8)

    ### expression(italic("prediction model"), "Random model", "Optimal model") ### expression(italic("prediction model"), italic("Random model"), italic("Optimal model"))
    # legend(x="bottomright", legend=c("prediction model", "Optimal model", "Worst model"), lty=c('dashed', 'solid', 'solid'), lwd=c(1.5, 1.5, 1.5, 1.5), col=c("blue", "black", "gray"), cex=0.6, merge=FALSE, bty = "n")
    legend(x="bottomright", legend=c("prediction model", "Optimal model", "random model"), lty=c('dashed', 'solid', 'solid'), lwd=c(1.5, 1.5, 1.5, 1.5), col=c("blue", "black", "gray"), cex=0.6, merge=FALSE, bty = "n")

    dev.off()
}
### AlbergCurve()

### description of the data set
DescriptiveAna <- function(data, xnames, file="desAna.csv")
{
	require(moments)
	desc <- matrix(0, nrow=ncol(data), ncol=9)
	rownames(desc) <- colnames(data)
	colnames(desc) <- c("Min", "25%", "Median", "Mean", "75%", "Max", "Std.Dev.", "Skewness", "Kurtosis")

	for (i in 1:ncol(data)) {
		si <- as.matrix(summary(data[,i]))
 		desc[i, 1:6] <- t(si)
		desc[i, 7] <- sd(data[,i])
		desc[i, 8] <- skewness(data[,i])
		desc[i, 9] <- kurtosis(data[,i])
	}

 	desc <- round(desc, 3)
	write.table(c(""), file=file, row.names=FALSE, col.names=FALSE, eol=",", append=FALSE)
	write.table(desc, file=file, sep=",", row.names=TRUE, col.names=TRUE, append=TRUE)
}

### correlation analysis (pearson or spearman)
CorrelationAna <- function(data, xnames, test="pearson", file="CorAna.csv")
{
	cor.e <- matrix(0, nrow=length(xnames), ncol=length(xnames))
	colnames(cor.e) <- xnames
	rownames(cor.e) <- xnames
	cor.p <- cor.e

	for (i in xnames) {
		for (j in xnames) {
			ct <- cor.test(data[, i], data[, j], method=test)
			cor.e[i, j] <- ct$estimate
			cor.p[i, j] <- ct$p.value
		}
	}
	cor.e <- round(cor.e, 3)
	cor.p <- round(cor.p, 3)
	write.table(c(""), file=file, row.names=FALSE, col.names=FALSE, eol=",", append=FALSE)
	write.table(cor.e, file=file, sep=",", row.names=TRUE, col.names=TRUE, append=TRUE)
	write.table(c(""), file=file, row.names=FALSE, col.names=FALSE, eol=",", append=TRUE)
	write.table(cor.p, file=file, sep=",", row.names=TRUE, col.names=TRUE, append=TRUE)
}

### PCA analysis
PrincipalComponentAna <- function(data, xnames, sort=FALSE, file="pca.csv", rotate="varimax")
{
	if (is.na(xnames)) { xnames <- colnames(data) }
	pca <- txnames <- NULL
	for (xname in xnames) {
		if (length(unique(data[, xname]))==1) { next } else { txnames <- c(txnames, xname) }
	}; xnames <- txnames

	for (incr in c(0:length(xnames))) {
		err_catch <- try(pca <- psych::principal(data[, xnames], nfactors=length(xnames)-incr, rotate=rotate, scores=TRUE, residuals=TRUE))
		if (class(err_catch)!="try-error") { break }
	}

	#### format similar with SPSS output
	mloadings <- unclass(pca$loadings) ### unclass loadings to matrix ### mloadings <- loadings(pca)[, 1:nrow(loadings(pca))] ### dloadings <- as.data.frame(loadings(pca)[, 1:dim(pca$loadings)[2]])
	eig <- colSums(mloadings*mloadings) ### egien vectors
	variance <- colSums(mloadings*mloadings)/nrow(mloadings) ### proportion variance ### variance <- eig/dim(pca$loadings)[1]

	df <- data.frame(SSloadings=eig, ProportionVar=variance, Variables=0) ### df <- data.frame(SSloadings=eig, ProportionVar=variance, CumulativeVar=cum_vars, Variables=0)
	rotated.pca <- rbind(t(df), mloadings)
	rotated.pca <- round(rotated.pca, 3)

	clusters <- apply(abs(mloadings), MARGIN=1, which.max) ### for each row, which one is max
	for (index in sort(unique(clusters))) {
		groupnames <- names(clusters)[clusters==index]
		df[index, "Variables"] <- paste(groupnames, collapse="+")
	}
	df <- df[df$Variables!=0, ]; df$CumulativeVar <- cumsum(df$ProportionVar)

	df <- df[, c("SSloadings", "ProportionVar", "CumulativeVar", "Variables")]
	df[, c("SSloadings", "ProportionVar", "CumulativeVar")] <- round(df[, c("SSloadings", "ProportionVar", "CumulativeVar")], 3)
	df$CumulativeVar <- cumsum(df$ProportionVar)
	write.table(df, file=file, sep=",", row.names=TRUE, col.names=TRUE, append=FALSE)
	write.table(rotated.pca, file=file, sep=",", row.names=TRUE, col.names=TRUE, append=TRUE)
	sink(file, append=TRUE); print(pca$loadings, cutoff=0.0); sink()
}

### normal R^2(linear model)
GlmR2.LinearRegression <- function(model)
{
	yrel <- model$y ### real y
	ybar <- mean(yrel) ### mean of real y
	ypre <- predict(model, type = "response") ### predicted value of y

	SST <- sum((yrel - ybar)^2)
	SSE <- sum((yrel - ypre)^2)
	R2  <- 1 - SSE/SST
	return(R2)
}

### Adjusted R^2 for linear regression
GlmR2Adjusted.LinearRegression <- function(model)
{
	n <- nrow(model$data)
	p <- length(attr(model$terms, "variables"))
	R2 <- 0
	b = model$y
	ybar <- mean(b)

	pre = predict( model, type = "response")

	SST <- sum((b - ybar)^2)
	SSE <- sum((b - pre)^2)

	R2adj <- 1 - (SSE/(n-p-1)) / (SST/(n-1))
	return(R2adj)
}

### Logistic regression model Pseudo R^2
LogisticGlmPseudoR2 <- function(model)
{
	null.model <- update(model, ~ 1, data=model$data)
	N <- length(model$y)
	lrt <- lmtest::lrtest(model, null.model)
	LR <- lrt$LogLik[1]; L0 <- lrt$LogLik[2]
	R2 <- (1-exp(-LR/N))/(1-exp(2*L0/N))
	return(R2)
}

### univariate logistic regression analysis
UnivariateLogisticAna <- function(model="lr", family=binomial, data, yname, xnames, method="normal", effort.name="LOC", number.name=yname, sampling=FALSE, weight.name="weight", remove.flag=FALSE, weight.flag=FALSE, cooksd.flag=FALSE, cookthreshold=1, file="ulr.csv", trace=FALSE)
{
	ulr.results <- data.frame(matrix(NA, nrow=length(xnames), ncol=10))
	rownames(ulr.results) <- xnames
	colnames(ulr.results) <- c("Constant", "Coefficient", "p-value", "delta_OR", "Std.err", "PseudoR2", "AdjR2", "N.EXCEPTION", "AUC", "p-value-adjust")

	if (remove.flag) {
		data <- RemoveConfoundingEffect(data, ynames=xnames, confeff.name=effort.name)
	}

	for (i in 1:length(xnames)) {
		xname <- xnames[i]
		if (trace) { cat("Independent name:", xname, "\n") } ### if (mode(data[, xname])=="complex") { next }

		tdata <- data
		fm <- paste(yname, xname, sep="~") ### ("%s ~ %s", yname, xname)
		ul <- NULL
		err_catch1 <- try(ul <- glm(as.formula(fm), data=tdata, family=binomial), silent=TRUE)
		if (class(err_catch1)=="try-error") {
			next
		} else {
			if (cooksd.flag) {
				cookd <- cooks.distance(ul)
				cookd[is.na(cookd)] <- cookthreshold + 1
				cookouts <- sum(cookd>=cookthreshold)
				# if (cookouts!=0) { break }
				# tdata <- tdata[-which(cookd>=cookthreshold), ]
				# if (trace) { cat("number of outliers: ", cookouts, "\n") }
			}
		}

		ul.sum <- summary(ul)
		pre <- predict(ul, data, type="response")
		if (tolower(method)=="divide") {
			valid.dt <- data.frame(REL=data[, yname], LOC=data[, effort.name], NUM=data[, number.name], PRE=pre/data[, effort.name])
		} else {
			valid.dt <- data.frame(REL=data[, yname], LOC=data[, effort.name], NUM=data[, number.name], PRE=pre)
		}

		N.EXCEP <- nrow(data) - nrow(ul$data)
		auc <- ComputeAUC(valid.dt)
		if (is.na(ul$coefficients[xname])) {
			ulr.results[i, "Constant"] <- ul$coefficients["(Intercept)"]
			ulr.results[i, "Coefficient"] <- NA
			ulr.results[i, "Std.err"] <- NA
			ulr.results[i, "p-value"] <- NA
			ulr.results[i, "R2"] <- NA
			ulr.results[i, "AdjR2"] <- NA
			ulr.results[i, "N.EXCEPTION"] <- N.EXCEP
			ulr.results[i, "AUC"] <- NA
			ulr.results[i, "delta_OR"] <- NA
		} else {
			ulr.results[i, "Constant"]    <- ul.sum$coefficients["(Intercept)", "Estimate"]
			ulr.results[i, "Coefficient"] <- ul.sum$coefficients[xname, "Estimate"]
			ulr.results[i, "Std.err"] <- ul.sum$coefficients[xname, "Std. Error"]
			ulr.results[i, "p-value"] <- ul.sum$coefficients[xname, "Pr(>|z|)"]
			err_catch2 <- try(lr.m <- rms::lrm(formula=as.formula(fm), data=ul$data), silent=TRUE)
			if (class(err_catch2)=="try-error") {
				ulr.results[i, "PseudoR2"] <- LogisticGlmPseudoR2(ul) ###
			} else {
				ulr.results[i, "PseudoR2"] <- lr.m$stats["R2"] ### LogisticGlmPseudoR2(ul) ###
			}

			ulr.results[i, "AdjR2"] <- GlmR2Adjusted.LinearRegression(ul)
			ulr.results[i, "N.EXCEPTION"] <- N.EXCEP
			ulr.results[i, "AUC"] <- auc

			coef <- ul$coefficients[xname] # coefficient
			asd  <- sd(tdata[, xname])
			ulr.results[i, "delta_OR"] <- exp(coef * asd) # delta_OR
		}
	}

	ulr.results[, "p-adj"] <- p.adjust(p=ulr.results[, "p-value"], method="BH")
	return(ulr.results)
}

UnivariateLogisticAna.CrossValidation <- function(data, yname, xnames, model="lr", family="binomial", effort.name="LOC", number.name="BUGS", totalFolds=10, totalRuns=10,  sampling=FALSE, remove.flag=FALSE, weight.flag=FALSE, cooksd.flag=FALSE, cookthreshold=1, weight.name="weight", file)
{
	auc.results <- matrix(0, nrow=totalFolds*totalRuns, ncol=length(xnames))
	colnames(auc.results) <- xnames
	auc.results <- as.data.frame(auc.results)
	erbpp.results <- erbce.results <- auc.results

	for (i in 1:length(xnames)) {
		cv.result <- CrossValidation(data=data, yname=yname, xnames=xnames[i], model=model, family=family, effort.name=effort.name, number.name=number.name, totalFolds=totalFolds, totalRuns=totalRuns, sampling=sampling, remove.flag=remove.flag, weight.flag=weight.flag, cooksd.flag=cooksd.flag, cookthreshold=cookthreshold, weight.name=weight.name)

		auc.results[ ,i]   <- cv.result[["AUC"]]
		erbpp.results[ ,i] <- cv.result[["ER"]][["BPP"]]
		erbce.results[ ,i] <- cv.result[["ER"]][["BCE"]]
	}

	file.auc   <- sprintf("ulr_cv_model(%s)_family(%s)_auc.csv", model, family)
	file.erbpp <- sprintf("ulr_cv_model(%s)_family(%s)_erbpp.csv", model, family)
	file.erbce <- sprintf("ulr_cv_model(%s)_family(%s)_erbce.csv", model, family)

	auc.results   <- round(auc.results, 3)
	erbpp.results <- round(erbpp.results, 3)
	erbce.results <- round(erbce.results, 3)

	write.table(auc.results,   file=file.auc,   row.names=FALSE, col.names=TRUE, sep=",")
	write.table(erbpp.results, file=file.erbpp, row.names=FALSE, col.names=TRUE, sep=",")
	write.table(erbce.results, file=file.erbce, row.names=FALSE, col.names=TRUE, sep=",")
}

RemoveConfoundingEffect <- function(data, ynames, confeff.name="LOC")
{
	#### confeff.name: name of the variable in data that each of ynames have confounding effect with
	for (yname in ynames) {
		if (yname==confeff.name) { next }

		formula <- as.formula(paste(yname, confeff.name, sep=" ~ "))
		model <- glm(formula=formula, data=data, family=gaussian)
		data[, yname] <- data[, yname]-model$fitted.values
	}
	data
}

### split each class in data into number folds
SplitData <- function(data, yname, folds, seedval)
{
	classes <- unique(data[, yname]) ### multiple classes

	set.seed(seedval)
	seedval <- sample(1:1000000, 1) ### global seed

	data$fold <- 0
	for (aclass in classes) { ### 对数据集中每一类都分类fols等份
		index <- which(data[, yname]==aclass)

		set.seed(seedval)
		seq <- sample(1:length(index))

		data$fold[index] <- rep(1:folds, length.out=length(index))[seq]
	}
	return(data)
}

### sampling metod, including "under", "over", and "both" sampling method
### only for data set with binary classes
Sampling <- function(data, yname, sampling, seedval)
{
	dataP <- data[data[, yname]==1, ]
	dataN <- data[data[, yname]==0, ]
	len.Pos <- nrow(dataP)
	len.Neg <- nrow(dataN)

	set.seed(seedval)
	seedval <- sample(1:1000000, 1)

	data.new <- NULL
	if (sampling==FALSE || len.Pos==len.Neg) {
		data.new <- data
	} else if (tolower(sampling)=="under") {
		### under majority class
		if (len.Pos < len.Neg) { ### len.Neg is large
			set.seed(seedval)
			seq.Neg <- sample(1:len.Neg, len.Pos)

			dataP.new <- dataP
			dataN.new <- dataN[seq.Neg,]
			data.new  <- rbind(dataP.new, dataN.new)
		} else {              ### len.Pos is large
			set.seed(seedval)
			seq.Pos <- sample(1:len.Pos, len.Neg)

			dataP.new <- dataP[seq.Pos,]
			dataN.new <- dataN
			data.new <- rbind(dataP.new, dataN.new)
		}
	} else if (tolower(sampling)=="over") {
		if (len.Pos < len.Neg) { ### increase dataPos
			set.seed(seedval)
			seq.Pos <- sample(1:len.Neg, len.Neg) %% len.Pos
			seq.Pos[which(seq.Pos==0)] <- len.Pos ### where is equal to zero then set to last element

			dataP.new <- dataP[seq.Pos,]
			dataN.new <- dataN
			data.new  <- rbind(dataP.new, dataN.new)
		} else {              ### increase dataNeg
			set.seed(seedval)
			seq.Neg <- sample(1:len.Pos, len.Pos) %% len.Neg
			seq.Neg[which(seq.Neg==0)] <- len.Neg

			dataP.new <- dataP
			dataN.new <- dataN[seq.Neg,]
			data.new <- rbind(dataP.new, dataN.new)
		}
	} else if (tolower(sampling)=="both") {
		n1 <- n2 <- (len.Pos + len.Neg) / 2

		if (len.Pos < len.Neg) { ### increase dataP and decrease dataN
			set.seed(seedval)
			seq.Pos <- sample(1:n1, n1) %% len.Pos
			seq.Pos[which(seq.Pos==0)] <- len.Pos

			set.seed(seedval)
			seq.Neg <- sample(1:len.Neg, n2)

			dataP.new <- dataP[seq.Pos, ]
			dataN.new <- dataN[seq.Neg, ]
			data.new <- rbind(dataP.new, dataN.new)
		} else { ### increase dataN and decrease dataP
			set.seed(seedval)
			seq.Pos <- sample(1:len.Pos, n1)

			set.seed(seedval)
			seq.Neg <- sample(1:n2, n2) %% len.Neg
			seq.Neg[which(seq.Neg==0)] <- len.Neg

			dataP.new <- dataP[seq.Pos, ]
			dataN.new <- dataN[seq.Neg, ]
			data.new <- rbind(dataP.new, dataN.new)
		}
	}

	return(data.new)
}

### this function used for validating whether a independent variable is significantly associated with the dependent variable
isSignificantVariable <- function(data, xname, yname, pthreshold=0.05, cooksd.flag=FALSE, weight.flag=FALSE, weight.name="weight", cookthreshold=1, trace=FALSE)
{
	udata <- data
	fm <- as.formula(paste(yname, xname, sep=" ~ "))
	while(TRUE) {
		if (weight.flag) {
			ul <- glm(formula=fm, family=binomial, data=udata, weights=weight)
		} else {
			ul <- glm(formula=fm, family=binomial, data=udata)
		}

		if (cooksd.flag==FALSE) { break }
		cookd <- cooks.distance(ul)
		cookd[is.na(cookd)] <- cookthreshold + 1
		cookouts <- sum(cookd>=cookthreshold)
		if (cookouts==0) { break }
		udata <- udata[which(cookd<1), ]
	}

	flag <- FALSE
	ul.sum <- summary(ul)
	if (nrow(ul.sum$coefficients)==1) {
		if (trace) { cat("\n", xname, " not enter to the model.") }
	} else {
		pvalue <- ul.sum$coefficients[2, 4]
		if (pvalue <= pthreshold) {
			flag <- TRUE
		} else {
			if (trace) { cat("\n", xname, " p-value:", pvalue, " not significant for ulr.") }
		}
	}
	return(flag)
}

### using univariate logistic analysis to select those independent variables that significantly associated with dependent variable
Univariate.Select <- function(data, yname, xnames, pthreshold=0.05, cooksd.flag=FALSE, weight.flag=FALSE, weight.name="weight", cookthreshold=1, trace=FALSE, parallel=FALSE)
{
	flag.list <- NULL
	if (parallel && .Platform$OS.type=="unix") {
		require(parallel)
		numWorkers <- 4
		flag.list <- mclapply(X=xnames, FUN=isSignificantVariable, data=data, yname=yname, pthreshold=pthreshold, cooksd.flag=cooksd.flag, weight.flag=weight.flag, weight.name=weight.name, cookthreshold=cookthreshold, trace=trace, mc.cores=numWorkers)
	} else {
		flag.list <-   lapply(X=xnames, FUN=isSignificantVariable, data=data, yname=yname, pthreshold=pthreshold, cooksd.flag=cooksd.flag, weight.flag=weight.flag, weight.name=weight.name, cookthreshold=cookthreshold, trace=trace)
	}

	res.xnames <- xnames[unlist(flag.list)]
	return(res.xnames) ### return a list of names which are significant
}

### stepwise variable selection
Stepwise.One <- function(data, yname, xnames, rule="aic", direction="backward", weight.flag=FALSE, weight.name="weight", trace=FALSE)
{
	cat("step", rule, direction, "\n")

	rule <- tolower(rule)
	direction <- tolower(direction)

	null.formula <- as.formula(sprintf("%s ~ 1", yname))
	full.formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse=" + ")))

	model <- null.model <- full.model <- NULL
	if (weight.flag) {
		null.model <- glm(null.formula, family=binomial, data=data, weights=weight)
		full.model <- glm(full.formula, family=binomial, data=data, weights=weight)
	} else {
		null.model <- glm(null.formula, family=binomial, data=data)
		full.model <- glm(full.formula, family=binomial, data=data)
	}

	if (rule=="pvalue") {
		repeat{
			if (trace) { cat("repeat...\n") }
			ce <- full.model$coefficients[attr(full.model$terms, "term.labels")]
			if (any(is.na(ce))) {
				rterm <- names(ce)[!is.na(ce)]
				full.formula <- sprintf("%s ~ %s", yname, paste(rterm, collapse="+"))
				if (weight.flag) {
					full.model <- glm(as.formula(full.formula), family=binomial, data=data, weights=weight)
				} else {
					full.model <- glm(as.formula(full.formula), family=binomial, data=data)
				}
			} else {
				break
			}
		}

		if (direction=="forward") {
			#### TOBE deleted:
			model <- Forward.Stepwise(null.model, full.model, trace=trace)
		} else {
			#### TOBE deleted:
			model <- Backward.Stepwise.Custormized(full.model, trace=trace)
		}
	} else {
		if (rule=="aic") { k <- 2 } else { k <- log(nrow(data)) }
		scope <- list(upper=full.model, lower=null.model)
		if (direction=="forward") {
			model <- step(null.model, scope=scope, direction="forward",  k=k, trace=trace)
		} else if (direction=="backward") {
			model <- step(full.model, scope=scope, direction="backward", k=k, trace=trace)
		} else if (direction=="forwardboth") {
			model <- step(null.model, scope=scope, direction="both",     k=k, trace=trace)
		} else if (direction=="backwardboth") {
			model <- step(full.model, scope=scope, direction="both",     k=k, trace=trace)
		}
	}

	return(model)
}

### stepwise logistic regression
# stepwise <- function(data, yname, xnames, rule="aic", direction="backward", confounding.effect.name="LOC", vif.threshold=10, univar.flag=FALSE, cooksd.flag=FALSE, remove.flag=FALSE, weight.flag=FALSE, weight.name="weight", pthreshold=0.05, cookthreshold=1, trace=FALSE, parallel=FALSE)
# {
# 	if (remove.flag != FALSE) {
# 		if (trace) { cat("Remove LOC confounding for variable selection.\n") }
# 		data <- RemoveConfoundingEffect(data, ynames=xnames, confeff.name=confounding.effect.name)
# 	}
#
# 	if (univar.flag) {
# 		if (trace) { cat("perform univariate logistic regression check.\n") }
# 		if (length(xnames) != 1) {
# 			txnames <- Univariate.Select(data=data, yname=yname, xnames=xnames, pthreshold=pthreshold, cooksd.flag=cooksd.flag, weight.flag=weight.flag, weight.name=weight, parallel=parallel)
# 			if (length(txnames) >= 1) {
# 				xnames <- txnames
# 			} else {
# 				xnames <- xnames
# 				if (trace) { cat("after univariate variable selection, number of variables is zero, use the previous variables names instead.")	}
# 			}
# 		}
# 	}
#
# 	if (!is.na(vif.threshold)) { xnames <- as.character(usdm::vifstep(data[, xnames], th=vif.threshold)@results$Variables) } ### remove variable while vif is large than 10. or using code # xnames <- vif_func(data[, xnames])
# 	model <- stepwise.model(data=data, yname=yname, xnames=xnames, rule=rule, direction=direction, weight.flag=weight.flag, weight.name=weight.name, trace=trace)
# 	return(model)
# }

### stepwise logistic regression variable selection
Stepwise <- function(data, yname, xnames, confounding.effect.name="LOC", rule="aic", direction="backward", vif.threshold=10, univar.flag=FALSE, pthreshold=0.05, cooksd.flag=FALSE, cookthreshold=1, weight.flag=FALSE, remove.flag=FALSE, weight.name="weight", trace=FALSE)
{
	if (length(xnames)==0) { stop("Length of independent variable is zero. Please check.\n") }

	data.new <- data
	xnames.new <- xnames
	if (remove.flag==TRUE) {
		if (trace) { cat("Removing confounding effect for independent variables before stepwise variable selection.\n") }
		data.new <- RemoveConfoundingEffect(data.new, ynames=xnames.new, confeff.name=confounding.effect.name)
	}

	if (univar.flag) {
		if (trace) { cat("Removing independent variables that are not significantly associated with dependent variable before stepwise variable selection.\n") }
		xnames.new <- Univariate.Select(data=data.new, yname=yname, xnames=xnames.new, pthreshold=pthreshold, cooksd.flag=cooksd.flag, weight.flag=weight.flag, weight.name=weight, parallel=parallel)
	}

	if (!is.null(vif.threshold)) {
		if (trace) { cat("Removing independent variables step by step with an vif threshold:", vif.threshold, "before stepwise variable selection.\n") }
		xnames.new <- as.character(usdm::vifstep(data.new[, xnames.new], th=vif.threshold)@results$Variables) ### remove variable while vif is large than vif.threshold or using code
	}

	if (length(xnames.new) == 0) {
		cat("Before stepwise variable selection all independent variables are removed, we therefore use all indpendent variables instead.\n")
		xnames.new <- xnames
	}

	model <- NULL
	while(TRUE) {
		model <- Stepwise.One(data=data.new, yname=yname, xnames=xnames.new, rule=rule, direction=direction, weight.flag=weight.flag, weight.name=weight.name, trace=trace)
		break_flag <- FALSE
		if (!is.null(vif.threshold)) {
			xnames.new <- subsets <- attr(model$terms, "term.labels")
			err_catch <- try(varvifs <- abs(car::vif(model)))
			if (class(err_catch)=="try-error") { break_flag <- TRUE	} else {
				if ((length(subsets)>=2) && (max(varvifs) > vif.threshold)) {
					xnames.new <- setdiff(xnames.new, subsets[which.max(varvifs)])
					break_flag <- FALSE
				} else { break_flag <- TRUE }
			}
		} else { break_flag <- TRUE }

		# if (cooksd.flag) {
		# 	cookd <- cooks.distance(model)
		# 	cookd[is.na(cookd)] <- cookthreshold + 1
		# 	num_cookouts <- sum(cookd>=cookthreshold)
		# 	data.new <- data.new[which(cookd<cookthreshold), ]
		# 	if (num_cookouts!=0) {
		# 		if (trace) { cat("Remove", num_cookouts, "(total remove:", nrow(data) - nrow(model$data), ") cases which cook distance larger than 1, redo variable selection.\n") }
		# 		break_flag <- FALSE
		# 	} else { break_flag <- TRUE }
		# } else { break_flag <- TRUE }

		if (break_flag) { break }
	}

	if (cooksd.flag) {
		cookd <- cooks.distance(model)
		cookd[is.na(cookd)] <- cookthreshold + 1
		num_cookouts <- sum(cookd>=cookthreshold)
		data.new <- data.new[which(cookd<cookthreshold), ]
		model <- glm(formula=model$formula, family=model$family, data=data.new)
		if ((num_cookouts!=0) & trace) { cat("Remove", num_cookouts, "(total remove:", nrow(data) - nrow(model$data), ") cases which cook distance larger than 1") }
	}

	# require(car) ### vif function
	# subsets <- attr(model$terms, "term.labels")
	# if (length(subsets) >= 2) {
	# 	vifs <- abs(car::vif(model))
	# 	while( (length(subsets) >= 2) && (max(vifs) > vif.threshold) ) {
	# 		if (length(subsets) >= 2) {
	# 			if (trace) { cat("Remove variable ", subsets[which.max(vifs)], " from model then redo variable selection.\n") }
	# 			xnames <- setdiff(xnames, subsets[which.max(vifs)])
	# 			model <- stepwise(data=data, yname=yname, xnames=xnames, rule=rule, direction=direction, vif.threshold=vif.threshold, univar.flag=univar.flag, cooksd.flag=cooksd.flag, weight.flag=weight.flag, weight.name="weight", trace=trace)
	# 			subsets <- attr(model$terms, "term.labels")
	# 			if (length(subsets) > 1) {
	# 				vifs <- abs(car::vif(model))
	# 			}
	# 		}
	# 	}
	# }

	return(model)
}

### this function sort the data frame base on the predicted value and second by other policies
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

#####################################################################
# AUC: area under the ROC curve.  a traditional classification criteria
# input:
# data.frame should has the following two columns
# column "PRE": the predict value
# column "REL": the real class, it's a binary value (have bugs or have no bugs)
#####################################################################
### optimized by Yibiao Yang on 2014-10-16
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

### compute area under curve for CE
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

##########################################################
### optimized by Yibiao Yang on 2014-09-11
### compute threshold or cutoff for a given (training) data set
### e.g. BPP(balance pf-pd)
###      BCE(balance classification error)
###      MFM(maximum f-measure)
##################################################
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

#####################################################################
### ER : effort reduction compared to random model (model evaluation criteria for cost-effectiveness)
#####################################################################
ComputeER <- function(data, threshold=NULL, cutoff=NULL, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
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

	totalLOC   <- sum(data$LOC)
	totalLOC.M <- sum(data$PREBIN*data$LOC)

	totalNUM   <- sum(data$NUM)
	totalNUM.R <- sum(data$PREBIN*data$REL*data$NUM)

	if (totalNUM==0) { stop("of this testing data set, total number of NUM is zero") }
	if (totalLOC==0) { stop("of this testing data set, total number of LOC is zero") }

	Effort.M <- totalLOC.M / totalLOC ### effort of model m
	Effort.R <- totalNUM.R / totalNUM ### effort of random model

	ER.TNE <- sum((1-data$PREBIN)*(1-data$REL)*data$LOC)/sum(data$LOC) ### TP/FP肯定要审查, FN部分在运行期间引发失效后也会被审查. 与全部审查模型相比, 模型m节省的工作量约简为 ER' = [Effort(all) - Effort(m)] / Effort(all) = SLOC(TN) / SLOC(TP+FP+FN+TN)
	ER.ABS.REL <- Effort.R - Effort.M ### 绝对的ER 相对于整个系统代码审查的百分比
	ER.ABS.LOC <- ER.ABS.REL*totalLOC ### 绝对的ER 相对于整个系统代码审查的百分比
	ER.REL.RAND <- (Effort.R - Effort.M) / Effort.R ### 相对于随机模型的ER
	ER.REL.SELF <- (Effort.R - Effort.M) / Effort.M ### 相对于自身模型的ER

	if (Effort.R==0) { ER.REL.RAND <- NA }
	if (Effort.M==0) { ER.REL.SELF <- NA }

	return(list(ER.TNE=ER.TNE, ER.ABS.REL=ER.ABS.REL, ER.ABS.LOC=ER.ABS.LOC, ER.REL.RAND=ER.REL.RAND, ER.REL.SELF=ER.REL.SELF, ER.EFFORT.SELF=Effort.M, ER.EFFORT.RAND=Effort.R))
}

#####################################################################
### ERAVG: optimization by Yibiao Yang 2014-09-11
### all possible ER
#####################################################################
ComputeERAVG <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE, allpercentcutoff=TRUE)
{
	if (!sorted) {
		data <- SortData(data=data, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
		sorted <- TRUE
	}

	len <- nrow(data)
	cumLOC <- cumsum(data$LOC)
	cumNUM <- cumsum(data$NUM)

	Effort.Ms.old <- cumLOC / cumLOC[len] ### 模型m的工作量
	Effort.Rs.old <- cumNUM / cumNUM[len] ### 随机模型的工作量

	Effort.Ms <- Effort.Ms.old
	Effort.Rs <- Effort.Rs.old
	if (!allpercentcutoff) {
		pos.equal <- which(data[1:(len-1), "PRE"]==data[2:len, "PRE"])
		if (length(pos.equal)!=0) {
			Effort.Ms <- Effort.Ms.old[-pos.equal]
			Effort.Rs <- Effort.Rs.old[-pos.equal]
		}
	}

	ERs.TNE <- cumsum((1-rev(data$REL))*data$LOC)/sum(data$LOC)
	ERs.ABS.REL <- Effort.Rs - Effort.Ms
	ERs.ABS.LOC <- ERs.ABS.REL*cumLOC[len]
	ERs.REL.RAND <- (Effort.Rs - Effort.Ms) / Effort.Rs
	ERs.REL.SELF <- (Effort.Rs - Effort.Ms) / Effort.Ms

	ERs.REL.RAND[Effort.Rs==0] <- NA
	ERs.REL.SELF[Effort.Ms==0] <- NA

	AVG.TNE <- mean(ERs.TNE)
	AVG.ABS.REL <- mean(ERs.ABS.REL)
	AVG.ABS.LOC <- mean(ERs.ABS.LOC)
	AVG.REL.RAND <- mean(ERs.REL.RAND, na.rm=TRUE)
	AVG.REL.SELF <- mean(ERs.REL.SELF, na.rm=TRUE)

	return(list(ER.AVG.TNE=AVG.TNE, ER.AVG.ABS.REL=AVG.ABS.REL, ER.AVG.ABS.LOC=AVG.ABS.LOC, ER.AVG.REL.RAND=AVG.REL.RAND, ER.AVG.REL.SELF=AVG.REL.SELF, ER.AVG.E.SELF=mean(Effort.Ms), ER.AVG.E.RAND=mean(Effort.Rs)))
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

### compute area under alberg curve
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

ComputePopt <- function(data, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
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

	Popt <- 1 - (area.opt - area.mdl)/(area.opt - area.wst)

	return(Popt)
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

### compute ers include ER-BPP ER-BCE ER-MFM ERAVG ER0.5
### fix by Yibiao Yang on 2013-07-08
ComputeERs <- function(train.dt, valid.dt, thresholds=NULL, effortaware=TRUE, sorted=FALSE, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE, allpercentcutoff=TRUE)
{
	if (!sorted) {
		train.dt <- SortData(data=train.dt, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
		valid.dt <- SortData(data=valid.dt, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
		sorted <- TRUE
	}

	if (is.null(thresholds)) {
		thresholds <- ComputeThreshold(data=train.dt, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)
	}

	ret <- NULL
	thresh.list.names <- c("V.BPP", "V.BCE", "V.MFM")
	for (thresh.list.name in thresh.list.names) {
		temp <- ComputeER(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, threshold=thresholds[[thresh.list.name]])
		names(temp) <- paste(names(temp), thresh.list.name, sep=".")
		ret <- c(ret, unlist(temp))
	}

	temp <- ComputeER(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, threshold=0.5)
	names(temp) <- paste(names(temp), "V.0.5", sep=".")
	ret <- c(ret, unlist(temp))

	cutoff.list.names <- c("P.BPP", "P.BCE", "P.MFM")
	for (cutoff.list.name in cutoff.list.names) {
		temp <- ComputeER(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, cutoff=thresholds[[cutoff.list.name]])
		names(temp) <- paste(names(temp), cutoff.list.name, sep=".")
		ret <- c(ret, unlist(temp))
	}

	cutoff.list.names <- c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1))
	for (cutoff.list.name in cutoff.list.names) {
		temp <- ComputeER(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, cutoff=as.numeric(cutoff.list.name))
		names(temp) <- paste(names(temp), "P", cutoff.list.name, sep=".")
		ret <- c(ret, unlist(temp))
	}

	temp <- ComputeERAVG(data=valid.dt, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)
	ret <- c(ret, unlist(temp))

	ret <- unlist(ret)
	return(ret)
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

Compute <- function(data.train, data.valid, yname, xnames, run=0, fold=0, model="lr", formula=NULL, family="binomial", effort.name="LOC", number.name="BUGS", cutoffs=c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1)), sampling=FALSE, cooksd.flag=FALSE, weight.flag=FALSE, cookthreshold=1, weight.name="weight", trace=FALSE)
{
	if (trace) { cat("cross validation at (run: ", run, ", fold: ", fold, ")", "\n") }

	data.new.valid <- data.valid
	data.new.train <- data.train

	### sampling parameter can be either FALSE or "under" or "over" or "both"
	if (sampling!=FALSE) {
		set.seed(fold)
		seedval <- sample(1:1000000, 1)
		data.new.train <- Sampling(data=data.new.train, yname=yname, sampling=sampling, seedval=seedval)
	}

	### Note: do not need to select variable agian(It was determined by previous selection
	### based on previous model formula and train dataset to obtain model, and valid it by valid dataset
	valid.dt.normal <- train.dt.normal <- valid.dt.divide <- train.dt.divide <- NULL
	valid.dt.normal.non_effort <- train.dt.normal.non_effort <- valid.dt.divide.non_effort <- train.dt.divide.non_effort <- NULL
	valid.pre <- train.pre <- NULL

	if (is.null(formula)) { formula <- as.formula(sprintf("%s ~ %s", yname, paste(xnames, collapse="+"))) }

	model <- tolower(model)
	if (model=="manualup") { ####
		valid.dt.normal <- data.frame(REL=data.new.valid[, yname], PRE=1/(data.new.valid[, xname] + 1), LOC=data.new.valid[, effort.name], NUM=data.new.valid[, number.name])
		train.dt.normal <- data.frame(REL=data.new.train[, yname], PRE=1/(data.new.train[, xname] + 1), LOC=data.new.train[, effort.name], NUM=data.new.train[, number.name])

		valid.dt.divide <- data.frame(REL=data.new.valid[, yname], PRE=1/(data.new.valid[, xname] + 1)/(data.new.valid[, effort.name] + 1), LOC=data.new.valid[, effort.name], NUM=data.new.valid[, number.name])
		train.dt.divide <- data.frame(REL=data.new.train[, yname], PRE=1/(data.new.train[, xname] + 1)/(data.new.train[, effort.name] + 1), LOC=data.new.train[, effort.name], NUM=data.new.train[, number.name])

		valid.dt.normal.non_effort <- data.frame(REL=data.new.valid[, yname], PRE=data.new.valid[, xname], LOC=1, NUM=data.new.valid[, number.name])
		train.dt.normal.non_effort <- data.frame(REL=data.new.train[, yname], PRE=data.new.train[, xname], LOC=1, NUM=data.new.train[, number.name])

		valid.dt.divide.non_effort <- valid.dt.normal.non_effort
		train.dt.divide.non_effort <- train.dt.normal.non_effort
	} else if (model=="manualdown") { #####
		valid.dt.normal <- data.frame(REL=data.new.valid[, yname], PRE=data.new.valid[, xname], LOC=data.new.valid[, effort.name], NUM=data.new.valid[, number.name])
		train.dt.normal <- data.frame(REL=data.new.train[, yname], PRE=data.new.train[, xname], LOC=data.new.train[, effort.name], NUM=data.new.train[, number.name])

		valid.dt.divide <- data.frame(REL=data.new.valid[, yname], PRE=data.new.valid[, xname]/(data.new.valid[, effort.name] + 1), LOC=data.new.valid[, effort.name], NUM=data.new.valid[, number.name])
		train.dt.divide <- data.frame(REL=data.new.train[, yname], PRE=data.new.train[, xname]/(data.new.train[, effort.name] + 1), LOC=data.new.train[, effort.name], NUM=data.new.train[, number.name])

		valid.dt.normal.non_effort <- data.frame(REL=data.new.valid[, yname], PRE=data.new.valid[, xname], LOC=1, NUM=data.new.valid[, number.name])
		train.dt.normal.non_effort <- data.frame(REL=data.new.train[, yname], PRE=data.new.train[, xname], LOC=1, NUM=data.new.train[, number.name])

		valid.dt.divide.non_effort <- valid.dt.normal.non_effort
		train.dt.divide.non_effort <- train.dt.normal.non_effort
	} else if (model=="lr") {
		mdl <- NULL
		if (weight.flag) {
			mdl <- glm(formula=formula, family=family, data=data.new.train, weights=data.new.train[, weight.name])
		} else {
			mdl <- glm(formula=formula, family=family, data=data.new.train)
		}

		# if (cooksd.flag) {
		# 	cookd <- cooks.distance(mdl)
		# 	cookd[is.na(cookd)] <- cookthreshold + 1
		# 	num_cookouts <- sum(cookd>=cookthreshold)
		# 	if (num_cookouts!=0) {
		# 		data.new.train <- data.new.train[which(cookd<cookthreshold), ]
		# 		mdl <- glm(formula=mdl$formula, family=mdl$family, data=data.new.train)
		# 	}
		# }

		valid.pre <- predict(mdl, data.new.valid, type="response")
		train.pre <- mdl$fitted.values
	} else if (model=="multinomial") {
		data.new.train[, yname] <- as.factor(data.new.train[, yname])
		mdl <- NULL
		if (weight.flag) {
			mdl <- multinom(formula=formula, data=data.new.train, weights=weight, trace=FALSE)
		} else {
			mdl <- multinom(formula=formula, data=data.new.train, trace=FALSE)
		}

		train.pre <- predict(mdl, data.new.train, type="class")
		valid.pre <- predict(mdl, data.new.valid, type="class")

		classes <- unique(data.new.valid[, yname]) ### multiple classes

		recall    <- 0
		precision <- 0
		for (aclass in sort(classes)) {
			index1     <- data.new.valid[, yname]==aclass
			sdata      <- data.new.valid[index1, ]
			svalid.pre <- valid.pre[index1]
			srecall <- sum(svalid.pre==aclass)/sum(index1)
			recall  <- recall + 1/length(classes)*srecall

			index2     <- valid.pre==aclass
			sdata      <- data.new.valid[index2, ]
			svalid.pre <- valid.pre[index2]
			sprecision <- sum(sdata[, yname]==aclass)/sum(index2)

			if (!is.na(sprecision)) {
				precision <- precision + 1/length(classes)*sprecision
			}

			cat(aclass, " precision:", sprecision, " recall:", srecall, "\n")
		}
		cat("\n")

		accuracy <- sum(valid.pre==data.new.valid[, yname]) / nrow(data.new.valid)

		return(list(ACCURACY=accuracy, RECALL=recall, PRECISION=precision))
	} else if (model=="glmnet") {
		mdl <- glmnet::glmnet(x=as.matrix(data.new.train[, xnames]), y=as.vector(data.new.train[, yname]), family="binomial")
		valid.pre <- predict(mdl, newx=as.matrix(data.new.valid[, c(xnames)]), s=0.001, type="response")[, "1"]
		train.pre <- predict(mdl, newx=as.matrix(data.new.train[, c(xnames)]), s=0.001, type="response")[, "1"]
	} else if (model=="ridge") {
		mdl <- ridge::logisticRidge(formula=formula, data=data.new.train[, c(xnames, yname)], lambda="automatic")
		valid.pre <- predict(mdl, data.new.valid, type="response")
		train.pre <- predict(mdl, data.new.train, type="response")
	} else if (model=="naivebayes") { # naive bayes classifier in package klaR
		# require(e1071)
		mdl <- e1071::naiveBayes(formula=formula, data=data.new.train[,c(xnames, yname)])
		valid.pre <- predict(mdl, data.new.valid, type="raw")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="raw")[, "1"]
	} else if (model=="nnet") { # single-hidden-layer neural network, possibly with skip-layer connections in package nnet
		# require(nnet)
		mdl <- nnet::nnet(formula=formula, data=data.new.train[, c(xnames, yname)], size=2)
		valid.pre <- predict(mdl, data.new.valid, type="raw")
		train.pre <- predict(mdl, data.new.train, type="raw")
	} else if (model=="neuralnet") { # neural networks using backpropagation in package neuralnet
		# require(neuralnet)
		tdata.new.train <- data.new.train
		tdata.new.train[,yname] <- as.numeric(as.character(tdata.new.train[, yname]))
		mdl <- neuralnet::neuralnet(formula=formula, data=tdata.new.train[, c(xnames, yname)], linear.output=FALSE)
		valid.pre <- neuralnet::compute(mdl, data.new.valid[, xnames])$net.result
		train.pre <- neuralnet::compute(mdl, data.new.train[, xnames])$net.result
	} else if (model=="mlp") { # multi-layer perceptron in RSNNS
		# require(RSNNS)
		tdata.new.train <- data.new.train
		tdata.new.train[,yname] <- as.numeric(as.character(tdata.new.train[,yname]))
		mdl <- RSNNS::mlp(x=as.matrix(tdata.new.train[,xnames]), y=as.matrix(tdata.new.train[,yname]))
		valid.pre <- predict(mdl, data.new.valid, type="prob")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="prob")[, "1"]
	} else if (model=="monmlp") { # multi-layer perceptron neural network in package monmlp
		# require(monmlp)
		tdata.new.train <- data.new.train
		tdata.new.train[, yname] <- as.numeric(as.character(tdata.new.train[, yname]))
		mdl <- monmlp::monmlp.fit(x=as.matrix(tdata.new.train[, xnames]), y=as.matrix(tdata.new.train[, yname]), hidden1=2, silent=TRUE)
		valid.pre <- monmlp.predict(x=as.matrix(data.new.valid[, xnames]), weights=mdl)
		train.pre <- monmlp.predict(x=as.matrix(data.new.train[, xnames]), weights=mdl)
	} else if (model=="rpart") { # Recursive Partitioning and Regression Trees in package rpart
		# require(rpart)
		mdl <- rpart::rpart(formula=formula, data=data.new.train[,c(xnames, yname)])
		valid.pre <- predict(mdl, newdata=data.new.valid, type="prob")[, "1"]
		train.pre <- predict(mdl, newdata=data.new.train, type="prob")[, "1"]
	} else if (model=="randomforest") { # Classi?cation and Regression with Random Forest in package randomForest
		# require(randomForest)
		mdl <- randomForest::randomForest(formula=formula, data = data.new.train[, c(xnames, yname)])
		valid.pre <- predict(mdl, newdata=data.new.valid, type="prob")[, "1"]
		train.pre <- predict(mdl, newdata=data.new.train, type="prob")[, "1"]
	} else if (model=="knn") { # k-Nearest Neighbour Classification in package class
		# require(class)
		valid.pre <- attr(class::knn(train=data.new.train[,c(xnames)], test=data.new.valid[,c(xnames)], cl=data.new.train[,c(yname)], k=10, prob=TRUE), "prob")
		train.pre <- attr(class::knn(train=data.new.train[,c(xnames)], test=data.new.train[,c(xnames)], cl=data.new.train[,c(yname)], k=10, prob=TRUE), "prob")
	} else if (model=="kknn") { # Weighted k-Nearest Neighbor Classifier in package kknn
		# require(kknn)
		valid.pre <- kknn::kknn(formula=formula, train = data.new.train[, c(xnames, yname)], test = data.new.valid[,c(xnames)], k=10)$prob[,"1"]
		train.pre <- kknn::kknn(formula=formula, train = data.new.train[, c(xnames, yname)], test = data.new.train[,c(xnames)], k=10)$prob[,"1"]
	} else if (model=="svm") { # support vector machine in package e1071
		# require(e1071)
		mdl <- e1071::svm(formula=formula, data=data.new.train[, c(xnames, yname)], probability=TRUE) ### type="C-classification", kernel="radial", cost=10, gamma=0.1)
		valid.pre <- attr(predict(mdl, data.new.valid, probability = TRUE), "probabilities")[, "1"]
		train.pre <- attr(predict(mdl, data.new.train, probability = TRUE), "probabilities")[, "1"]
	} else if (model=="ksvm") { # Support Vector Machines in package kernlab
		# require(kernlab)
		mdl <- kernlab::ksvm(formula=formula, data = data.new.train[,c(xnames, yname)], type="C-bsvc", kernel="rbfdot", kpar=list(sigma=0.1), C=10, prob.model=TRUE)
		valid.pre <- predict(mdl, data.new.valid, type="probabilities")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probabilities")[, "1"]
	} else if (model=="tree") { # Fit a Classi?cation or Regression Tree in pakcage tree
		# require(tree)
		mdl <- tree::tree(formula=formula, data = data.new.train[,c(xnames, yname)])
		valid.pre <- predict(mdl, data.new.valid)[,"1"]
		train.pre <- predict(mdl, data.new.train)[,"1"]
	}

	##############################
	####
	#### RWeka classifiers
	####
	##############################
	#### function learners
	else if (model=="linearregression") {
		# require(RWeka)
		mdl <- RWeka::LinearRegression(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="logistic") {
		# require(RWeka)
		mdl <- RWeka::Logistic(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="smo") {
		# require(RWeka)
		mdl <- RWeka::SMO(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	}

	#### lazy learners
	else if (model=="ibk") {
		# require(RWeka)
		mdl <- RWeka::IBk(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="lbr") { #### cannot handle numeric attribute
		# require(RWeka)
		mdl <- RWeka::LBR(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	}

	#### meta learners
	else if (model=="adaboostm1") {
		# require(RWeka)
		mdl <- RWeka::AdaBoostM1(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="bagging") {
		# require(RWeka)
		mdl <- RWeka::Bagging(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="logitboost") {
		# require(RWeka)
		mdl <- RWeka::LogitBoost(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="multiboostab") {
		# require(RWeka)
		mdl <- RWeka::MultiBoostAB(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="stacking") {
		# require(RWeka)
		mdl <- RWeka::Stacking(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="costsensitiveclassifier") {
		# require(RWeka)
		mdl <- RWeka::CostSensitiveClassifier(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	}

	#### rule learners
	else if (model=="jrip") {
		# require(RWeka)
		mdl <- RWeka::JRip(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="m5rules") {
		# require(RWeka)
		mdl <- RWeka::M5Rules(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="oner") {
		# require(RWeka)
		mdl <- RWeka::OneR(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="part") {
		# require(RWeka)
		mdl <- RWeka::PART(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	}

	#### tree learners
	else if (model=="j48") {
		# require(RWeka)
		formula <- as.formula(sprintf("as.factor(%s) ~ %s", yname, paste(xnames, collapse="+")))
		mdl <- RWeka::J48(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="lmt") {
		# require(RWeka)
		mdl <- RWeka::LMT(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="m5p") {
		# require(RWeka)
		mdl <- RWeka::M5P(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	} else if (model=="decisionstump") {
		# require(RWeka)
		mdl <- RWeka::DecisionStump(formula=formula, data=data.new.train[,c(xnames,yname)])
		valid.pre <- predict(mdl, data.new.valid, type="probability")[, "1"]
		train.pre <- predict(mdl, data.new.train, type="probability")[, "1"]
	}

	### normal (effort-aware)
	valid.dt.normal.effort <- data.frame(REL=data.new.valid[,yname], PRE=valid.pre, LOC=data.new.valid[, effort.name], NUM=data.new.valid[, number.name])
	train.dt.normal.effort <- data.frame(REL=data.new.train[,yname], PRE=train.pre, LOC=data.new.train[, effort.name], NUM=data.new.train[, number.name])

	### divide (effort-aware)
	valid.dt.divide.effort <- data.frame(REL=data.new.valid[,yname], PRE=valid.pre/(data.new.valid[, effort.name] + 1), LOC=data.new.valid[, effort.name], NUM=data.new.valid[, number.name])
	train.dt.divide.effort <- data.frame(REL=data.new.train[,yname], PRE=train.pre/(data.new.train[, effort.name] + 1), LOC=data.new.train[, effort.name], NUM=data.new.train[, number.name])

	### normal (non effort-aware)
	valid.dt.normal.non_effort <- data.frame(REL=data.new.valid[,yname], PRE=valid.pre, LOC=1, NUM=data.new.valid[, number.name])
	train.dt.normal.non_effort <- data.frame(REL=data.new.train[,yname], PRE=train.pre, LOC=1, NUM=data.new.train[, number.name])

	### divide (non effort-aware)
	valid.dt.divide.non_effort <- data.frame(REL=data.new.valid[,yname], PRE=valid.pre/(data.new.valid[, effort.name] + 1), LOC=1, NUM=data.new.valid[, number.name])
	train.dt.divide.non_effort <- data.frame(REL=data.new.train[,yname], PRE=train.pre/(data.new.train[, effort.name] + 1), LOC=1, NUM=data.new.train[, number.name])

	if (!is.numeric(valid.dt.normal$REL)) {
		train.dt.normal.effort$REL <- as.numeric(as.character(train.dt.normal.effort$REL))
		valid.dt.normal.effort$REL <- as.numeric(as.character(valid.dt.normal.effort$REL))

		train.dt.divide.effort$REL <- as.numeric(as.character(train.dt.divide.effort$REL))
		valid.dt.divide.effort$REL <- as.numeric(as.character(valid.dt.divide.effort$REL))

		train.dt.normal.non_effort$REL <- as.numeric(as.character(train.dt.normal.non_effort$REL))
		valid.dt.normal.non_effort$REL <- as.numeric(as.character(valid.dt.normal.non_effort$REL))

		train.dt.divide.non_effort$REL <- as.numeric(as.character(train.dt.divide.non_effort$REL))
		valid.dt.divide.non_effort$REL <- as.numeric(as.character(valid.dt.divide.non_effort$REL))
	}

	if (any(is.na(valid.pre))) { stop("valid.pre is NA") }
	if (any(is.na(train.pre))) { stop("train.pre is NA") }

	subcompute <- function(train.dt, valid.dt, worstcase=FALSE, bestcase=FALSE, LOCUP=FALSE)
	{
		train.dt.nea <- SortData(data=train.dt, effortaware=FALSE, sorted=FALSE) ### non effort-aware
		valid.dt.nea <- SortData(data=valid.dt, effortaware=FALSE, sorted=FALSE) ### non effort-aware

		train.dt.ea <- SortData(data=train.dt, effortaware=TRUE, sorted=FALSE) ### effort-aware
		valid.dt.ea <- SortData(data=valid.dt, effortaware=TRUE, sorted=FALSE) ### effort-aware

		sorted <- TRUE
		allpercentcutoff <- TRUE

		effortaware <- TRUE
		ARCE <- ComputeCEs(data=valid.dt.ea, cutoffs=cutoffs, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP)
		thresholds.ea <- ComputeThreshold(data=train.dt.ea, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)
		ER <- ComputeERs(train.dt=train.dt.ea, valid.dt=valid.dt.ea, thresholds=thresholds.ea, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)

		effortaware <- FALSE
		thresholds.nea <- ComputeThreshold(data=train.dt.nea, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)
		CP <- ComputeCPs(train.dt=train.dt.nea, valid.dt=valid.dt.nea, thresholds=thresholds.nea, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)

		AUC <- ComputeAUC(data=valid.dt.nea, effortaware=effortaware, sorted=sorted, worstcase=worstcase, bestcase=bestcase, LOCUP=LOCUP, allpercentcutoff=allpercentcutoff)

		RECALL <- CP$RECALL
		PRECISION <- CP$PRECISION
		F1 <- CP$F1
		ACCURACY <- CP$ACCURACY

		ER <- c(ARCE$ER, ER)
		RECALL <- c(ARCE$RECALL, RECALL)
		PRECISION <- c(ARCE$PRECISION, PRECISION)
		ACCURACY <- c(ARCE$ACCURACY, ACCURACY)
		F1 <- c(ARCE$F1, F1)

		AUC <- rbind(AUC)
		PRECISION <- rbind(PRECISION)
		RECALL <- rbind(RECALL)
		ACCURACY <- rbind(ACCURACY)
		F1 <- rbind(F1)
		ER <- rbind(ER)
		CE <- rbind(ARCE$CE)

		return(list(AUC=AUC, PRECISION=PRECISION, RECALL=RECALL, ACCURACY=ACCURACY, F1=F1, ER=ER, CE=CE))
	}

	normal.effort <- subcompute(train.dt=train.dt.normal.effort, valid.dt=valid.dt.normal.effort) ### column LOC with different values
	divide.effort <- subcompute(train.dt=train.dt.divide.effort, valid.dt=valid.dt.divide.effort) ### column LOC with different values

	normal.non_effort <- subcompute(train.dt=train.dt.normal.non_effort, valid.dt=valid.dt.normal.non_effort) ### column LOC with same value 1
	divide.non_effort <- subcompute(train.dt=train.dt.divide.non_effort, valid.dt=valid.dt.divide.non_effort) ### column LOC with same value 1

	return(list(normal.effort=normal.effort, divide.effort=divide.effort, normal.non_effort=normal.non_effort, divide.non_effort=divide.non_effort))
}

### one fold of cross-validation
OneFold <- function(data, yname, xnames, run=0, fold=0, model="lr", formula=NULL, family="binomial", effort.name="LOC", number.name="BUGS", cutoffs=c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1)), sampling=FALSE, cooksd.flag=FALSE, remove.at.fold.flag=FALSE, remove.flag=FALSE, weight.flag=FALSE, cookthreshold=1, weight.name="weight", trace=FALSE)
{
	data.valid <- data[data$fold==fold, ]
	data.train <- data[data$fold!=fold, ]

	if (remove.at.fold.flag) { ### default is FALSE
		if (remove.flag=="separate") {
			if (length(xnames)!=0) {
				data.valid <- RemoveConfoundingEffect(data.valid, ynames=xnames, confeff.name=effort.name)
				data.train <- RemoveConfoundingEffect(data.train, ynames=xnames, confeff.name=effort.name)
			}
		} else if (remove.flag=="together") {
			if (length(xnames)!=0) {
				data.new <- RemoveConfoundingEffect(data, ynames=xnames, confeff.name=effort.name)
				data.valid <- data.new[data.new$fold==fold, ]
				data.train <- data.new[data.new$fold!=fold, ]
			}
		}
	}

	Compute(data.train=data.train, data.valid=data.valid, fold=fold, model=model, formula=formula, family=family, yname=yname, xnames=xnames, effort.name=effort.name, number.name=number.name, run=run, sampling=sampling, cutoffs=cutoffs, cooksd.flag=cooksd.flag, cookthreshold=1, weight.flag=weight.flag, weight.name=weight.name, trace=trace)
}

### one time of XXX-fold cross-validation
OneRun <- function(run, data, yname, xnames, totalRuns=10, totalFolds=10, model="lr", formula=NULL, family="binomial", effort.name="LOC", number.name="BUGS", cutoffs=c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1)), sampling=FALSE, cooksd.flag=FALSE, remove.at.fold.flag=FALSE, remove.flag=FALSE, weight.flag=FALSE, cookthreshold=1, weight.name="weight", parallel=FALSE, trace=FALSE)
{
	seed <- totalRuns * 100 + totalFolds * 10 + run
	set.seed(seed)
	seedval  <- sample(1:1000000, 1)

	data.new <- SplitData(data=data, yname=yname, folds=totalFolds, seedval=seedval)

	if (trace) { cat("cross validation at run: ", run, "(total runs:", totalRuns, ") for each of fold(total folds:", totalFolds, ").\n") }

	resList.FOLDS <- NULL
	if (parallel && (.Platform$OS.type=="unix")) {
		require(parallel) ### this is for multi-threading
		numWorkers <- 4
		resList.FOLDS <- mclapply(c(1:totalFolds), FUN=OneFold, data=data.new, model=model, formula=formula, family=family, yname=yname, xnames=xnames, effort.name=effort.name, number.name=number.name, run=run, remove.at.fold.flag=remove.at.fold.flag, remove.flag=remove.flag, sampling=sampling, cutoffs=cutoffs, cooksd.flag=cooksd.flag, cookthreshold=cookthreshold, weight.flag=weight.flag, weight.name=weight.name, trace=trace, mc.cores=numWorkers)
	} else {
		resList.FOLDS <-   lapply(c(1:totalFolds), FUN=OneFold, data=data.new, model=model, formula=formula, family=family, yname=yname, xnames=xnames, effort.name=effort.name, number.name=number.name, run=run, remove.at.fold.flag=remove.at.fold.flag, remove.flag=remove.flag, sampling=sampling, cutoffs=cutoffs, cooksd.flag=cooksd.flag, cookthreshold=cookthreshold, weight.flag=weight.flag, weight.name=weight.name, trace=trace)
	}

	names(resList.FOLDS) <- paste("FOLD_", 1:totalFolds, sep="")
	return(resList.FOLDS)
}

###########################################################################
### cross-validation for machine learning models or logistic regression model
### remove flag can be "separate", "together", FALSE which is implemented in the 'Compute' function
######## parameters description ########
### data: data frame {column names with yname(dependent variable name) and xnames(indpendent variable names)}
### yname: dependent variable name, e.g. yname="buggy" indicate whether a module is buggy or not
### xnames: independent variables names, e.g. xnames=c("x1", "x2", "x3")
### effort.name: the proxy for computing cost-effectiveness, e.g. effort.name="LOC" indicate the number lines of code (effort.name should be a column in the 'data' data.frame)
### number.name: the number of bugs, e.g. number.name="bugs" indicate the number of bugs in a class (number.name should in the 'data' data.frame)
### totalFolds: number of folds, e.g. if m*n folds, the totalFolds indicate n
### totalRuns: number of runs, e.g. if m*n folds, the totalRuns is m
### others could be set as default
###########################################################################
CrossValidation <- function(data, yname, xnames, totalFolds=10, totalRuns=10, model="lr", formula=NULL, family="binomial", effort.name="LOC", number.name="BUGS", cutoffs=c(0.01, 0.05, 0.1, 0.15, seq(from=0.2, to=1.0, by=0.1)), sampling=FALSE, cooksd.flag=FALSE, remove.at.fold.flag=FALSE, remove.flag=FALSE, weight.flag=FALSE, cookthreshold=1, weight.name="weight", parallel=FALSE, trace=FALSE)
{
	### runs n times m-fold cross-validation
	if (remove.at.fold.flag==FALSE) {
		### if not remove at fold thus remove confounding effect for the entire data set
		if ((remove.flag!=FALSE)&&(length(xnames)!=0)) {
			data <- RemoveConfoundingEffect(data, ynames=xnames, confeff.name=effort.name)
		}
	}

	if (trace) { cat("Cross Validation for each run.\n") }

	resList.RUNS <- NULL
	if (parallel && (.Platform$OS.type=="unix")) {
		### this is for multi-threading
		require(parallel)
		numWorkers <- 6
		resList.RUNS <- mclapply(c(1:totalRuns), FUN=OneRun, data=data, model=model, formula=formula, family=family, yname=yname, xnames=xnames, effort.name=effort.name, number.name=number.name, totalRuns=totalRuns, totalFolds=totalFolds, remove.at.fold.flag=remove.at.fold.flag, remove.flag=remove.flag, sampling=sampling, cutoffs=cutoffs, cooksd.flag=cooksd.flag, cookthreshold=cookthreshold, weight.flag=weight.flag, weight.name=weight.name, parallel=parallel, trace=trace, mc.cores=numWorkers)
	} else {
		resList.RUNS <-   lapply(c(1:totalRuns), FUN=OneRun, data=data, model=model, formula=formula, family=family, yname=yname, xnames=xnames, effort.name=effort.name, number.name=number.name, totalRuns=totalRuns, totalFolds=totalFolds, remove.at.fold.flag=remove.at.fold.flag, remove.flag=remove.flag, sampling=sampling, cutoffs=cutoffs, cooksd.flag=cooksd.flag, cookthreshold=cookthreshold, weight.flag=weight.flag, weight.name=weight.name, parallel=parallel, trace=trace)
	}

	names(resList.RUNS) <- paste("RUN_", 1:totalRuns, sep="")

	subunlist <- function(theList, type.name="normal") {
		all.theList <- sapply(sapply(theList, '['), '[', type.name)

		RECALLs <- PRECISIONs <- F1s <- ACCURACYs <- AUCs <- CEs <- ERs <- NULL
		for (i in seq(all.theList)) {
			RECALLs <- rbind(RECALLs, all.theList[[i]][["RECALL"]])
			PRECISIONs <- rbind(PRECISIONs, all.theList[[i]][["PRECISION"]])
			ACCURACYs <- rbind(ACCURACYs, all.theList[[i]][["ACCURACY"]])
			F1s <- rbind(F1s, all.theList[[i]][["F1"]])
			CEs <- rbind(CEs, all.theList[[i]][["CE"]])
			ERs <- rbind(ERs, all.theList[[i]][["ER"]])
			AUCs <- rbind(AUCs, all.theList[[i]][["AUC"]])
		}

		return(list(AUC=AUCs, ER=ERs, PRECISION=PRECISIONs, RECALL=RECALLs, F1=F1s, ACCURACY=ACCURACYs, CE=CEs))
	}

	normal.effort <- subunlist(theList=resList.RUNS, type.name="normal.effort")
	divide.effort <- subunlist(theList=resList.RUNS, type.name="divide.effort")

	normal.non_effort <- subunlist(theList=resList.RUNS, type.name="normal.non_effort")
	divide.non_effort <- subunlist(theList=resList.RUNS, type.name="divide.non_effort")

	return(list(normal.effort=normal.effort, divide.effort=divide.effort, normal.non_effort=normal.non_effort, divide.non_effort=divide.non_effort))
}
