require(RWeka)

source("util.R")

root <- "(final)data"

method = "DOWN"

files <- list.files(path = root, full.names = TRUE)
files <- files[!R.utils::isDirectory(files)]

all.F1 <- all.MCC <- all.AUC <- NULL

df.stat <- NULL
for (file in files) {
  cat(file, "\n")
  
  #?????ļ??ĵ?һ????SLOC??????һ????bug??Ϣ
  data <- read.csv(file, sep = ",", stringsAsFactors = FALSE)
  
  data$NUM <- as.numeric(data[, ncol(data)])
  data$LOC <- as.numeric(data[, 1])
  data$BUG <- 1 * (as.numeric(data$NUM >= 1))
  
  data[is.na(data)] <- 0
  
  
  df.stat <-
    rbind(df.stat, c(
      nrow(data),
      sum(data$BUG),
      sum(data$BUG) / nrow(data),
      sum(data$NUM)
    ))
  
  cname <- "LOC"
  
  valid.dt <-
    data.frame(
      NUM = data$NUM,
      REL = data$BUG,
      LOC = data$LOC,
      PRE = data[, cname]
    )
  
  sorted           <- FALSE
  worstcase        <- TRUE
  bestcase         <- FALSE
  LOCUP            <- FALSE
  allpercentcutoff <- FALSE
  
  
  AUC <-
    ComputeAUC(
      sorted = sorted,
      worstcase = worstcase,
      bestcase = bestcase,
      LOCUP = LOCUP,
      allpercentcutoff = allpercentcutoff,
      data = valid.dt
    )
  cp <-
    ComputeCPs(
      sorted = sorted,
      worstcase = worstcase,
      bestcase = bestcase,
      LOCUP = LOCUP,
      allpercentcutoff = allpercentcutoff,
      train.dt = valid.dt,
      valid.dt = valid.dt
    )
  
  
  aF1 <- cp[["F1"]][paste("F1.P.", seq(0.1, 1, 0.1), sep = "")]
  
  aMCC <- cp[["MCC"]][paste("MCC.P.", seq(0.1, 1, 0.1), sep = "")]
  aAUC <- AUC
  
  
  all.F1 <- rbind(all.F1, aF1)
  all.MCC <- rbind(all.MCC, aMCC)
  all.AUC <- rbind(all.AUC, aAUC)
}

path_res <- file.path(root, "result")
if (!dir.exists(path_res)) {
  dir.create(path_res)
}

colnames(df.stat) <- c("N", "#BUGGY", "#buggyrate", "#bugs")
rownames(df.stat) <- basename(files)
write.csv(df.stat, file = file.path(path_res, "summary-all.csv"))


rownames(all.F1) <-
  rownames(all.MCC) <- rownames(all.AUC) <- basename(files)

### Classification
write.csv(round(all.F1, 3),
          file = file.path(path_res, paste("all-F1-", method, ".csv", sep = "")),
          row.names = TRUE)

write.csv(round(all.MCC, 3),
          file = file.path(path_res, paste("all-MCC-", method, ".csv", sep = "")),
          row.names = TRUE)
write.csv(round(all.AUC, 3),
          file = file.path(path_res, paste("all-AUC-", method, ".csv", sep = "")),
          row.names = TRUE)

 
