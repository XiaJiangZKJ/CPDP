library(effsize)
# 读取CSV文件
data <- read.csv("3-(TSE23Li)AUCandMCC.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

pvalue_MDvsDTB = wilcox.test(data$manualDown_AUC, data$DSSDPP_AUC, paired = TRUE)
effsize_MDvsDTB = cliff.delta(data$manualDown_AUC, data$DSSDPP_AUC)
show(pvalue_MDvsDTB)
show(effsize_MDvsDTB) 

pvalue_MDvsBiLO = wilcox.test(data$manualDown_MCC,data$DSSDPP_MCC, paired = TRUE)
effsize_MDvsBiLO = cliff.delta(data$manualDown_MCC,data$DSSDPP_MCC)
show(pvalue_MDvsBiLO) 
show(effsize_MDvsBiLO)
