library(effsize)
# 读取CSV文件
data <- read.csv("1and2-(ASE20Li-ICSE20Li)AUC.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

pvalue_MDvsDTB = wilcox.test(data$manualDown, data$DTB_RF_OptADPT, paired = TRUE)
effsize_MDvsDTB = cliff.delta(data$manualDown, data$DTB_RF_OptADPT)
show(pvalue_MDvsDTB)
show(effsize_MDvsDTB) 

pvalue_MDvsBiLO = wilcox.test(data$manualDown,data$BiLO_CPDP, paired = TRUE)
effsize_MDvsBiLO = cliff.delta(data$manualDown,data$BiLO_CPDP)
show(pvalue_MDvsBiLO) 
show(effsize_MDvsBiLO)
