library(effsize)
# 读取CSV文件
data <- read.csv("4-(TR24Abdu)AUCandF1.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

pvalue_MDvsDTB = wilcox.test(data$manualDown_AUC, data$DH_CNN_AUC, paired = TRUE)
effsize_MDvsDTB = cliff.delta(data$manualDown_AUC, data$DH_CNN_AUC)
show(pvalue_MDvsDTB)
show(effsize_MDvsDTB) 

pvalue_MDvsBiLO = wilcox.test(data$manualDown_F1,data$DH_CNN_F1, paired = TRUE)
effsize_MDvsBiLO = cliff.delta(data$manualDown_F1,data$DH_CNN_F1)
show(pvalue_MDvsBiLO) 
show(effsize_MDvsBiLO)
