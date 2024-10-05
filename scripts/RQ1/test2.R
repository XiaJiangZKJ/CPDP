library(effsize)
# 读取CSV文件
data <- read.csv("6-(ESA24Nikravesh)AUCandF1.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

pvalue_MDvsDTB = wilcox.test(data$manualDown_AUC, data$RF_RS__AUC, paired = TRUE)
effsize_MDvsDTB = cliff.delta(data$manualDown_AUC, data$RF_RS__AUC)
show(pvalue_MDvsDTB)
show(effsize_MDvsDTB) 

pvalue_MDvsBiLO = wilcox.test(data$manualDown_F1,data$RF_DEPT_M2__F1, paired = TRUE)
effsize_MDvsBiLO = cliff.delta(data$manualDown_F1,data$RF_DEPT_M2__F1)
show(pvalue_MDvsBiLO) 
show(effsize_MDvsBiLO)
