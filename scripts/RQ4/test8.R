library(effsize)
# 读取CSV文件
data <- read.csv("8.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

pvalue_MDvsDTB = wilcox.test(data$ManualDown_F1, data$DP_CLME_F1, paired = TRUE)
effsize_MDvsDTB = cliff.delta(data$ManualDown_F1, data$DP_CLME_F1)
show(pvalue_MDvsDTB)
show(effsize_MDvsDTB) 
