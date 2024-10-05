# 读取CSV文件
data <- read.csv("1and2-(ASE20Li-ICSE20Li)AUC.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

# 使用原始数据绘制箱型图
boxplot(data$manualDown, data$DTB_RF_OptADPT, data$BiLO_CPDP,
        names=c("ManualDown", "DTB-RF-OptADPT", "BiLO-CPDP"),
        main="Comparison of Techniques" , ylab="AUC")

# 添加图例说明，并调整位置以避免遮挡图像
legend( legend=c("ManualDown", "DTB-RF-OptADPT", "BiLO-CPDP"), 
       fill=c("lightblue", "lightgreen", "pink"),
       "bottomright",cex=0.565) # cex调整图例大小