# 读取CSV文件
data <- read.csv("6-(ESA24Nikravesh)AUCandF1.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

# 使用原始数据绘制箱型图
boxplot(data$manualDown_F1, data$RF_DEPT_M2__F1,
        names=c("ManualDown", "DEPT-RF"),
        main="Comparison of Techniques" , ylab="F1")

# 添加图例说明，并调整位置以避免遮挡图像
legend( legend=c("ManualDown", "DEPT-RF"), 
        fill=c("lightblue", "lightgreen", "pink"),
        "bottomright",cex=0.565) # cex调整图例大小