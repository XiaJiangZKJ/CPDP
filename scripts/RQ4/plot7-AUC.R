# 读取CSV文件
data <- read.csv("7.csv", header = TRUE)

# 检查数据框的前几行，确保数据正确导入
head(data)

# 替换列名中的点为下划线
colnames(data) <- gsub("\\.", "_", colnames(data))

# 使用原始数据绘制箱型图
boxplot(data$ManualDown_AUC, data$SSE_AUC,
        names=c("ManualDown", "SSE"),
        main="Comparison of Techniques" , ylab="AUC")

# 添加图例说明，并调整位置以避免遮挡图像
legend( legend=c("ManualDown", "SSE"), 
        fill=c("lightblue", "lightgreen", "pink"),
        "bottomright",cex=0.565) # cex调整图例大小

# 使用原始数据绘制箱型图
boxplot(data$ManualDown_F1, data$SSE_F1,
        names=c("ManualDown", "SSE"),
        main="Comparison of Techniques" , ylab="F1")

# 添加图例说明，并调整位置以避免遮挡图像
legend( legend=c("ManualDown", "SSE"), 
        fill=c("lightblue", "lightgreen", "pink"),
        "bottomright",cex=0.565) # cex调整图例大小

# 使用原始数据绘制箱型图
boxplot(data$ManualDown_MCC, data$SSE_MCC,
        names=c("ManualDown", "SSE"),
        main="Comparison of Techniques" , ylab="MCC")

# 添加图例说明，并调整位置以避免遮挡图像
legend( legend=c("ManualDown", "SSE"), 
        fill=c("lightblue", "lightgreen", "pink"),
        "bottomright",cex=0.565) # cex调整图例大小
