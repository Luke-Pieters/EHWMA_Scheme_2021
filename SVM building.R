# Biulding SVM 
# Helper packages
library(ggplot2)
library(cowplot)
# Modeling packages
library(caret)
library(kernlab)

#data
library(readxl)
Full_data <- read_excel("SVM Training Data p=2.xlsx")
#View(Full_data)

summary(Full_data)
Full_data$Type <- as.factor(Full_data$Type)
Full_data <- Full_data[,c(1,2,4)]

set.seed(2021)
k <- sample(nrow(Full_data),nrow(Full_data)*0.7)
Training_data <- Full_data[k,]
Test_data <- Full_data[-k,]
Test_true_type <- Test_data$Type
Test_data <- Test_data[,-ncol(Test_data)]

#Test_data <- Test_data[,-3]

#fit SMV
#model 1

Rad_SMV <- train(Type ~ .,
                data = Full_data,
                method = "svmRadial",               
                trControl = trainControl(method = "cv", number = 5),
                tuneLength = 10)

Lin_SMV <- train(Type ~ .,
                 data = Full_data,
                 method = "svmLinear",               
                 trControl = trainControl(method = "cv", number = 5),
                 tuneLength = 10)


Rad_SMV$results
Rad_SMV$bestTune
Rad_SMV$finalModel
Rad_SMV$terms

Rad_SMV$control$indexOut

Lin_SMV$results
Lin_SMV$bestTune
Lin_SMV$finalModel
Lin_SMV$terms

predicted_type <- predict(object = Rad_SMV,Test_data)
confusionMatrix(predicted_type,Test_true_type)

predicted_type <- predict(object = Lin_SMV,Test_data)
confusionMatrix(predicted_type,Test_true_type)



#plot resulting models


#make grid
xgrid <- vector()
xmin <- c(min(Full_data$V1),min(Full_data$V2))
xmax <- c(max(Full_data$V1),max(Full_data$V2))

x1seq <- seq(xmin[1]-0.02,xmax[1]+0.02,0.005)
x2seq <- seq(xmin[2]-0.05,xmax[2]+0.05,0.02)

xgrid <- cbind(rep(x1seq,each=length(x2seq)),rep(x2seq,length(x1seq)))

colours <- c("1"="chartreuse3","2"="deepskyblue")

pnts <- c(28.05685,45.18514,"2")
pnts <- rbind(pnts,c(28.05685,46.00514,"1"))
pnts <- as.data.frame(pnts)
names(pnts) <- c("V1","V2","Type")
pnts$Type <- as.factor(pnts$Type)
pnts$V1 <- as.numeric(pnts$V1)
pnts$V2 <- as.numeric(pnts$V2)

#make plot for Radial results
ygrid <- as.numeric(predict(Rad_SMV,xgrid))
plot_rad <- ggplot(mapping = aes(x=xgrid[,1],y=xgrid[,2]))+
  geom_point(aes(x=pnts$V1,y=pnts$V2,col=pnts$Type),size=2,alpha=0.9)+
  geom_point(col=ygrid+2,pch=16,size=2) +
  scale_color_manual(values = colours)+
  labs(title = "Prediction Regions using Radial SVM",colour = "Prediction")+
  xlab("Charateristic 1")+
  ylab("Charateristic 2")+
  theme_minimal()

plot_rad

#make plot for Linear results
ygrid <- as.numeric(predict(Lin_SMV,xgrid))
plot_lin <- ggplot(mapping = aes(x=xgrid[,1],y=xgrid[,2]))+
  geom_point(col=ygrid+2,pch=16,size=2) +
  geom_point(aes(x=Full_data$V1,y=Full_data$V2,col=Full_data$Type),
             size=2,alpha=0.9)+
  scale_color_manual(values = colours)+
  labs(title = "Prediction Regions using Linear SVM",
       colour = "True Classififcation")+
  xlab("Charateristic 1")+
  ylab("Charateristic 2")+
  theme_minimal()

#print plot
ggdraw()+
  draw_plot(plot_rad, x = 0.05, y = 0.2, width = .4, height = .5) +
  draw_plot(plot_lin, x = 0.45, y = 0.2, width = .5, height = .5)



