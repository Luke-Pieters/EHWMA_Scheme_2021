#SVM for real data

# Helper packages
library(ggplot2)
# Modeling packages
library(caret)
library(kernlab)

#data
library(readxl)
Real_data <- read_excel("Real data.xlsx", 
                        sheet = "Data", col_names = FALSE)
View(Real_data)

n <- 5
p <- 2
N <- nrow(Real_data)

phi1 <- 0.25
phi2 <- 0.05
h <- 10.30

mu0 <- c(28.29, 45.85) #true mean of process
sigma0 <- matrix(c(0.0035,-0.0046,
                   -0.0046, 0.0226),
                 nrow = 2,ncol = 2,byrow = T)#true var of process

#test if OOC function
test.ooc <- function(x,cl){
  if ((x>=cl)) {
    TRUE
  }else{
    FALSE
  }#ifelse
}#test function

x1_bar <- rowMeans(Real_data[,1:n])#average of each sample for X1
x2_bar <- rowMeans(Real_data[,(n+1):10])#average of each sample for X2

Yt_1 <- mu0
Yt_bar <- mu0
t <- 1
signal <- FALSE
OOC <- vector()
T2_vec <- vector()
Y_sum <- 0
for (i in 1:N) {
  Yt <- c(x1_bar[i],x2_bar[i])
  
  if (t==1) {
    var_eh <- (phi1**2)*sigma0/n
  }else{
    var_eh <- (phi1**2) + ((1-phi1-(t-2)*phi2)/(t-1))**2 + 
      (((1-phi1+phi2)/(t-1))**2)*(t-2)
    var_eh <- var_eh*(sigma0)/n
  }#ifelse t=1
  
  EH <- phi1*Yt -phi2*Yt_1 +(1-phi1+phi2)*Yt_bar
  T2 <- t(EH-mu0)%*%solve(var_eh)%*%(EH-mu0)
  
  T2_vec <- append(T2_vec,T2)
  
  signal <- test.ooc(T2,h)
  
  if (signal) {
    OOC <-append(OOC,TRUE)
    Yt_1 <- mu0
    Yt_bar <- mu0
    t <- 1
    signal <- FALSE
    Y_sum <- 0
  }else{
    OOC <- append(OOC,FALSE)
    Y_sum <- Y_sum+Yt
    Yt_bar <- Y_sum/t
    t= t+1
    Yt_1 <- Yt
    
  }#else
}#for

#plot charting T2 stats for each sample
colours <- c("TRUE"="red","FALSE"="blue")
Chart_data <- data.frame(as.factor(c(1:N)),T2_vec,OOC)
ggplot(Chart_data,aes(Chart_data[,1],Chart_data[,2],col=Chart_data[,3],group=1))+
  geom_path()+
  geom_point()+
  geom_hline(yintercept = h,col="red",linetype="longdash")+
  scale_color_manual(values = colours)+
  labs(title = "T2 values for each sample",colour = "OOC")+
  xlab("Sample")+
  ylab("T2")

#build table for output and SVM predictions
OOC_table <- data.frame(x1_bar,x2_bar,T2_vec,OOC)
names(OOC_table) <- c("V1","V2","T2","OOC")
SVM_data <- OOC_table[OOC_table$OOC,]#take only the OOC samples
names(SVM_data) <- c("V1","V2","T2","OOC")

OOC_var_Prediction <- predict(Rad_SMV,SVM_data)#make classification prediction

SVM_data <- cbind(SVM_data,OOC_var_Prediction)

#plot prediction regions and sample predictions
colours <- c("1"="green4","2"="blue4")
ygrid <- as.numeric(predict(Rad_SMV,xgrid))
ggplot(mapping = aes(x=xgrid[,1],y=xgrid[,2]))+
  geom_point(col=ygrid+2,pch=16,size=2) +
  geom_point(aes(x=SVM_data$V1,y=SVM_data$V2,
                 col=SVM_data$OOC_var_Prediction),size=2,alpha=0.9)+
  scale_color_manual(values = colours)+
  labs(title = "Classification using Radial SVM",colour = "Classififcation")+
  xlab("Charateristic 1")+
  ylab("Charateristic 2")+
  theme_minimal()
