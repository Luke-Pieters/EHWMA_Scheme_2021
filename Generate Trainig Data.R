#Generating Training data 
library(xlsx) #write to excel
set.seed(1234)

phi1 <- 0.25
phi2 <- 0.05
h <- 10.34

p <- 2

datapoints <- 30 #how many samples for each shift size

shifts <- seq(0.25,3,0.25)
n <- 5

counter <- 0
type <- 0

mu0 <- c(28.29, 45.85)#true mean of process
sigma0 <- matrix(c(0.0035,-0.0046,
                   -0.0046, 0.0226),
                 nrow = 2,ncol = 2,byrow = T)#true var of process

output0 <- vector()

#test if OOC function
test.ooc <- function(x,cl){
  if ((x>=cl)) {
    TRUE
  }else{
    FALSE
  }#ifelse
}#test function

total <- p*length(shifts)*datapoints 
progress_counter <- 0
pb <- winProgressBar(title = "progress bar",label = "% done", min = 0,
                     max = total, width = 300)

#make IC data points
# while (counter<datapoints) {
#   t <- 0
#   signal <- FALSE
#   Yt <- vector()
#   Yt_1 <- mu0
#   Yt_bar <- mu0
#   while ((t<100000)&(signal==FALSE)&(counter<datapoints)) {
#     t <- t+1
#     z <- MASS::mvrnorm(n,mu = mu0,Sigma = sigma0)
# 
#     if (t==1) {
#       var_eh <- (phi1**2)*sigma0/n
#     }else{
#       var_eh <- (phi1**2) + ((1-phi1-(t-2)*phi2)/(t-1))**2 +
#         (((1-phi1+phi2)/(t-1))**2)*(t-2)
#       var_eh <- var_eh*(sigma0)/n
#     }#ifelse t=1
# 
#     Yt <- cbind(Yt,colMeans(z))
# 
#     EH <- phi1*Yt[,t] -phi2*Yt_1 +(1-phi1+phi2)*Yt_bar
#     T2 <- t(EH-mu0)%*%solve(var_eh)%*%(EH-mu0)
#     signal <- test.ooc(T2,h)
#     Yt_1 <- Yt[,t]
#     Yt_bar <- rowMeans(Yt)
# 
#     if (signal==F ) {
#       output0 <- rbind(output0,t(c(Yt[,t],T2,type)))
#       counter <- nrow(output0)
#     }
#   }#while IC
# }#while

#final_output <- output0

final_output <- vector()
output0 <- vector()

#for OOC types

for (d in 1:p) {
  type <- d
  for (k in 1:length(shifts)) {
    delta <- rep(0,p)
    delta[d] <- shifts[k]*sqrt(sigma0[d,d])
    signchange <- -1
    
    for (i in 1:datapoints) {
      progress_counter <- progress_counter +1
      setWinProgressBar(pb, progress_counter,
                        label =paste( round(progress_counter/total*100, 0),
                                      "% done"))
      signchange <- -1*signchange
      mean_shifted <- mu0 + signchange*delta 
      t <- 0
      signal <- FALSE
      Yt <- vector()
      Yt_1 <- mu0
      Yt_bar <- mu0
      while ((t<100000)&(signal==FALSE)) {
        t <- t+1
        z <- MASS::mvrnorm(n,mu = mean_shifted,Sigma = sigma0)
        z_bar <- colMeans(z)
        
        if (t==1) {
          var_eh <- (phi1**2)*sigma0/n
        }else{
          var_eh <- (phi1**2) + ((1-phi1-(t-2)*phi2)/(t-1))**2 + 
            (((1-phi1+phi2)/(t-1))**2)*(t-2)
          var_eh <- var_eh*(sigma0)/n
        }#ifelse t=1
        
        Yt <- cbind(Yt,z_bar)
        
        EH <- phi1*Yt[,t] -phi2*Yt_1 +(1-phi1+phi2)*Yt_bar
        T2 <- t(EH-mu0)%*%solve(var_eh)%*%(EH-mu0)
        signal <- test.ooc(T2,h)
        Yt_1 <- Yt[,t]
        Yt_bar <- rowMeans(Yt)
      }#while IC  
      newrow <- c(z_bar,T2,type)
      output0 <- rbind(output0,t(newrow))
    }#for datapoints
  }#for shifts
}#for p


final_output <- as.data.frame(rbind(final_output,output0))
colnames(final_output) <- c("V1","V2","T2","Type")

file_name <- "SVM Training Data p=2.xlsx"
#print results to excel
write.xlsx(final_output,
           sheetName = "data",
           file = file_name,row.names=FALSE,col.names = TRUE)
close(pb)

summary(final_output)

#make 3d plot and 2d plot of generated data
library(car)
library(rgl)
 
scatter3d(x=final_output$V1,
          y=final_output$V2,
          z=final_output$T2,
          groups = as.factor(final_output$Type),
           surface = FALSE)
 
plot(x=final_output$V1,y=final_output$V2,col=as.factor(final_output$Type))
