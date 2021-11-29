#Robustness UNIVARIATE CASES
library(xlsx) #write to excel
library(dplyr)
set.seed(1234)

#Set up
  sims <- 5000
  n <- 1 #sample size
  phi1 <- c(0.1,0.5)
  phi2 <- matrix(c(0.05,0.09,
                   0.05,0.25),byrow = T,nrow = 2,ncol = 2)
  
  L <- matrix(c(2.543,2.6,
               2.81,2.804),byrow = T,nrow = 2,ncol = 2)
  
  percentiles <- c(0.05,0.25,0.5,0.75,0.95)
  
  dists <- c("N(0,1)","t(10)","t(100)","t(1000)",
             "GAM(1,1)","GAM(10,1)","LogNorm(0,1)","X2(30)")
  
  means <- c(0,0,0,0,1,10,exp(0.5),30)
  vars <- c(1,1.25,50/49,1000/998,1,10,exp(1)*(exp(1)-1),60)
  
  output <- vector()
  
  #test if OOC function
  test.ooc <- function(x,ucl,lcl){
    if ((x>=ucl)|(x<=lcl)) {
      TRUE
    }else{
      FALSE
    }#ifelse
  }#test function
  
  #random sample functions
  rdist.functions <- function(index,n){
    case_when(index == 1 ~ rnorm(n),
              index == 2 ~ rt(n,10),
              index == 3 ~ rt(n,100),
              index == 4 ~ rt(n,1000),
              index == 5 ~ rgamma(n,1,1),
              index == 6 ~ rgamma(n,10,1),
              index == 7 ~ rlnorm(n,0,1),
              index == 8 ~ rchisq(n,30),
              )
  }#rdist function
  
  
  
  
#run sims
  for (i in 1:length(phi1)) {
    for (j in 1:ncol(phi2)) {
      for (d in 1:length(dists)) {
        count_vec <- vector()
        mu0 <- means[d]
        sigma0 <- sqrt(vars[d])
        for (c in 1:sims) {
          t <- 0
          signal <- FALSE
          Xt <- vector()
          Xt_1 <- mu0
          Xt_bar <- mu0
          while ((t<100000)&(signal==FALSE)) {
            t <- t+1
            y <- rdist.functions(d,n)
            
            if (t==1) {
              var_eh <- (phi1[i]**2)*(sigma0**2)/n
            }else{
              var_eh <- (phi1[i]**2) + ((1-phi1[i]-(t-2)*phi2[i,j])/(t-1))**2 + 
                (((1-phi1[i]+phi2[i,j])/(t-1))**2)*(t-2)
              var_eh <- var_eh*(sigma0**2)/n
            }#ifelse t=1
            
            ucl <- mu0 + L[i,j]*sqrt(var_eh)
            lcl <- mu0 - L[i,j]*sqrt(var_eh)
            
            Xt[t] <- mean(y)
            
            EH <- phi1[i]*Xt[t] -phi2[i,j]*Xt_1 +(1-phi1[i]+phi2[i,j])*Xt_bar
            signal <- test.ooc(EH,ucl,lcl)
            Xt_1 <- Xt[t]
            Xt_bar <- mean(Xt)
          }#while IC
          count_vec[c] <- t
        }#for sims
        newrow <- c(phi1[i],phi2[i,j],L[i,j],
                    dists[d],mean(count_vec),sd(count_vec),
                    quantile(count_vec,probs = percentiles))
        output <- cbind(output,newrow)
      }#for dists
    }#for phi2
  }#for phi1
  
  output <- as.data.frame(t(output))
  names(output) <- c("phi1","phi2","L",
                     "Dist","ARL","SDRL","P5","P25","P50","P75","P95")
  
  #print output to excel
  file_name <- "Robustness Univariate case.xlsx"
  write.xlsx(output,
             sheetName = "Univariate",
             file = file_name,row.names=FALSE,col.names = TRUE)
  
  
  