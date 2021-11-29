#MULTIVARIATE CASES
library(xlsx) #write to excel
set.seed(1234)

#Set up
  sims <- 2000
  n <- 1 #sample size
  phi1 <- c(0.1,0.25)
  phi2 <- matrix(c(0.01,0.05,0.09,
                   0.05,0.1,0.2),
                 byrow = T,nrow = 2,ncol = 3)#3 levels of l2 for each l1
  p <- c(2,3,4)
  h1 <- matrix(c(9.00,9.15,9.44,
                 11.09,11.32,11.73,
                13.10,13.29,13.62),
              byrow = T,nrow = 3,ncol = 3)
  h2 <- matrix(c(10.34,10.37,10.38,
                 12.70,12.70,12.62,
                 14.60,14.62,14.68),
               byrow = T,nrow = 3,ncol = 3)
  h <- array(c(h1,h2),dim = c(3,3,2))#p,phi2,phi1
  
  shifts <- seq(0,3,0.25)
  #shifts <- c(0)
  
  #test if OOC function
  test.ooc <- function(x,cl){
    if ((x>=cl)) {
      TRUE
    }else{
      FALSE
    }#ifelse
  }#test function
  
  results_ARL <- vector() #avg
  results_SDRL <- vector() #stdv
  results_MRL <- vector() #median
  
  
#start sims
  for (k in 1:length(shifts)) {
    ARL_vec <- vector()
    SDRL_vec <- vector()
    MRL_vec <- vector()
    for (i in 1:length(phi1)) {
      for (d in 1:length(p)) {
        mu0 <- rep(0,p[d])
        sigma0 <- diag(p[d])
        delta <- sqrt((shifts[k]**2)/p[d])
        mean_shifted <- mu0 + sigma0%*%(rep(delta,p[d]))
        for (j in 1:ncol(phi2)) {
          count_vec <- vector()
          for (c in 1:sims) {
            t <- 0
            signal <- FALSE
            Yt <- vector()
            Yt_1 <- mu0
            Yt_bar <- mu0
            while ((t<100000)&(signal==FALSE)) {
              t <- t+1
              z <- MASS::mvrnorm(n,mu = mean_shifted,Sigma = sigma0)
              
              if (t==1) {
                var_eh <- (phi1[i]**2)*sigma0
              }else{
                var_eh <- (phi1[i]**2) + 
                  ((1-phi1[i]-(t-2)*phi2[i,j])/(t-1))**2 + 
                  (((1-phi1[i]+phi2[i,j])/(t-1))**2)*(t-2)
                var_eh <- var_eh*(sigma0)
              }#ifelse t=1
              
              Yt <- cbind(Yt,z)
              
              EH <- phi1[i]*Yt[,t] -
                    phi2[i,j]*Yt_1 +(1-phi1[i]+phi2[i,j])*Yt_bar
              T2 <- t(EH-mu0)%*%solve(var_eh)%*%(EH-mu0)
              signal <- test.ooc(T2,h[d,j,i])
              Yt_1 <- Yt[,t]
              Yt_bar <- rowMeans(Yt)
              
            }#while IC
            count_vec[c] <- t
          }#for sims
          ARL_vec <- append(ARL_vec,mean(count_vec))
          SDRL_vec <- append(SDRL_vec,sd(count_vec))
          MRL_vec <- append(MRL_vec,median(count_vec))
        }#for phi2
      }#for p
    }#for phi1
    results_ARL <- cbind(results_ARL,ARL_vec)
    results_SDRL <- cbind(results_SDRL,SDRL_vec)
    results_MRL <- cbind(results_MRL,MRL_vec)
  }#for shifts

#prep results
  results_ARL <- as.data.frame(t(results_ARL))
  results_SDRL <- as.data.frame(t(results_SDRL))
  results_MRL <- as.data.frame(t(results_MRL))
  
  newrows <- vector()
  for (i in 1:length(phi1)) {
    for (d in 1:length(p)) {
      for (j in 1:ncol(phi2)) {
        newrows <- cbind(newrows,c(phi1[i],p[d],phi2[i,j],h[d,j,i]))
      }
    }
  }
  
  results_ARL <- rbind(newrows,results_ARL)
  results_SDRL <- rbind(newrows,results_SDRL)
  results_MRL <- rbind(newrows,results_MRL)
  
  rownames(results_ARL) <- c("phi1","p","phi2","h",shifts)
  rownames(results_SDRL) <- c("phi1","p","phi2","h",shifts)
  rownames(results_MRL) <- c("phi1","p","phi2","h",shifts)
  
  file_name <- paste("Mulrivariate case n=",n,".xlsx",sep = "")
  #print results to excel
  write.xlsx(results_ARL,
             sheetName = "ARL",
             file = file_name,row.names=TRUE,col.names = FALSE)
  write.xlsx(results_SDRL,
             sheetName = "SDRL",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(results_MRL,
             sheetName = "MRL",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
