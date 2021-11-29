#UNIVARIATE CASES
library(xlsx) #write to excel
set.seed(1234)

#Set up
  sims <- 10000
  n <- 1 #sample size
  phi1 <- c(0.1,0.25,0.5,0.9)
  phi2 <- matrix(c(0.01,0.05,0.09,
                      0.05,0.1,0.2,
                      0.05,0.1,0.25,
                      0.05,0.1,0.25),
                 byrow = T,nrow = 4,ncol = 3)#3 levels of l2 for each l1
  
  mu0 <- 0
  sigma0 <- 1
  L <- matrix(c(2.516,2.540,2.600,
                2.772,2.763,2.762,
                2.804,2.803,2.794,
                2.804,2.809,2.801),
              byrow = T,nrow = 4,ncol = 3)#L for each phi2 phi1 pair
  
  shifts <- seq(0,3,0.25)
  
  #test if OOC function
  test.ooc <- function(x,ucl,lcl){
    if ((x>=ucl)|(x<=lcl)) {
      TRUE
    }else{
      FALSE
    }#ifelse
  }#test function

  results_ARL <- vector() #avg
  results_SDRL <- vector() #stdv
  results_MRL <- vector() #median
  
  total <- ncol(phi2)*length(phi1)*length(shifts)
  progress_counter <- 0
  pb <- winProgressBar(title = "progress bar",label = "%", min = 0,
                       max = total, width = 300)
  
#Start Sims
    for (k in 1:length(shifts)) {
      mean_shifted <- mu0 + shifts[k]*sigma0
      ARL_vec <- vector()
      SDRL_vec <- vector()
      MRL_vec <- vector()
      for (i in 1:length(phi1)) {
        for (j in 1:ncol(phi2)) {
          count_vec <- vector()
          progress_counter <- progress_counter +1
          setWinProgressBar(pb, progress_counter,
                            label =paste( round(progress_counter/total*100, 0),
                                          "% done"))
          for (c in 1:sims) {
            t <- 0
            signal <- FALSE
            Xt <- vector()
            Xt_1 <- mu0
            Xt_bar <- mu0
            
            while ((t<100000)&(signal==FALSE)) {
              t <- t +1
              y <- rnorm(n,mean = mean_shifted,sd = sigma0)
              
              if (t==1) {
                var_eh <- (phi1[i]**2)*(sigma0**2)/n
              }else{
                var_eh <- (phi1[i]**2) + 
                          ((1-phi1[i]-(t-2)*phi2[i,j])/(t-1))**2 + 
                          (((1-phi1[i]+phi2[i,j])/(t-1))**2)*(t-2)
                var_eh <- var_eh*(sigma0**2)/n
              }#ifelse t=1
              
              ucl <- mu0 + L[i,j]*sqrt(var_eh)
              lcl <- mu0 - L[i,j]*sqrt(var_eh)
              
              Xt[t] <- mean(y)
              
              EH <- phi1[i]*Xt[t] -
                    phi2[i,j]*Xt_1 +(1-phi1[i]+phi2[i,j])*Xt_bar
              signal <- test.ooc(EH,ucl,lcl)
              Xt_1 <- Xt[t]
              Xt_bar <- mean(Xt)
              
            }#while IC
            count_vec[c] <- t
          }#for sims
          ARL_vec <- append(ARL_vec,mean(count_vec))
          SDRL_vec <- append(SDRL_vec,sd(count_vec))
          MRL_vec <- append(MRL_vec,median(count_vec))
        }#for phi 2
      }#for phi 1
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
      for (j in 1:ncol(phi2)) {
        newrows <- cbind(newrows,c(phi1[i],phi2[i,j]))
      }
    }
    newrows <- rbind(newrows,as.vector(t(L)))
    
    results_ARL <- rbind(newrows,results_ARL)
    results_SDRL <- rbind(newrows,results_SDRL)
    results_MRL <- rbind(newrows,results_MRL)
    
    rownames(results_ARL) <- c("phi1","phi2","L",shifts)
    rownames(results_SDRL) <- c("phi1","phi2","L",shifts)
    rownames(results_MRL) <- c("phi1","phi2","L",shifts)
    
    file_name <- paste("Univariate case n=",n,".xlsx",sep = "")
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
    
    close(pb)