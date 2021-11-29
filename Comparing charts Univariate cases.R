#comparing charts univariate cases

library(xlsx) #write to excel
library(dplyr)
set.seed(1234)

#Set up
  sims <- 10000
  n <- 1 #sample size
  phi <- c(0.1,0.25,0.5)
  
  phi2 <- matrix(c(0.01,0.05,0.09,
                   0.05,0.1,0.2,
                   0.05,0.1,0.25),
                 byrow = T,nrow = 4,ncol = 3)#3 levels of l2 for each l1
  
  
  L_hw <- c(2.514,2.767,2.807)
  L_ew <- c(2.483,2.687,2.778)
  L_eew <- matrix(c(2.399,2.105,2.310,
                2.523,2.395,2.527,
                2.705,2.651,2.625),byrow = T,nrow = 3,ncol = 3)
  
  
  mu0 <- 0
  sigma0 <- 1 #sdv
  
  shifts <- seq(0,3,0.25)
  
  #test if OOC function
  test.ooc <- function(x,ucl,lcl){
    if ((x>=ucl)|(x<=lcl)) {
      TRUE
    }else{
      FALSE
    }#ifelse
  }#test function
  
  results_EW <- vector()
  names_EW <- vector()
  results_EEW <- vector()
  names_EEW <- vector()
  results_HW <- vector()
  names_HW <- vector()
  
  #progress bar
  total <- length(phi)*length(shifts)*3
  progress_counter <- 0
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = total, width = 300)
  
#start sims: EWMA chart======
#============================
  for (k in 1:length(shifts)) {
    mean_shifted <- mu0 + shifts[k]*sigma0
    vec <- vector()
    for (i in 1:length(phi)) {
        count_vec <- vector()
        progress_counter <- progress_counter +1 #progress bar counter
        setWinProgressBar(pb, progress_counter,
                            label = paste( round(progress_counter/total*100, 0),
                                       "% done"))
        for (c in 1:sims) {
          t <- 0
          signal <- FALSE
          Xt <- vector()
          Zt_1 <- mu0
          
          while ((t<100000)&(signal==FALSE)) {
            t <- t +1
            y <- rnorm(n,mean = mean_shifted,sd = sigma0)
            
            var_z <- (phi[i]/(2-phi[i]))*(1-((1-phi[i])**(2*t)))*(sigma0**2)/n
            
            ucl <- mu0 + L_ew[i]*sqrt(var_z)
            lcl <- mu0 - L_ew[i]*sqrt(var_z)
            
            Xt[t] <- mean(y)
            
            Zt <- phi[i]*Xt[t] +(1-phi[i])*Zt_1
            signal <- test.ooc(Zt,ucl,lcl)
            Zt_1 <- Zt
            
          }#while IC
          count_vec[c] <- t
        }#for sims
        vec <- append(vec,c(mean(count_vec),sd(count_vec),median(count_vec)))
    }#for phi 
    results_EW <- rbind(results_EW,t(vec))
  }#for shifts
  
  
#start sims: HWMA chart======
#============================
  for (k in 1:length(shifts)) {
    mean_shifted <- mu0 + shifts[k]*sigma0
    vec <- vector()
    for (i in 1:length(phi)) {
      count_vec <- vector()
      progress_counter <- progress_counter +1
      setWinProgressBar(pb, progress_counter,
                        title=paste( round(progress_counter/total*100, 0),
                                     "% done"))
      for (c in 1:sims) {
        t <- 0
        signal <- FALSE
        Xt <- vector()
        Xt_bar <- mu0
        
        while ((t<100000)&(signal==FALSE)) {
          t <- t +1
          y <- rnorm(n,mean = mean_shifted,sd = sigma0)
          
          if (t==1) {
            var_h <- (phi[i]**2)*(sigma0**2)/n
          }else{
            var_h <- (phi[i]**2) + (((1-phi[i])**2)/(t-1))
            var_h <- var_h*(sigma0**2)/n
          }#ifelse t=1
          
          ucl <- mu0 + L_hw[i]*sqrt(var_h)
          lcl <- mu0 - L_hw[i]*sqrt(var_h)
          
          Xt[t] <- mean(y)
          
          Ht <- phi[i]*Xt[t] +(1-phi[i])*Xt_bar
          signal <- test.ooc(Ht,ucl,lcl)
          Xt_bar <- mean(Xt)
          
        }#while IC
        count_vec[c] <- t
      }#for sims
      vec <- append(vec,c(mean(count_vec),sd(count_vec),median(count_vec)))
    }#for phi 
    results_HW <- rbind(results_HW,t(vec))
  }#for shifts 

  
  #start sims: Extended EWMA chart======
  #============================
  for (k in 1:length(shifts)) {
    mean_shifted <- mu0 + shifts[k]*sigma0
    vec <- vector()
    for (i in 1:length(phi)) {
      progress_counter <- progress_counter +1
      setWinProgressBar(pb, progress_counter,
                        label = paste( round(progress_counter/total*100, 0),
                                       "% done"))
      for (j in 1:ncol(phi2)) {
        count_vec <- vector()
        for (c in 1:sims) {
          t <- 0
          signal <- FALSE
          Xt <- vector()
          EZt_1 <- mu0
          Xt_1 <-  mu0
          
          while ((t<100000)&(signal==FALSE)) {
            t <- t +1
            y <- rnorm(n,mean = mean_shifted,sd = sigma0)
            
            var_ez <- (phi[i]**2 + 
                      phi2[i,j]**2)*((1-
                      (1-phi[i]+phi2[i,j])^(2*t))/(2*(phi[i]-phi2[i,j]) -
                                                     (phi[i]-phi2[i,j])^2))
            var_ez <- var_ez -
                      2*(1-phi[i] +
                      phi2[i,j])*(phi[i]*phi2[i,j])*((1-
                      (1-phi[i]+phi2[i,j])^(2*(t-1)))/(2*(phi[i]-phi2[i,j]) - 
                                                         (phi[i]-phi2[i,j])^2))
            var_ez <- var_ez*(sigma0**2)/n
            
            ucl <- mu0 + L_eew[i,j]*sqrt(var_ez)
            lcl <- mu0 - L_eew[i,j]*sqrt(var_ez)
            
            Xt[t] <- mean(y)
            
            EZt <- phi[i]*Xt[t] - phi2[i,j]*Xt_1 +(1-phi[i])*EZt_1
            signal <- test.ooc(EZt,ucl,lcl)
            EZt_1 <- EZt
            Xt_1 <- Xt[t]
            
          }#while IC
          count_vec[c] <- t
        }#for sims
        vec <- append(vec,c(mean(count_vec),sd(count_vec),median(count_vec)))
      }#for phi2
    }#for phi 
    results_EEW <- rbind(results_EEW,t(vec))
  }#for shifts  

  
#prep output
  for (i in 1:length(phi)) {
    newnames1 <- cbind(rep(phi[i],3),rep(L_ew[i],3))
    names_EW <- cbind(names_EW,t(newnames1))
    
    newnames2 <- cbind(rep(phi[i],3),rep(L_hw[i],3))
    names_HW <- cbind(names_HW,t(newnames2))
    
    for (j in 1:ncol(phi2)) {
      newnames3 <- cbind(rep(phi[i],3),rep(phi2[i,j],3),rep(L_eew[i,j],3))
      names_EEW <- cbind(names_EEW,t(newnames3))
    }
  }
  
  results_EW <- as.data.frame(rbind(names_EW,
                                    t(rep(c("ARL","SDRL","MRL"),3)),
                                    results_EW))
  rownames(results_EW) <- c("Phi","L","Measure",shifts)
  
  results_HW <- as.data.frame(rbind(names_HW,
                                    t(rep(c("ARL","SDRL","MRL"),3)),
                                    results_HW))
  rownames(results_HW) <- c("Phi","L","Measure",shifts)
  
  results_EEW <- as.data.frame(rbind(names_EEW,
                                     t(rep(c("ARL","SDRL","MRL"),3*3)),
                                     results_EEW))
  rownames(results_EEW) <- c("Phi","Phi2","L","Measure",shifts)  
  
  file_name <- "Comparing Charts Univariate case.xlsx"
  #print results to excel
  write.xlsx(results_EW,
             sheetName = "EWMA",
             file = file_name,row.names=TRUE,col.names = FALSE)
  write.xlsx(results_HW,
             sheetName = "HWMA",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(results_EEW,
             sheetName = "Extended EWMA",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  close(pb)
  