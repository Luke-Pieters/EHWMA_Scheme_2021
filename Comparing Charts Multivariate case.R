#Comparing MULTIVARIATE CASES
library(xlsx) #write to excel
set.seed(1234)

#Set up
  sims <- 5000
  n <- 1 #sample size
  phi <- c(0.1,0.25)
  phi2 <- matrix(c(0.01,0.05,0.09,
                   0.05,0.1,0.2),
                 byrow = T,nrow = 2,ncol = 3)#3 levels of l2 for each l1
  
  p <- c(2,3,4)
  
  h_ew <- matrix(c(8.92,9.977,
                   11.03,12.20,
                   12.95,14.10),nrow = length(p),ncol = length(phi),byrow = T)
  
  h_hw <- matrix(c(8.96,10.5,
                   11.09,12.60,
                   13.10,14.65),nrow = length(p),ncol = length(phi),byrow = T)
  
  h1 <- matrix(c(8.852,8.749,8.959,
                 10.972,10.997,11.072,
                 12.906,13.07,12.984),
               byrow = T,nrow = 3,ncol = 3)
  h2 <- matrix(c(10.023,10.116,10.189,
                 12.243,12.348,12.165,
                 14.281,14.395,14.460),
               byrow = T,nrow = 3,ncol = 3)
  h_eew <- array(c(h1,h2),dim = c(3,3,2))#p,phi2,phi1
  
  shifts <- seq(0,3,0.25)
  
  #test if OOC function
  test.ooc <- function(x,cl){
    if ((x>=cl)) {
      TRUE
    }else{
      FALSE
    }#ifelse
  }#test function
  
  results_EW <- list()
  results_EEW <- list()
  results_HW <- list()
  
  #progress bar
  total <- length(p)*length(phi)*length(shifts)*3
  progress_counter <- 0
  pb <- winProgressBar(title = "progress bar",label = "%", min = 0,
                       max = total, width = 300)
  
#Start sims MEWMA========
#========================
  for (d in 1:length(p)) {
    output <- vector()
    newnames <- vector()
    mu0 <- rep(0,p[d])
    sigma0 <- diag(p[d])
    for (k in 1:length(shifts)) {
      newrow <- vector()
      delta <- sqrt((shifts[k]**2)/p[d])
      mean_shifted <- mu0 + sigma0%*%(rep(delta,p[d]))
      for (i in 1:length(phi)) {
        count_vec <- vector()
        progress_counter <- progress_counter +1 #progress bar counter
        setWinProgressBar(pb, progress_counter,
                          label =paste( round(progress_counter/total*100, 0),
                                        "% done"))
        for (c in 1:sims) {
          t <- 0
          signal <- FALSE
          Xt <- vector()
          Zt_1 <- mu0
          while ((t<100000)&(signal==FALSE)) {
            t <- t+1
            x <- MASS::mvrnorm(n,mu = mean_shifted,Sigma = sigma0)
            
            Xt <- cbind(Xt,x)
            
            var_z <- (phi[i]/(2-phi[i]))*(1-((1-phi[i])**(2*t)))*sigma0
            
            Zt <- phi[i]*Xt[,t] +(1-phi[i])*Zt_1
            T2 <- t(Zt-mu0)%*%solve(var_z)%*%(Zt-mu0)
            signal <- test.ooc(T2,h_ew[d,i])
            Zt_1 <- Zt
            
          }#while IC
          count_vec[c] <- t
        }#for sims
        newrow <- cbind(newrow,
                        t(c(mean(count_vec),sd(count_vec),median(count_vec))))
      }#for phi
      output <- rbind(output,newrow)
    }#for shifts
    newnames <- t(cbind(rep(phi,each=3),
                        rep(h_ew[d,],each=3),rep(c("ARL","SDRL","MRL"),2)))
    results_EW[[d]] <- as.data.frame(rbind(newnames,output))
  }#for p
  
  
  #Start sims MHWMA===============
  #===============================
  for (d in 1:length(p)) {
    output <- vector()
    newnames <- vector()
    mu0 <- rep(0,p[d])
    sigma0 <- diag(p[d])
    for (k in 1:length(shifts)) {
      newrow <- vector()
      delta <- sqrt((shifts[k]**2)/p[d])
      mean_shifted <- mu0 + sigma0%*%(rep(delta,p[d]))
      for (i in 1:length(phi)) {
        count_vec <- vector()
        progress_counter <- progress_counter +1 #progress bar counter
        setWinProgressBar(pb, progress_counter,
                          label =paste( round(progress_counter/total*100, 0),
                                       "% done"))
        for (c in 1:sims) {
          t <- 0
          signal <- FALSE
          Xt <- vector()
          Xt_bar <- mu0
          
          while ((t<100000)&(signal==FALSE)) {
            t <- t +1
            y <- MASS::mvrnorm(n,mu = mean_shifted,Sigma = sigma0)
            
            if (t==1) {
              var_h <- (phi[i]**2)*(sigma0)
            }else{
              var_h <- (phi[i]**2) + (((1-phi[i])**2)/(t-1))
              var_h <- var_h*(sigma0)
            }#ifelse t=1
            
            Xt <- cbind(Xt,y)
            
            Ht <- phi[i]*Xt[,t] +(1-phi[i])*Xt_bar
            T2 <- t(Ht-mu0)%*%solve(var_h)%*%(Ht-mu0)
            signal <- test.ooc(T2,h_hw[d,i])
            Xt_bar <- rowMeans(Xt)
            
          }#while IC
          count_vec[c] <- t
        }#for sims
        newrow <- cbind(newrow,
                        t(c(mean(count_vec),sd(count_vec),median(count_vec))))
      }#for phi
      output <- rbind(output,newrow)
    }#for shifts
    newnames <- t(cbind(rep(phi,each=3),
                        rep(h_hw[d,],each=3),rep(c("ARL","SDRL","MRL"),2)))
    results_HW[[d]] <- as.data.frame(rbind(newnames,output))
  }#for p
 
#Start sims extended MEWMA========
#=================================
  for (d in 1:length(p)) {
    output <- vector()
    newnames <- vector()
    mu0 <- rep(0,p[d])
    sigma0 <- diag(p[d])
    for (k in 1:length(shifts)) {
      newrow <- vector()
      delta <- sqrt((shifts[k]**2)/p[d])
      mean_shifted <- mu0 + sigma0%*%(rep(delta,p[d]))
      for (i in 1:length(phi)) {
        progress_counter <- progress_counter +1 #progress bar counter
        setWinProgressBar(pb, progress_counter,
                          label =paste( round(progress_counter/total*100, 0),
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
              t <- t+1
              x <- MASS::mvrnorm(n,mu = mean_shifted,Sigma = sigma0)
              
              Xt <- cbind(Xt,x)
              
              var_ez <- (phi[i]**2 + phi2[i,j]**2)*((1-
                        (1-phi[i]+phi2[i,j])^(2*t))/(2*(phi[i]-phi2[i,j]) -
                                                       (phi[i]-phi2[i,j])^2))
              var_ez <- var_ez - 2*(1-phi[i] + 
                                      phi2[i,j])*(phi[i]*phi2[i,j])*((1-
                        (1-phi[i]+phi2[i,j])^(2*(t-1)))/(2*(phi[i]-phi2[i,j]) -
                                                          (phi[i]-phi2[i,j])^2))
              var_ez <- var_ez*sigma0
              
              EZt <- phi[i]*Xt[,t] -phi2[i,j]*Xt_1 +(1-phi[i] +phi2[i,j])*EZt_1
              T2 <- t(EZt-mu0)%*%solve(var_ez)%*%(EZt-mu0)
              signal <- test.ooc(T2,h_eew[d,j,i])
              EZt_1 <- EZt
              Xt_1 <- Xt[,t]
              
            }#while IC
            count_vec[c] <- t
          }#for sims
          newrow <- cbind(newrow,
                          t(c(mean(count_vec),sd(count_vec),median(count_vec))))
        }#for phi2
      }#for phi
      output <- rbind(output,newrow)
    }#for shifts
    newnames <- t(cbind(rep(phi,each=3*3),
                        rep(phi2,each=3),
                        rep(h_eew[d,,],each=3),rep(c("ARL","SDRL","MRL"),2)))
    results_EEW[[d]] <- as.data.frame(rbind(newnames,output))
  }#for p 

  
#prep output
  for (d in 1:length(p)) {
    rownames(results_EW[[d]]) <- c("Phi","h","Measure",shifts)
    rownames(results_EEW[[d]]) <- c("Phi","Phi2","h","Measure",shifts)
    rownames(results_HW[[d]]) <- c("Phi","h","Measure",shifts)
  }
  
  file_name <- "Comparing Charts Multivariate cases.xlsx"
  #print results to excel
  #print MEWMA results
  write.xlsx(results_EW[[1]],
             sheetName = "MEWMA p=2",
             file = file_name,row.names=TRUE,col.names = FALSE)
  write.xlsx(results_EW[[2]],
             sheetName = "MEWMA p=3",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(results_EW[[3]],
             sheetName = "MEWMA p=4",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  #print extended MEWMA results
  write.xlsx(results_EEW[[1]],
             sheetName = "Extended MEWMA p=2",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(results_EEW[[2]],
             sheetName = "Extended MEWMA p=3",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(results_EEW[[3]],
             sheetName = "Extended MEWMA p=4",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  #print MHWMA results  
  write.xlsx(results_HW[[1]],
             sheetName = "MHWMA p=2",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(results_HW[[2]],
             sheetName = "MHWMA p=3",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(results_HW[[3]],
             sheetName = "MHWMA p=4",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  close(pb)
  