#Robustness MULTIVARIATE CASES
library(xlsx) #write to excel
library(mvtnorm)
set.seed(1234)

#Set up
  sims <- 2000
  n <- 1 #sample size
  phi1 <- c(0.1,0.5)
  phi2 <- matrix(c(0.05,0.09,
                   0.05,0.25),byrow = T,nrow = 2,ncol = 2)
  p <- c(2,4,10)
  
  h <- list()
  h[[1]] <- matrix(c(9.14,9.44,
                10.37,10.38),byrow = T,nrow = 2,ncol = 2)
  h[[2]] <- matrix(c(13.29,13.62,
                 14.60,14.60),byrow = T,nrow = 2,ncol = 2)
  h[[3]] <- matrix(c(23.3,23.8,
                     24.95,24.97),byrow = T,nrow = 2,ncol = 2)
  
  percentiles <- c(0.05,0.25,0.5,0.75,0.95)
  
  dists <- c("N(0,1)","t(10)","t(100)","t(1000)",
             "GAM(1,1)","GAM(10,1)","LogNorm(0,1)","X2(30)")
  
  means <- c(0,0,0,0,1,10,exp(0.5),30)
  vars <- c(1,1.25,50/49,1000/998,1,10,exp(1)*(exp(1)-1),60)
  
  output <- list()
  output[[1]] <- matrix()
  output[[2]] <- matrix()
  output[[3]] <- matrix()
  
  #test if OOC function
  test.ooc <- function(x,cl){
    if ((x>=cl)) {
      TRUE
    }else{
      FALSE
    }#ifelse
  }#test function
  
  #random sample functions
  rdist.functions <- function(index,p){
    case_when(index == 1 ~ rnorm(n=p,0,1),
              index == 2 ~ rt(n=p,df = 10),
              index == 3 ~ rt(n=p,df = 100),
              index == 4 ~ rt(n=p,df = 1000),
              index == 5 ~ rgamma(n=p,shape = 1,rate=1),
              index == 6 ~ rgamma(n=p,shape = 10,rate=1),
              index == 7 ~ rlnorm(n=p,meanlog = 0,sdlog = 1),
              index == 8 ~ rchisq(n=p,df = 30)
    )
  }#rdist function
  total <- length(p)*length(phi1)*ncol(phi2)*length(dists)
  progress_counter <- 0
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = total, width = 300)
  
  
  
#run sims
  for (d in 1:length(p)) {
    outputbuild <- vector()
    newnames <- vector()
    for (i in 1:length(phi1)) {
      for (j in 1:ncol(phi2)) {
        newblocks <- vector()
        newnames <- cbind(newnames,
                          matrix(rep(c(phi1[i],phi2[i,j],h[[d]][i,j]),7),
                                 nrow = 3,ncol = 7,byrow = FALSE))
        for (m in 1:length(dists)) {
          count_vec <- vector()
          mu0 <- rep(means[m],p[d])
          sigma0 <- diag(rep(vars[m],p[d]))
          
          progress_counter <- progress_counter +1
          setWinProgressBar(pb, progress_counter,
                            title=paste( round(progress_counter/total*100, 0),
                                                               "% done"))
          for (c in 1:sims) {
            t <- 0
            signal <- FALSE
            Yt <- vector()
            Yt_1 <- mu0
            Yt_bar <- mu0
            while ((t<100000)&(signal==FALSE)) {
              t <- t+1
              z <- rdist.functions(index = m,p = p[d])
              
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
              signal <- test.ooc(T2,h[[d]][i,j])
              Yt_1 <- Yt[,t]
              Yt_bar <- rowMeans(Yt)
              
            }#while IC
            count_vec[c] <- t
          }#for sims
          
          newrow <- c(mean(count_vec),sd(count_vec),
                      quantile(count_vec,probs = percentiles))
          newblocks <- cbind(newblocks,newrow)
        }#for dists
        
        outputbuild <- cbind(outputbuild,t(newblocks))
      }#for phi2
    }#for phi1
    output[[d]] <- as.data.frame(rbind(newnames,outputbuild))
    rownames(output[[d]]) <- c("Phi1","Phi2","h",dists)
  }#for p
  
  file_name <- "Robustness Mulrivariate case.xlsx"
  
  #print results to excel
  write.xlsx(output[[1]],
             sheetName = "P=2",
             file = file_name,row.names=TRUE,col.names = FALSE)
  write.xlsx(output[[2]],
             sheetName = "P=4",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  write.xlsx(output[[3]],
             sheetName = "P=10",
             file = file_name,row.names=TRUE,col.names = FALSE,append = T)
  close(pb)
  
  