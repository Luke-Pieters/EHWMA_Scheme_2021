
# EWMA chart for steady-state case


require(matrixStats)
m0=0
s=1
delta=0.75
m1=m0+delta*s
rep=50000
tau=10
subgroups=10000
n=1
w=0.10
L=2.813
#UCL=
#LCL=
#ARL_sim_EWMA_chart<-function(rep,subgroups,n,l,L,w,UCL,LCL)
ARL_sim_EWMA_chart<-function(rep,subgroups,n,l,L,w)
{
  
  
  z<-rep(0,subgroups)
  xi<-rep(0,subgroups)
  UCL<-rep(0,subgroups)
  LCL<-rep(0,subgroups)
  
  for(j in 1:subgroups)
  {
    UCL[j]<-m0+L*(s/sqrt(n))*sqrt((w/(2-w))*(1-(1-w)^(2*j)))
    #UCL[j]<-m0+L*(s/sqrt(n))*sqrt((w/(2-w)))
    LCL[j]<--UCL[j]
  }
  
  
  #vector RL values for each repetition
  RL<-rep(0,rep)
  set.seed(1098)
  for(i in 1:rep)
  {
    total<-subgroups*n
    
    
    sub1a<-matrix(rnorm(n*tau,m0,s),ncol=n,nrow=tau)
    sub1b<-matrix(rnorm(total-(n*tau),m1,s),ncol=n,nrow=subgroups-tau)
    sub2<-rbind(sub1a,sub1b)
    
    
    
    
    for(j in 1:subgroups)
    {
      
      #calculation of x
      xi[j]<-mean(sub2[j,])
      
      
      #calculation of z
      
      if(j==1)
      {
        z[1]<-w*xi[1]+(1-w)*m0
      }
      else
      {
        z[j]<-w*xi[j]+(1-w)*z[j-1]
      }
      
      
      
      if(z[j]>=UCL[j] && j>tau|z[j]<=LCL[j] && j>tau)
      {
        RL[i]=j-tau;
        break;
      }
      else
      {
        RL[i]=0;
      }
      
      
      if(j==subgroups && RL[i]==0)
      {
        RL[i]<-subgroups
        break
        
      }
      
      
    }
    
  }
  
  #ARL
  
  ARL<-mean(RL)
  SDRL<-sd(RL)
  QUA<-quantile(RL, prob = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
  
  
  #return(ARL)
  return(list(ARL=ARL,SDRL=SDRL,QUA=QUA))
  
}
ARL_sim_EWMA_chart(rep,subgroups,n,l,L,w)
