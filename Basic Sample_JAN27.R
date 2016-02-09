set.seed(10072015)

n <- 10000              ## Number of subjects
sigma.rw <- 0.1         ## Outpatient Variability
BETA.oi <- c(-6.5, 0.1) ## Outpatient -> Inpatinet.Ill Coefficients    
BETA.ir <- c(-0.2, -0.1)## Inpatient.Ill -> Inpatient.Recovering Coefficients
BETA.ri <- c(-2, -0.05) ## Inpatient.Recovering -> Inpatient.Ill Coefficients
BETA.ro <- c(0, -0.005) ## Inpatient.Recovering -> Outpatient Coefficients
BETA.ev <- c(-14, 0.25) ## Probability of Event Coeffeicients
Beta.obs <- c(2, -0.1)  ## Probability of Observation Coefficients
dtO <- 1                ## Time between of Outpatient random walk 
dtI <- 1                ## Time between Inpatient observations
mu.inp <- 0.75          ## Inpatient.Ill Drift
mu.rec <- -0.8          ## Inpatient.Recovering Drift
sigma.inp <- 1.2        ## Inpatient.Ill Variability
sigma.rec <- 1.25       ## Inpatient.Recovering Variability
BETA.ev2 <- NULL        ## Competing Event's Event Coeffeicients 
FollowTime <- 800       ## Time in Simulation
lambda <- 5             ##  
B <- c(3.8,-1.9) ## Outpatient Visits  ARE Roughly 40 Days apart ##



EHRDat <- vector("list", length = n)
for(i in 1:n){
  Ev <- 0; State <- "Outpatient";
  BMVALUES <- STATE <- TIME <- NULL
  BM <- 0
  Time <- 0; 
  while(Ev == 0 & tail(Time,1) <= FollowTime){
    if(State == "Outpatient"){
      BM <- BM 
      BM <- BM + rnorm(1,0,sigma.rw)
      # BMvalues <- c(BMvalues, BM)
      X <- c(1,BM)
      ###Probability of Switch to Inpatient
      Poi <- 1/(1 + exp(-X%*%BETA.oi))
      Switch <- rbinom(1,1,Poi)
      State <- ifelse(Switch == 1, "Inpatient.Ill", State)
    } 
    else if(State == "Inpatient.Ill"){
      BM <- BM + rnorm(1,mu.inp,sigma.inp)
      # BMvalues <- c(BMvalues, BM)
      X <- cbind(1,BM)
      ###Probability switch to Inpatient Recovery
      Pir <- 1/(1 + exp(-X%*%BETA.ir))
      Switch <- rbinom(1,1,Pir)
      State <- ifelse(Switch == 1, "Inpatient.Recovery", State)
    } 
    else if(State == "Inpatient.Recovery"){
      BM <- BM + rnorm(1,mu.rec,sigma.rec) 
      #  BMvalues <- c(BMvalues, BM)
      X <- cbind(1,BM)
      ###Probability switch recovery to outpatient
      Pro <- exp(X%*%BETA.ro)/(1 + (exp(X%*%BETA.ro) + exp(X%*%BETA.ri)))
      ###Probability switch recovery to illness
      Pri <- exp(X%*%BETA.ri)/(1 + (exp(X%*%BETA.ro) + exp(X%*%BETA.ri)))
      Switch <- rmultinom(1,1,c(Pro,Pri,(1- Pro - Pri)))
      State <- c("Outpatient","Inpatient.Ill","Inpatient.Recovery")[which(Switch == 1)]
    }
    
    STATE <- c(STATE, State)
    BMVALUES <-c(BMVALUES,BM)
    Time <- Time + ifelse(State == "Outpatient", dtO, dtI)
    TIME <- c(TIME, Time)
    
    Pev <- 1/(1 + exp(-X%*%BETA.ev))
    Ev <- rbinom(1,1,Pev)
    
    
  }
  EHRDat[[i]] <- data.frame(BMvalues = BMVALUES, State = STATE, Time = TIME, Event = Ev)
}

#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################


## The rest of the code is for figures 


### INFORMED PRESENCE ##
EHRDat.Obs <- vector("list", length = n)
### Inpatient Encounters ###
for(i in 1:n){
  State <- EHRDat[[i]]$State
  Inp <- which(State %in% c("Inpatient.Ill", "Inpatient.Recovery"))
  if(length(Inp)>0){
    Inp.Start <- Inp[c(1,which(diff(Inp) > 1) + 1)] ###Identify the beginning of Inpatient State
    Inpr <- rev(Inp)
    Inp.End <- rev(Inpr[c(1,which(-diff(Inpr) > 1) + 1)]) ###Identify the end of Inpatient State
    Inp.State <- cbind(Inp.Start,Inp.End)
    x <- t(rbind(rep(1,length(EHRDat[[i]]$BMvalues['Time'=Inp.Start])),EHRDat[[i]]$BMvalues['Time'=Inp.Start])) #*#
    pInp <- 1/(1 + exp(-x%*%Beta.obs)) ### probability that visit is observed
    Inp.Obs <- rbinom(nrow(Inp.State), 1, pInp)
    Inp.State.r <- Inp.State[Inp.Obs == 1,] ###Get the indices for the inpatient states to KEEP
    if(length(Inp.State.r)!=0)
    {
      if (length(Inp.State.r)!=2)
      {
        InpState.Keep <- do.call("c", as.list(apply(Inp.State.r, 1, function(x)seq(x[1],x[2]))))
      } else{
        InpState.Keep <- seq(Inp.State.r[1],Inp.State.r[2])
      }
    } else{
      InpState.Keep <- NULL
    }
  }else{
    InpState.Keep <- NULL
  }
  ## Outpatient Encounters ##
  Out <- which(State %in% c("Outpatient"))
  if(length(Out)>0){  
    t <- Out[1]
    Out.Keep <- t
    x <- t(rbind(rep(1,length(EHRDat[[i]]$BMvalues['Time'=Out])),EHRDat[[i]]$BMvalues['Time'=Out])) #*#
    
    while(t < nrow(x)){
      
      x[,2] <- ifelse(x[,2]<0,-0.75,1.25) 
      lambda <- exp(x[t,]%*%B)
      
      #  lambda <- ifelse(lambda<0,40,lambda)
      t <- t + rpois(1,lambda)
      if (t < nrow(x)){
        if (EHRDat[[i]]$State['Time'=t]=='Outpatient'){
          Out.Keep <- cbind(Out.Keep,t)
        }
      }
    }    
    if(Out.Keep[length(Out.Keep)]>Out[length(Out)]){Out.Keep <-Out.Keep[-length(Out.Keep)]}
    EHRDat.Obs[[i]] <- EHRDat[[i]][c(InpState.Keep,Out.Keep),]
  }
}




all <- c()
for(i in 1:n){
  if(NCOL(EHRDat.Obs[[i]])>1){
    EHRDat.Obs[[i]]$PATID <- i
    all <- rbind(all,EHRDat.Obs[[i]])
    print(i)
  }
}
all<- all[order(all$PATID,all$Time),]
library(dplyr)
all <- all %>% distinct(PATID, Time)

set.seed(101)
all$BMvalues <- all$BMvalues + rbinom(nrow(all),1,0.25)*runif(nrow(all),-25,100)
## ADD in variability

save(all, file = "C:/Users/mp282/Dropbox/EHRSimulation/Dat/sim.RData")



