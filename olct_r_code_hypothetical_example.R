
###############################################################################
#### Below simulation is based on the hypothetical example
###############################################################################



library(MCPMod)
library(ggplot2)

# Generate the data from the underlying true curve
## reference code: https://onlinelibrary.wiley.com/doi/10.1002/bimj.201700276
simdatar <- function(ssize=100,doses=c(0,100,200,300,400),model="linear",sigma=1.5, e0, emax, e1, ed50, delta, b1, b2){
  dose <- rep(doses,ssize)
  ndose=length(doses)
  
  if (model=="linear"){
    mn <- e0+delta*doses
  }
  
  if (model=="emax"){
    mn <-emax(doses,e0,emax,ed50)
  }
  
  if (model=="quadratic"){
    mn<-quadratic(doses,e0,b1,b2)
  }
  
  if (model=="exponential"){
    mn<-exponential(doses,e0,e1,delta)
  }
  
  if (model=="logistic"){
    mn<-logistic(doses,e0, emax, ed50, delta)
    
  }
  
  mnVec <- rep(mn, ssize)
  out<-mnVec + rnorm(length(mnVec),0,sigma)
  dat <- as.data.frame(cbind(dose,out))
  colnames(dat)<-c("dose","resp")
  return(dat)
}





sim_power<-function(n0, n1, n2, n3, n4, 
                    targetdose, mcid, sd, truemodel, sim_n, 
                    candidate_models, doselevel, doserange){
  
  
  olct=olct_est_correct=olct_est_low=olct_est_high=rep(0, sim_n)
  olct_estci_correct=olct_estci_low=olct_estci_high=rep(0, sim_n)
  mcpmod=mcp_range_correct=mcp_correct=mcp_low=mcp_high=mcp_range_low=mcp_range_high=rep(0, sim_n)
  
  for (i in 1:sim_n){
    
    ### simulating data
   
    set.seed(i)
    
    if (truemodel=="linear") simdata<-simdatar(ssize=100,doses=c(0,100,200,300,400),model="linear",sigma=1.5, e0=0.2, delta=0.002) 
    if (truemodel=="emax") simdata<-simdatar(ssize=100,doses=c(0,100,200,300,400),model="emax",sigma=1.5, e0=0.2, ed50=200, emax=0.8)
    if (truemodel=="quadratic") simdata<-simdatar(ssize=100,doses=c(0,100,200,300,400),model="quadratic",sigma=1.5, e0=0.2, b1=0.006, b2=-0.000015)
    if (truemodel=="logistic") simdata<-simdatar(ssize=100,doses=c(0,100,200,300,400),model="logistic",sigma=1.5, e0=0.2, emax=0.8, ed50=254.7411, delta=111.6221)
    if (truemodel=="exponential") simdata<-simdatar(ssize=100,doses=c(0,100,200,300,400),model="exponential",sigma=1.5, e0=0.2, e1=0.2, delta=236.0445)

    
    dose<-simdata$dose
    resp=simdata$resp
    
    ### calculate mean and standard deviation
    mean=aggregate(simdata$resp, list(dose=simdata$dose), mean)
    std=aggregate(simdata$resp, list(dose=simdata$dose), sd)
    
    ### pre-defined model list
    models <- candidate_models
    
    
    # Fit the MCP-Mod
    fit.MCPMod = MCPMod(simdata, models, alpha = 0.05, pVal = TRUE, selModel = "maxT", doseEst = "MED2", clinRel = mcid, dePar = 0.05)
    
    # check if there is a significant model
    if (fit.MCPMod$signf=="TRUE") {
      
      mcpmod[i]=1
      
      # check if there is a dose recommended by the model
      if (names(unlist(fit.MCPMod)[length(unlist(fit.MCPMod))])=="tdose.MED2,90%"){
        
        # dose recommended by the model
        estdose<-unlist(fit.MCPMod)$`tdose.MED2,90%`  
        # select the minimum effective dose among the available doses
        mod_dos<-dose[dose>=estdose][1] 
        
        ## check if dose is non-null
        if (is.na(estdose)=="FALSE"){
          if (mod_dos==targetdose) mcp_correct[i]=1
          if (mod_dos<targetdose) mcp_low[i]=1
          if (mod_dos>targetdose) mcp_high[i]=1
         
         if (length(doselevel[doselevel==estdose])==1){
            if (estdose==targetdose) mcp_range_correct[i]=1
            if (estdose<targetdose) mcp_range_low[i]=1
            if (estdose>targetdose) mcp_range_high[i]=1
          }
          
          else if (length(doselevel[doselevel==estdose])==0){
            if (((estdose+doserange)>=targetdose) & (estdose<=targetdose)) mcp_range_correct[i]=1
            if (estdose+doserange<targetdose) mcp_range_low[i]=1
            if (estdose>targetdose) mcp_range_high[i]=1
        }

        if (is.na(estdose)=="TRUE"){
          mcp_correct[i]=0
          mcp_range_correct[i]=0
        }
        
      }
      
      if (names(unlist(fit.MCPMod)[length(unlist(fit.MCPMod))])!="tdose.MED2,90%"){
        mcp_correct[i]=0
        mcp_range_correct[i]=0
      }
      
      
    }
    
    if (fit.MCPMod$signf=="FALSE") {
      
      mcpmod[i]=mcp_correct[i]=0
      mcp_range_correct[i]=0
      
    }
    
    
    
    ## OLCT test
    fit.aovmean = lm(resp~as.factor(dose)-1, data=simdata)
    EstMean.aov = coef(fit.aovmean)
    aov.sigma = summary(fit.aovmean)$sigma
    Contrast.Coef = c(-2,-1,0,1,2)
    PoC.Contrast=sum(Contrast.Coef*EstMean.aov)/(aov.sigma*sqrt(sum(Contrast.Coef^2/c(n0, n1, n2, n3, n4))))>
      qt(1-0.05,df=n0+n1+n2+n3+n4-5)
    if (PoC.Contrast=="TRUE")  {
      olct[i]=1
      # dose recommended by the mean of dose > MCID + mean of placebo
      olct_dos<- mean[(mean$x>=(mcid+mean[mean$dose==0, 2]))&(mean$dose>0), 1][1]
      # dose recommended by the mean of dose > MCID + mean of placebo and 90% lower bound CI >= placebo effect
      olct2_dos<-mean[(mean$x>=(mcid+mean[mean$dose==0, 2]))&(mean$dose>0)&((mean$x-qnorm(0.95)*(std$x)/sqrt(n0))>=mean[mean$dose==0, 2]), 1][1]
      
    }
    if (PoC.Contrast=="FALSE") {
      olct[i]=0
      olct_dos<-NA
      olct2_dos<-NA
    }
    
    
    if (is.na(olct_dos)=="FALSE"){
      if (olct_dos==targetdose) olct_est_correct[i]=1
      if (olct_dos<targetdose) olct_est_low[i]=1
      if (olct_dos>targetdose) olct_est_high[i]=1
    }
    if (is.na(olct_dos)=="TRUE"){olct_est_correct[i]=0}
    
    if (is.na(olct2_dos)=="FALSE"){
      if (olct2_dos==targetdose) olct_estci_correct[i]=1
      if (olct2_dos<targetdose) olct_estci_low[i]=1
      if (olct2_dos>targetdose) olct_estci_high[i]=1
    }
    if (is.na(olct2_dos)=="TRUE"){ olct_estci_correct[i]=0}
    
    
    
  }
  
  olctpower<-mean(olct)
  mcpmodpower<-mean(mcpmod)
  prob_est_olct<-mean(olct_est_correct)
  prob_estci_olct<-mean(olct_estci_correct)
  prob_mcp<-mean(mcp_correct)
  prob_mcp_range<-mean(mcp_range_correct)
  
  prob_est_olct_low<-mean(olct_est_low)
  prob_est_olct_high<-mean(olct_est_high)
  prob_estci_olct_low<-mean(olct_estci_low)
  prob_estci_olct_high<-mean(olct_estci_high)
  prob_mcp_low<-mean(mcp_low)
  prob_mcp_high<-mean(mcp_high)
  prob_mcp_range_low<-mean(mcp_range_low)
  prob_mcp_range_high<-mean(mcp_range_high)
  
  
  c(doserange=doserange, 
    # power
    olctpower=olctpower,mcpmodpower=mcpmodpower, 
    # prob of correct dose
    prob_est_olct=prob_est_olct, prob_estci_olct=prob_estci_olct, 
    prob_mcp=prob_mcp, prob_mcp_range=prob_mcp_range, 
    # prob of low and high
    prob_est_olct_low=prob_est_olct_low, prob_est_olct_high=prob_est_olct_high,
    prob_estci_olct_low=prob_estci_olct_low, prob_estci_olct_high=prob_estci_olct_high, 
    prob_mcp_low=prob_mcp_low, prob_mcp_high=prob_mcp_high, 
    prob_mcp_range_low=prob_mcp_range_low, prob_mcp_range_high=prob_mcp_range_high)
  
}



#####################################################################
######## example 
#####################################################################


#### linear model
## including linear model
sim_result1<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   200, 0.3, 1.5, "linear", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          exponential =guesst (d= 200, p=0.3 , model ="exponential", Maxd = 400),
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))


# excluding linear model
sim_result2<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   200, 0.3, 1.5, "linear", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))




### emax model
#including emax model
sim_result_emax1<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   200, 0.3, 1.5, "emax", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          exponential =guesst (d= 200, p=0.3 , model ="exponential", Maxd = 400),
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))

# excluding emax model
sim_result_emax2<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   200, 0.3, 1.5, "emax", 1000, 
                   list ( linear=NULL,
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          exponential =guesst (d= 200, p=0.3 , model ="exponential", Maxd = 400),
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))



### quad model
### including quadratic model
sim_result_quad1<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   100, 0.3, 1.5, "quadratic", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          exponential =guesst (d= 200, p=0.3 , model ="exponential", Maxd = 400),
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))

# excluding quadratic model
sim_result_quad2<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   100, 0.3, 1.5, "quadratic", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          exponential =guesst (d= 200, p=0.3 , model ="exponential", Maxd = 400)),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))




### logistic model
# including logistic model
sim_result_log1<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   300, 0.3, 1.5, "logistic", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))

# excluding logistic model
sim_result_log2<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   300, 0.3, 1.5, "logistic", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))




### exponential model
# including exponential model
sim_result_exp1<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   300, 0.3, 1.5, "exponential", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          exponential =guesst (d= 200, p=0.3 , model ="exponential", Maxd = 400),
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))

# excluding exponential model
sim_result_exp2<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(100, 100, 100, 100, 100,
                   300, 0.3, 1.5, "exponential", 1000, 
                   list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
                          linear=NULL,
                          logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
                          quadratic=guesst(d = 200, p = 1, "quadratic")),
                   c(0, 100, 200, 300, 400), 10*i)
  c(result=result)
}))



###############################################################################
####################### underlying curve parameters ###########################
###############################################################################

### emax model
### model function: f(d,theta)=E0+Emax d/(ED50 + d)
### r function: emax(dose, e0, eMax, ed50)
emax=guesst(d=200, p=0.5, model ="emax")
emax_resp=emax(c(0,100,200,300,400),0.2, 0.8, 200)
emax_curve=list(emax=guesst(d= 200, p=0.5, model ="emax")) 
plotModels(emax_curve, doses=c(0, 100, 200, 300, 400), base = 0.2,
           maxEff=0.8)


### linear model
### model function: f(d,theta)=E0+delta d
### r function: linear(dose, e0, delta)
linear_resp=linear(c(0,100,200,300,400), 0.2, 0.002)
linear_curve=list(linear =NULL)
plotModels(linear_curve, doses=c(0, 100, 200, 300, 400), base = 0.2,
           maxEff=0.8)


### Quadratic model
### model function: f(d,theta)=E0+beta1 d+beta2 d^2
### r function: quadratic(dose, e0, b1, b2)
### delta=b2/abs(b1)
quadratic=guesst(d = 200, p = 1, "quadratic")
quadratic_resp=quadratic(c(0,100,200,300,400),0.2,0.006,-0.000015)
quad_curve=list(quadratic=guesst(d = 200, p = 1, "quadratic")) 
plotModels(quad_curve, doses=c(0, 100, 200, 300, 400), base = 0.2,
           maxEff=0.8)


### logistic model
### model function: f(d,theta)=E0+Emax/(1 + exp((ED50-d)/delta)).
### r function: logistic(dose, e0, eMax, ed50, delta)
logistic=guesst(d=c(100, 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400)
logistic_resp=logistic(c(0,100,200,300,400), 0.2, 0.8, 254.7411, 111.6221)
logistic_curve=list(logistic=guesst(d=c(100, 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400))
plotModels(logistic_curve, doses=c(0, 100, 200, 300, 400), base = 0.2,
           maxEff=0.8)


### exponential model
### model function: f(d,theta)=E0+E1 (exp(d/delta)-1).
### r function: exponential(dose, e0, e1, delta) e1	E1 parameter
exponential=guesst(d=200, p=0.3, model ="exponential", Maxd = 400)
exp_resp=exponential(c(0, 100, 200, 300, 400), 0.2, 0.2, 236.0445)
exp_curve=list(exponential=guesst(d=200, p=0.3, model ="exponential", Maxd = 400))
plotModels(exp_curve, doses=c(0, 100, 200, 300, 400), base = 0.2,
           maxEff=0.8)




#####################################################################
######## shape of candidate models 
#####################################################################

modlist=list ( emax =guesst (d= 200, p=0.5 , model ="emax"),
               linear =NULL, 
               exponential =guesst (d= 200, p=0.3, model ="exponential", Maxd = 400),
               logistic = guesst (d=c(100 , 300), p=c (0.2, 0.6) ,"logistic", Maxd = 400), 
               quadratic=guesst(d = 200, p = 1, "quadratic"))
plotModels(modlist, doses=c(0, 100, 200, 300, 400), base = 0.2,
           maxEff=0.8)


