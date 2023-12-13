
###############################################################################
#### Below simulation is based on the data from below published clinical trial 
## ref: Shah, SJ., Voors, AA.,McMurrary, JJV., Kitzman, DW., Viethen, T., Wirtz, AB., Huang, E., Pap, AFP., Solomon, SD. (2019). 
##### Effect of Neladenoson Bialanate on Exercise Capacity Among Patients With Heart Failure With Preserved Ejection Fraction A Randomized Clinical Trial. 
##### Journal of the American Medical Association. 321 (21): 2101-2112.
###############################################################################

library(MCPMod)
library(ggplot2)
library(ggpubr)


# Generate the data from the underlying true curve
## reference code: https://onlinelibrary.wiley.com/doi/10.1002/bimj.201700276
simdatar <- function(ssize=60,doses=c(0,5,10,20,30,40),model="linear",sigma=80, e0, emax, e1, ed50, delta, b1, b2, h){
  dose <- rep(doses,ssize)
  ndose=length(doses)
  
  if (model=="linear"){
    mn <- e0+delta*doses
  }
  
  if (model=="emax"){
    mn <-emax(doses,e0,emax,ed50)
  }
  
  if (model=="sigEmax"){
    eMax=emax
    mn<-sigEmax(doses, e0, eMax, ed50, h)
    
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





sim_power<-function(n0, n1, n2, n3, n4, n5, 
                    targetdose, mcid, sd, truemodel, sim_n, 
                    candidate_models, doselevel, doserange){
  
  
  olct=olct_est_correct=olct_est_low=olct_est_high=rep(0, sim_n)
  olct_estci_correct=olct_estci_low=olct_estci_high=rep(0, sim_n)
  mcpmod=mcp_range_correct=mcp_correct=mcp_low=mcp_high=mcp_range_low=mcp_range_high=rep(0, sim_n)
  
  for (i in 1:sim_n){
    
    ### simulating data
    
    set.seed(i)
    
    if (truemodel=="linear") simdata<-simdatar(ssize=n0,doses=c(0,5,10,20,30,40),model="linear",sigma=80, e0=0, delta=1) 
    if (truemodel=="emax") simdata<-simdatar(ssize=n0,doses=c(0,5,10,20,30,40),model="emax",sigma=80, e0=0, ed50=1.25, emax=41.25)
    if (truemodel=="quadratic") simdata<-simdatar(ssize=n0,doses=c(0,5,10,20,30,40),model="quadratic",sigma=80, e0=0, b1=2.667, b2=-0.044)
    if (truemodel=="sigemax1") simdata<-simdatar(ssize=n0,doses=c(0,5,10,20,30,40),model="sigEmax",sigma=80, e0=0, emax=40.1, ed50=9, h=4)
    if (truemodel=="sigemax2") simdata<-simdatar(ssize=n0,doses=c(0,5,10,20,30,40),model="sigEmax",sigma=80, e0=0, emax=45, ed50=20, h=3)
    
    
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
          # if (abs(estdose-targetdose)<=doserange) mcp_range_correct[i]=1
          # if (estdose+doserange<targetdose) mcp_range_low[i]=1
          # if (estdose-doserange>targetdose) mcp_range_high[i]=1
          
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
    Contrast.Coef = c(-5, -3, -1, 1, 3, 5)
    PoC.Contrast=sum(Contrast.Coef*EstMean.aov)/(aov.sigma*sqrt(sum(Contrast.Coef^2/c(n0, n1, n2, n3, n4, n5))))>
      qt(1-0.05,df=n0+n1+n2+n3+n4+n5-6)
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

sigemax1=guesst(d= c(9, 40), p=c(0.5, 40/40.1), model ="sigEmax")
sigemax2=guesst(d= c(20, 40), p=c(0.5, 40/45),  model ="sigEmax")



#### linear model
## including linear model
sim_result1_linear<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   20, 15, 80, "linear", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"),
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result1_linear


#### linear model
## excluding linear model
sim_result0_linear<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   20, 15, 80, "linear", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          quadratic=guesst(d = 30, p = 1, "quadratic"),
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result0_linear





#### emax model
## including emax model
sim_result1_emax<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   5, 15, 80, "emax", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"),
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result1_emax



sim_result0_emax<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   5, 15, 80, "emax", 1000, 
                   list ( linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"),
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result0_emax




#### quadratic model
## including quadratic model
sim_result1_quad<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   10, 15, 80, "quadratic", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"),
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result1_quad




## excluding quadratic model
sim_result0_quad<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   10, 15, 80, "quadratic", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result0_quad



#### sigEmax1 model
## including sigemax1 model
sim_result1_sigemax1<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   10, 15, 80, "sigemax1", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"),
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result1_sigemax1




## excluding sigemax1 model
sim_result0_sigemax1<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   10, 15, 80, "sigemax1", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"), 
                          sigEmax=sigemax2),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result0_sigemax1



#### sigEmax2 model
## including sigemax2 model
sim_result1_sigemax2<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   20, 15, 80, "sigemax2", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"),
                          sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2)),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result1_sigemax2


## excluding sigemax2 model
sim_result0_sigemax2<-do.call(rbind, lapply(1:9, function(i) {
  result=sim_power(60, 60, 60, 60, 60, 60,
                   20, 15, 80, "sigemax2", 1000, 
                   list ( emax=guesst(d= 1.25, p=0.5, model ="emax"),
                          linear=NULL,
                          quadratic=guesst(d = 30, p = 1, "quadratic"), 
                          sigEmax=sigemax1),
                   c(0,5,10,20,30,40), i)
  c(result=result)
}))

sim_result0_sigemax2





#######################################
########### Curve #####################
#######################################
modlist=list(linear=NULL,
             sigEmax=matrix(c(sigemax1, sigemax2), byrow=T, ncol=2),
             emax=guesst(d= 1.25, p=0.5, model ="emax"),
             quadratic=guesst(d = 30, p = 1, "quadratic"))
plotModels(modlist, doses=c(0,5,10,20,30,40), base = 0, maxEff=40,
           xlab="Dose (mg)", ylab="6MWD absolute change from baseline (m)")




###############################################################################
####################### underlying curve parameters ###########################
###############################################################################

### emax model
### model function: f(d,theta)=E0+Emax d/(ED50 + d)
### r function: emax(dose, e0, eMax, ed50)
emax_resp=emax(c(0,5,10,20,30,40),0, 41.25, 1.25)


#### sigemax model
sigEmax1=sigEmax(c(0,5,10,20,30,40), 0, 40.1, 9, 4)
sigEmax2=sigEmax(c(0,5,10,20,30,40), 0, 45, 20, 3)



### linear model
### model function: f(d,theta)=E0+delta d
### r function: linear(dose, e0, delta)
linear_resp=linear(c(0,5,10,20,30,40), 0, 1)


### Quadratic model
### model function: f(d,theta)=E0+beta1 d+beta2 d^2
### r function: quadratic(dose, e0, b1, b2)
### delta=b2/abs(b1)
quadratic_resp=quadratic(c(0,5,10,20,30,40),0,2.667,-0.044)













########################################################
###########  correct dose figure ######################
#######################################################

plot_dose_correct<-function(data){
  
  
  ############ correct dose selection
  ## reshape to longitudinal data structure
  plot_data_correct<-data.frame(doserange=rep(data[, 1], 3),
                                method=c(rep("olct", 9), rep("mcpmod", 9), rep("mcpmod_range", 9)),
                                prob=c(data[, 5],data[, 6],data[, 7]))
  
  #plot prob by method
  figure_dose_correct=ggplot(plot_data_correct, aes(x=doserange, y=prob)) + 
    geom_line(aes(color=method)) +
    scale_color_manual(name='method', labels=c('MCPMod', 'MCPMod range', 'OLCT' ),
                       values=c('red', 'blue', 'green'))+
    labs("Probability of recommending target dose correctly",
         y= "Probability of selecting dose correctly", 
         x = "Dose Recomendation Range (+/-)")+
    scale_x_continuous(breaks=seq(1,9,1), limits=c(1,9))+
    scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1))
  
  figure_dose_correct
  
}



figure_correct=ggarrange(
  ggarrange(plot_dose_correct(sim_result1_linear)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_correct(sim_result0_linear)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2), 
  ggarrange(plot_dose_correct(sim_result1_emax)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_correct(sim_result0_emax)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_correct(sim_result1_quad)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_correct(sim_result0_quad)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_correct(sim_result1_sigemax1)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_correct(sim_result0_sigemax1)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_correct(sim_result1_sigemax2)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_correct(sim_result0_sigemax2)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  nrow=5)

annotate_figure(figure_correct,
                bottom = text_grob("Dose Range (+)"),
                left = text_grob("Probability of selecting dose correctly", rot = 90))






########################################################
###########  low dose figure ######################
#######################################################

plot_dose_low<-function(data){
  
  
  ############ correct dose selection
  ## reshape to longitudinal data structure
  plot_data_low<-data.frame(doserange=rep(data[, 1], 3),
                            method=c(rep("olct", 9), rep("mcpmod", 9), rep("mcpmod_range", 9)),
                            prob=c(data[, 10],data[, 12],data[, 14]))
  
  #plot prob by method
  figure_dose_low=ggplot(plot_data_low, aes(x=doserange, y=prob)) + 
    geom_line(aes(color=method)) +
    scale_color_manual(name='method', labels=c('MCPMod', 'MCPMod range', 'OLCT' ),
                       values=c('red', 'blue', 'green'))+
    labs("Probability of selecting dose lower than true target dose",
         y= "Probability of selecting dose lower than true target dose", 
         x = "Dose Recomendation Range (+/-)")+
    scale_x_continuous(breaks=seq(1,9,1), limits=c(1,9))+
    scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1))
  
  figure_dose_low
  
}



figure_low=ggarrange(
  ggarrange(plot_dose_low(sim_result1_linear)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_low(sim_result0_linear)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2), 
  ggarrange(plot_dose_low(sim_result1_emax)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_low(sim_result0_emax)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_low(sim_result1_quad)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_low(sim_result0_quad)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_low(sim_result1_sigemax1)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_low(sim_result0_sigemax1)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_low(sim_result1_sigemax2)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_low(sim_result0_sigemax2)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  nrow=5)

annotate_figure(figure_low,
                bottom = text_grob("Dose Range (+)"),
                left = text_grob("Probability of selecting dose lower than true target dose", rot = 90))






########################################################
###########  high dose figure ######################
#######################################################

plot_dose_high<-function(data){
  
  
  ############ correct dose selection
  ## reshape to longitudinal data structure
  plot_data_high<-data.frame(doserange=rep(data[, 1], 3),
                             method=c(rep("olct", 9), rep("mcpmod", 9), rep("mcpmod_range", 9)),
                             prob=c(data[, 11],data[, 13],data[, 15]))
  
  #plot prob by method
  figure_dose_high=ggplot(plot_data_high, aes(x=doserange, y=prob)) + 
    geom_line(aes(color=method)) +
    scale_color_manual(name='method', labels=c('MCPMod', 'MCPMod range', 'OLCT' ),
                       values=c('red', 'blue', 'green'))+
    labs("Probability of selecting dose higher than true target dose",
         y= "Probability of selecting dose higher than true target dose", 
         x = "Dose Recomendation Range (+/-)")+
    scale_x_continuous(breaks=seq(1,9,1), limits=c(1,9))+
    scale_y_continuous(breaks=seq(0,0.6,0.2), limits=c(0,0.6))
  
  figure_dose_high
  
}



figure_high=ggarrange(
  ggarrange(plot_dose_high(sim_result1_linear)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_high(sim_result0_linear)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2), 
  ggarrange(plot_dose_high(sim_result1_emax)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_high(sim_result0_emax)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_high(sim_result1_quad)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_high(sim_result0_quad)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_high(sim_result1_sigemax1)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_high(sim_result0_sigemax1)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  ggarrange(plot_dose_high(sim_result1_sigemax2)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            plot_dose_high(sim_result0_sigemax2)+ 
              theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
            ncol = 2),
  nrow=5)

annotate_figure(figure_high,
                bottom = text_grob("Dose Range (+)"),
                left = text_grob("Probability of selecting dose higher than true target dose", rot = 90))






