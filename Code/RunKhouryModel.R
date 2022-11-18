###RunKhouryModel



##Import table
SummaryTable_Efficacy_NeutRatio_SD_SEM<-read.csv(paste0(currentdirectory,"/Data/SummaryTable_Efficacy_NeutRatio_SD_SEM.csv"),fileEncoding="UTF-8-BOM")


###Tidying up Summary Table.
SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)
SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_LB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)-1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)
SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_UB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)+1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)
SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName[SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName==" rAd26-S+rAd5-S"]="Gam-COVID-Vac"



####Logistic Model from Khoury et al.
ProbRemainUninfected=function(logTitre,logk,C50){1/(1+exp(-exp(logk)*(logTitre-C50)))}

LogisticModel_PercentUninfected=function(mu_titre,sig_titre,logk,C50){
  NumInteration<-max(length(mu_titre),length(sig_titre),length(logk),length(C50))
  Output<-NULL
  
  if (length(C50)==1) {
    C50=rep(C50,NumInteration)
  }
  
  if (length(logk)==1) {
    logk=rep(logk,NumInteration)
  }
  
  if (length(sig_titre)==1) {
    sig_titre=rep(sig_titre,NumInteration)
  }
  
  if (length(mu_titre)==1) {
    mu_titre=rep(mu_titre,NumInteration)
  }
  
  for (i in 1:NumInteration) {
    Step=sig_titre[i]*0.001
    IntegralVector=seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i]=sum(ProbRemainUninfected(IntegralVector,logk[i],C50[i])*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}


### Logistic model for Raw Efficacy Counts (Khoury et al.)
FittingLogistic_Raw<-function(logRisk0,logk,C50,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector){
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(N_C))
  }
  
  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  if (length(C50)==1) {
    C50=rep(C50,length(N_C))
  }
  
  if (length(logk)==1) {
    logk=rep(logk,length(N_C))
  }
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk[i],C50[i]))))
  }
  LL
}


### Fitting Model from Khoury et al - in order to get hessian/covariance matrix
LogisticEstimate=c("logk"=log(2.7),"C50"=log10(0.5))
FittedLogistic_RawEfficacy_MeanRept_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                                                SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)

### Confidence bounds (Khoury et al.):
Cov<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian)[9:10,9:10]


#####Evaluating the Indiviual Risk curve from Khoury et al.
NeutValuelog10New<-seq(log10(0.005),log10(64),by=0.001)
EstimatedEfficacy<-data.frame("NeutValuelog10New"=NeutValuelog10New,"ProbRemainUninfected_var"=ProbRemainUninfected(NeutValuelog10New,tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2]))


####Confindence Intervals of individual risk curve from Khoury et al.
EstimatedEfficacy$Lower<-NA
EstimatedEfficacy$Upper<-NA

N=10000
TempParameters<-rmvnorm(N,mean=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)),sigma=Cov)


for (i in 1:nrow(EstimatedEfficacy)) {
  tempEfficacyRange<-ProbRemainUninfected(EstimatedEfficacy$NeutValuelog10New[i],TempParameters[,1],TempParameters[,2])  
  EstimatedEfficacy$Lower[i]<-quantile(tempEfficacyRange,0.025)
  EstimatedEfficacy$Upper[i]<-quantile(tempEfficacyRange,0.975)
}
