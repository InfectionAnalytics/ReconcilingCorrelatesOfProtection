
### Plot figure 3

##Importing Bergwerk data
BergwerkData<-read.csv(paste0(currentdirectory,"/Data/Bergwerk_FalseData.csv"),fileEncoding="UTF-8-BOM")

##Efficacy estimate
Efficacy_Chodick<-0.94 #Chodick paper efficacy around the same time period

####BergwerkNormalisation
###Chodick study

#Estimate LOD for Bergwerk, and group all titres at the LOD together at >0 (for plotting purposes)
LODestimate<-min(c(BergwerkData$PeakNeutTiter[BergwerkData$PeakNeutTiter>0],BergwerkData$PeriInfectionNeutTiter[BergwerkData$PeriInfectionNeutTiter>0]),na.rm=TRUE)
BergwerkData$PeakNeutTiter[BergwerkData$PeakNeutTiter==0]=LODestimate
BergwerkData$PeriInfectionNeutTiter[BergwerkData$PeriInfectionNeutTiter==0]=LODestimate


ConversionFactorToConvalescence=(10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Pfizer"])

PeriPfizer_Control<-10^(mean(log10(BergwerkData$PeriInfectionNeutTiter[BergwerkData$Infected==0]),na.rm=TRUE))
BergwerkData_Plot<-BergwerkData
BergwerkData_Plot$PeriInfectionNeutTiter_NormControl<-BergwerkData_Plot$PeriInfectionNeutTiter/(PeriPfizer_Control/ConversionFactorToConvalescence)

Bergwerkdata_subset_Control_Plot<-BergwerkData_Plot[c("PeriInfectionNeutTiter_NormControl","Infected")]
Bergwerkdata_subset_Control_Plot$countcolumn<-1

Bergwerkplot_Control_Plot<-aggregate(.~PeriInfectionNeutTiter_NormControl,data=Bergwerkdata_subset_Control_Plot,FUN=function(x){sum(x)})
Bergwerkplot_Control_Plot$FractionInfected<-(1-Efficacy_Chodick)*(Bergwerkplot_Control_Plot$Infected/sum(Bergwerkplot_Control_Plot$Infected))/((Bergwerkplot_Control_Plot$countcolumn-Bergwerkplot_Control_Plot$Infected)/sum(Bergwerkplot_Control_Plot$countcolumn-Bergwerkplot_Control_Plot$Infected))

Bergwerkplot_Control_Plot$Control<-Bergwerkplot_Control_Plot$countcolumn-Bergwerkplot_Control_Plot$Infected
standarderrorfromdata=max(SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)
MeanRandom=rnorm(N,mean=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Pfizer"],sd=sqrt(standarderrorfromdata^2))
ModelParamtemp=rmvnorm(N,mean=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)),sigma=Cov)
SDrandom=rnorm(N,mean=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],sd=SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1])

TableBootstrapEstimates_Plot<-data.frame("index"=rep(c(1:N),each=nrow(Bergwerkplot_Control_Plot)),
                                         "PeriInfectionNeutTiter_NormControl"=rep(Bergwerkplot_Control_Plot$PeriInfectionNeutTiter_NormControl,N),
                                         "Estimate"=NA,
                                         "Estimate_noerrorefficacy"=NA
)

for (i in 1:N) {
  tempdraw<-sample(Bergwerkplot_Control_Plot$PeriInfectionNeutTiter_NormControl,
                   sum(Bergwerkplot_Control_Plot$Infected),
                   replace=TRUE,
                   prob=binconf(Bergwerkplot_Control_Plot$Infected,sum(Bergwerkplot_Control_Plot$Infected),0.05)[,3])
  tempdraw_control<-sample(Bergwerkplot_Control_Plot$PeriInfectionNeutTiter_NormControl,
                           sum(Bergwerkplot_Control_Plot$Control),
                           replace=TRUE,
                           prob=binconf(Bergwerkplot_Control_Plot$Control,sum(Bergwerkplot_Control_Plot$Control),0.05)[,3])
  tempfreq<-tabulate(match(tempdraw,Bergwerkplot_Control_Plot$PeriInfectionNeutTiter_NormControl))/length(tempdraw)
  tempfreq_control<-tabulate(match(tempdraw_control,Bergwerkplot_Control_Plot$PeriInfectionNeutTiter_NormControl))/length(tempdraw_control)
  tempEstimates_noerrorinefficacy=(1-Efficacy_Chodick)*(tempfreq/tempfreq_control)
  TableBootstrapEstimates_Plot$Estimate_noerrorefficacy[TableBootstrapEstimates_Plot$index==i]=tempEstimates_noerrorinefficacy
}


Bergwerkplot_Control_Plot$FractionInfected_Upper_noerroreff=aggregate(Estimate_noerrorefficacy~PeriInfectionNeutTiter_NormControl,data=TableBootstrapEstimates_Plot,FUN=function(x){quantile(x,c(0.975))})[,2]
Bergwerkplot_Control_Plot$FractionInfected_Lower_noerroreff=aggregate(Estimate_noerrorefficacy~PeriInfectionNeutTiter_NormControl,data=TableBootstrapEstimates_Plot,FUN=function(x){quantile(x,c(0.025))})[,2]

TableBootstrapEstimates_NaN0_Plot=TableBootstrapEstimates_Plot
TableBootstrapEstimates_NaN0_Plot$Estimate_noerrorefficacy[is.nan(TableBootstrapEstimates_NaN0_Plot$Estimate_noerrorefficacy)]=0

Bergwerkplot_Control_Plot$FractionInfected_Upper_NaN0_noerroreff=aggregate(Estimate_noerrorefficacy~PeriInfectionNeutTiter_NormControl,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.975))})[,2]
Bergwerkplot_Control_Plot$FractionInfected_Lower_NaN0_noerroreff=aggregate(Estimate_noerrorefficacy~PeriInfectionNeutTiter_NormControl,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.025))})[,2]
Bergwerkplot_Control_Plot$FractionInfected_Lower_Combined_noerroreff=pmin(Bergwerkplot_Control_Plot$FractionInfected_Lower_noerroreff,Bergwerkplot_Control_Plot$FractionInfected_Lower_NaN0_noerroreff)
Bergwerkplot_Control_Plot$FractionInfected_Upper_Combined_noerroreff=pmax(Bergwerkplot_Control_Plot$FractionInfected_Upper_noerroreff,Bergwerkplot_Control_Plot$FractionInfected_Upper_NaN0_noerroreff)

xticks=c(0.0078,0.0156,0.0313,0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
yticks=100*c(0,0.2,0.4,0.6,0.8,1)

Figure3A<-ggplot(data=EstimatedEfficacy,aes(x=(10^(NeutValuelog10New))))+
  
  geom_errorbar(data=Bergwerkplot_Control_Plot,aes(x=PeriInfectionNeutTiter_NormControl,ymin=100*(1-FractionInfected_Upper_Combined_noerroreff),ymax=100*(1-FractionInfected_Lower_Combined_noerroreff),alpha=countcolumn),width=0.1) +
  geom_point(data=Bergwerkplot_Control_Plot,aes(x=PeriInfectionNeutTiter_NormControl,y=100*(1-(FractionInfected)),alpha=countcolumn),shape=21,fill="dodgerblue") +
  geom_line(aes(y=(100*ProbRemainUninfected_var)),linetype=1,color="red") +
  geom_ribbon(aes(ymin=100*(Lower),ymax=100*Upper),color=NA,fill="red",alpha=0.15) +
  scale_x_log10(breaks=xticks,labels=xticks) +
  scale_y_continuous(breaks=yticks) +
  coord_cartesian(ylim=c(0,105),xlim=c(0.0156/1.2,1.2*max(xticks)),expand=FALSE) +
  labs(y="Individual protection (%)",
       x="Neutralisation\n(fold of convalescence)",
       alpha="Number individuals") +
  scale_alpha(breaks=c(1,10,30,100,200)) +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-0.5, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        axis.text.x=element_text(angle=45,hjust=1),
        # legend.spacing.y = unit(0, 'cm'),
        # legend.key = element_rect(size = 1),
        # legend.key.size = unit(0.45, "cm"),
        # legend.margin = margin(t=0,b=0.7,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
  ) 


ggsave(plot=Figure3A,filename="Figure3A_FALSEDATA.pdf",height=4,width=5.5)









###FengData
# EfficacyFeng_overall=1-((55/4372)/(196/4194)) 
EfficacyFeng_Live=EstimatePopulationEfficacyBreakthroughModel(log10(FengMean_Live),SDFeng_Live,exp(ApproximateFeng_Live$estimate)) #Approx. from model  
EfficacyFeng_Pseudo=EstimatePopulationEfficacyBreakthroughModel(log10(FengMean_Pseudo_IU),SDFeng_Pseudo,exp(ApproximateFeng_Pseudo$estimate)) #Approx. from model  
FoldDilutions2=10*(2^c(0:11))
binningBy2fold_symptom<-cut(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="ID50"],FoldDilutions2)
ForNumerics_symptom<-split(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="ID50"],binningBy2fold_symptom)
FractionPer_symptom=sapply(ForNumerics_symptom,FUN=length)/length(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="ID50"])
binningBy2fold_control<-cut(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"],FoldDilutions2)
ForNumerics_control<-split(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"],binningBy2fold_control)
FractionPer_control=sapply(ForNumerics_control,FUN=length)/length(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"])

RawEfficacy_Feng=1-((1-EfficacyFeng_Pseudo)*(FractionPer_symptom/FractionPer_control))


binningBy2fold_symptom_NF50<-cut(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="NF50"],FoldDilutions2)
ForNumerics_symptom_NF50<-split(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="NF50"],binningBy2fold_symptom_NF50)
FractionPer_symptom_NF50=sapply(ForNumerics_symptom_NF50,FUN=length)/length(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="NF50"])
binningBy2fold_control_NF50<-cut(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"],FoldDilutions2)
ForNumerics_control_NF50<-split(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"],binningBy2fold_control_NF50)
FractionPer_control_NF50=sapply(ForNumerics_control_NF50,FUN=length)/length(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"])
RawEfficacy_Feng_NF50=1-((1-EfficacyFeng_Live)*(FractionPer_symptom_NF50/FractionPer_control_NF50))

FengPlot=data.frame("NeutLevel_Norm"=c(10^(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"]+(0.5*(log10(FoldDilutions2[1:(length(FoldDilutions2)-1)])+log10(FoldDilutions2[2:length(FoldDilutions2)]))))/FengMean_Pseudo,
                                       10^(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"]+(0.5*(log10(FoldDilutions2[1:(length(FoldDilutions2)-1)])+log10(FoldDilutions2[2:length(FoldDilutions2)]))))/FengMean_Live),
                    "Assay"=rep(c("Pseudo","Live"),each=length(FoldDilutions2)-1),
                    "Efficacy"=c(RawEfficacy_Feng,RawEfficacy_Feng_NF50))

FengPlot$Efficacy[FengPlot$Assay=="Pseudo"]=RawEfficacy_Feng
FengPlot$Efficacy[FengPlot$Assay=="Live"]=RawEfficacy_Feng_NF50
FengPlot$NumberIndividuals[FengPlot$Assay=="Pseudo"]=sapply(ForNumerics_control,FUN=length)+sapply(ForNumerics_symptom,FUN=length)
FengPlot$NumberIndividuals[FengPlot$Assay=="Live"]=sapply(ForNumerics_control_NF50,FUN=length)+sapply(ForNumerics_symptom_NF50,FUN=length)

NormFeng_Case_ID50<-nlm(function(p){-sum(log(dnorm(log10(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="ID50" & FengData$NeutLevel>FengLOD_Pseudo]),p[1],exp(p[2]))))-sum((FengData$Group=="Primary" & FengData$Assay=="ID50" & FengData$NeutLevel<FengLOD_Pseudo)*log(pnorm(log10(FengLOD_Pseudo),p[1],exp(p[2]))))},c(mean(log10(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="ID50"])),log(sd(log10(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="ID50"])))))
NormFeng_Control_ID50<-nlm(function(p){-sum(log(dnorm(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50" & FengData$NeutLevel>FengLOD_Pseudo]),p[1],exp(p[2]))))-sum((FengData$Group=="NAAT Negative" & FengData$Assay=="ID50" & FengData$NeutLevel<FengLOD_Pseudo)*log(pnorm(log10(FengLOD_Pseudo),p[1],exp(p[2]))))},c(mean(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"])),log(sd(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"])))))
NormFeng_Case_NF50<-nlm(function(p){-sum(log(dnorm(log10(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="NF50" & FengData$NeutLevel>FengLOD_Live]),p[1],exp(p[2]))))-sum((FengData$Group=="Primary" & FengData$Assay=="NF50" & FengData$NeutLevel<FengLOD_Live)*log(pnorm(log10(FengLOD_Live),p[1],exp(p[2]))))},c(mean(log10(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="NF50"])),log(sd(log10(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="NF50"])))))
NormFeng_Control_NF50<-nlm(function(p){-sum(log(dnorm(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50" & FengData$NeutLevel>FengLOD_Live]),p[1],exp(p[2]))))-sum((FengData$Group=="NAAT Negative" & FengData$Assay=="NF50" & FengData$NeutLevel<FengLOD_Live)*log(pnorm(log10(FengLOD_Live),p[1],exp(p[2]))))},c(mean(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"])),log(sd(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"])))))


TableBootstrapEstimates_Plot<-data.frame("index"=rep(c(1:N),each=nrow(FengPlot)),
                                         "NeutLevel"=rep(FengPlot$NeutLevel_Norm,N),
                                         "Assay"=rep(FengPlot$Assay,N),
                                         "Estimate"=NA
)


for (i in 1:N) {
  tempdraw_log<-rnorm(length(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="ID50"]),NormFeng_Case_ID50$estimate[1],exp(NormFeng_Case_ID50$estimate[2]))
  tempdraw_log_NF50<-rnorm(length(FengData$NeutLevel[FengData$Group=="Primary" & FengData$Assay=="NF50"]),NormFeng_Case_NF50$estimate[1],exp(NormFeng_Case_NF50$estimate[2]))
  tempdraw_control_log<-rnorm(length(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"]),NormFeng_Control_ID50$estimate[1],exp(NormFeng_Control_ID50$estimate[2]))
  tempdraw_control_log_NF50<-rnorm(length(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"]),NormFeng_Control_NF50$estimate[1],exp(NormFeng_Control_NF50$estimate[2]))
  tempdraw=10^tempdraw_log
  tempdraw_NF50=10^tempdraw_log_NF50
  tempdraw_control=10^tempdraw_control_log
  tempdraw_control_NF50=10^tempdraw_control_log_NF50
  tempsort<-split(tempdraw,cut(tempdraw,FoldDilutions2))
  tempsort_NF50<-split(tempdraw_NF50,cut(tempdraw_NF50,FoldDilutions2))
  tempfreq<-sapply(tempsort,FUN=length)/length(tempdraw)
  tempfreq_NF50<-sapply(tempsort_NF50,FUN=length)/length(tempdraw_NF50)
  tempsort_control<-split(tempdraw_control,cut(tempdraw_control,FoldDilutions2))
  tempsort_control_NF50<-split(tempdraw_control_NF50,cut(tempdraw_control_NF50,FoldDilutions2))
  tempfreq_control<-sapply(tempsort_control,FUN=length)/length(tempdraw_control)
  tempfreq_control_NF50<-sapply(tempsort_control_NF50,FUN=length)/length(tempdraw_control_NF50)
  tempEstimates=(1-EfficacyFeng_Pseudo)*(tempfreq/tempfreq_control)
  tempEstimates_NF50=(1-EfficacyFeng_Live)*(tempfreq_NF50/tempfreq_control_NF50)
  TableBootstrapEstimates_Plot$Estimate[TableBootstrapEstimates_Plot$index==i]=c(tempEstimates,tempEstimates_NF50)
}

temp<-aggregate(Estimate~NeutLevel+Assay,data=TableBootstrapEstimates_Plot,FUN=function(x){quantile(x,c(0.975))})
temp$Lower<-1-temp$Estimate
colnames(temp)[colnames(temp)=="NeutLevel"]<-"NeutLevel_Norm"
FengPlot<-join(FengPlot,temp[,c("Assay","NeutLevel_Norm","Lower")],by=c("Assay","NeutLevel_Norm"))
temp<-aggregate(Estimate~NeutLevel+Assay,data=TableBootstrapEstimates_Plot,FUN=function(x){quantile(x,c(0.025))})
temp$Upper<-1-temp$Estimate
colnames(temp)[colnames(temp)=="NeutLevel"]<-"NeutLevel_Norm"
FengPlot<-join(FengPlot,temp[,c("Assay","NeutLevel_Norm","Upper")],by=c("Assay","NeutLevel_Norm"))




TableBootstrapEstimates_NaN0_Plot=TableBootstrapEstimates_Plot
TableBootstrapEstimates_NaN0_Plot$Estimate[is.nan(TableBootstrapEstimates_NaN0_Plot$Estimate)]=1


temp<-aggregate(Estimate~NeutLevel+Assay,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.975))})
temp$Lower_NaN0<-1-temp$Estimate
colnames(temp)[colnames(temp)=="NeutLevel"]<-"NeutLevel_Norm"
FengPlot<-join(FengPlot,temp[,c("Assay","NeutLevel_Norm","Lower_NaN0")],by=c("Assay","NeutLevel_Norm"))
temp<-aggregate(Estimate~NeutLevel+Assay,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.025))})
temp$Upper_NaN0<-1-temp$Estimate
colnames(temp)[colnames(temp)=="NeutLevel"]<-"NeutLevel_Norm"
FengPlot<-join(FengPlot,temp[,c("Assay","NeutLevel_Norm","Upper_NaN0")],by=c("Assay","NeutLevel_Norm"))

FengPlot$Lower_NaN0=1-aggregate(Estimate~NeutLevel+Assay,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.975))})[,3]
FengPlot$Upper_NaN0=1-aggregate(Estimate~NeutLevel+Assay,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.025))})[,3]

FengPlot$Lower_Combined=pmin(FengPlot$Lower_NaN0,FengPlot$Lower)
FengPlot$Upper_Combined=pmax(FengPlot$Upper_NaN0,FengPlot$Upper)

xticks=c(0.0078,0.0156,0.0313,0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
yticks=100*c(0,0.2,0.4,0.6,0.8,1)

PopulationEfficacyEstimates$Assay[PopulationEfficacyEstimates$Study=="Feng Live"]="Live"
PopulationEfficacyEstimates$Assay[PopulationEfficacyEstimates$Study=="Feng Pseudo"]="Pseudo"


Figure3CD<-ggplot(data=EstimatedEfficacy,aes(x=(10^(NeutValuelog10New))))+
  geom_errorbar(data=FengPlot[FengPlot$NumberIndividuals>0,],aes(x=NeutLevel_Norm,ymin=100*Lower_Combined,ymax=100*Upper_Combined,alpha=NumberIndividuals),width=0.1) +
  geom_point(data=FengPlot[FengPlot$NumberIndividuals>0,],aes(x=NeutLevel_Norm,y=100*Efficacy,alpha=NumberIndividuals),shape=21,fill="dodgerblue") +
  geom_line(aes(y=(100*ProbRemainUninfected_var)),linetype=1,color="red") +
  geom_ribbon(aes(ymin=100*(Lower),ymax=100*Upper),color=NA,fill="red",alpha=0.15) +
  geom_line(data=PopulationEfficacyEstimates[grepl("Feng",PopulationEfficacyEstimates$Study) & PopulationEfficacyEstimates$Pop_Indv=="Individual",], aes(x=10^log10NeutLevel,y=100*Efficacy),alpha=1) + 
  scale_x_log10(breaks=xticks,labels=xticks) +
  scale_y_continuous(breaks=yticks) +
  coord_cartesian(ylim=c(0,105),xlim=c(0.0156/1.2,1.2*max(xticks)),expand=FALSE) +
  labs(y="Individual protection (%)",
       x="Neutralisation\n(fold of convalescence)",
       alpha="Number individuals") +
  scale_alpha(breaks=c(1,10,30,100,200)) +
  facet_rep_wrap(~Assay,ncol=2) +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # panel.spacing.y = unit(-0.5, "lines"),
        # panel.spacing.x = unit(-1, "lines"),
        axis.text.x=element_text(angle=45,hjust=1),
        # legend.spacing.y = unit(0, 'cm'),
        # legend.key = element_rect(size = 1),
        # legend.key.size = unit(0.45, "cm"),
        # legend.margin = margin(t=0,b=0.7,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
  ) 


ggsave(plot=Figure3CD,filename="Figure3CD.pdf",height=4,width=8.8)






###Moderna Extracted
EfficacyModerna=0.928 #Figure 4 panel c of Gilbert et al.
FoldDilutions2=(2^c(0:14))
GilbertCases=GilbertcID50Raw$cases
GilbertCases=GilbertCases[!is.na(GilbertCases)]
GilbertControl=GilbertcID50Raw$Control
binningBy2fold_symptom<-cut(GilbertCases,FoldDilutions2)
ForNumerics_symptom<-split(GilbertCases,binningBy2fold_symptom)
FractionPer_symptom=sapply(ForNumerics_symptom,FUN=length)/length(GilbertCases)
binningBy2fold_control<-cut(GilbertControl,FoldDilutions2)
ForNumerics_control<-split(GilbertControl,binningBy2fold_control)
FractionPer_control=sapply(ForNumerics_control,FUN=length)/length(GilbertControl)
RawEfficacy_Gilbert=1-((1-EfficacyModerna)*(FractionPer_symptom/FractionPer_control))


xvector_foldConv=(10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])*(FoldDilutions2[1:(length(FoldDilutions2)-1)])/ModernaMean_IC50
GilbertRawDataEfficacy=data.frame("Titre"=xvector_foldConv,
                                  "Efficacy"=RawEfficacy_Gilbert,
                                  "NumberIndividuals"=FractionPer_symptom*length(GilbertCases)+FractionPer_control*length(GilbertControl))









TableBootstrapEstimates_Plot<-data.frame("index"=rep(c(1:N),each=length(xvector_foldConv)),
                                         "NeutLevel"=rep(xvector_foldConv,N),
                                         "Estimate"=NA
)


for (i in 1:N) {
  tempdraw_log<-rnorm(length(GilbertCases),NormGilbert_Case_ID50$estimate[1],exp(NormGilbert_Case_ID50$estimate[2]))
  tempdraw_log[tempdraw_log<log10(LODGilbertcID50)]=log10(LODGilbertcID50/2)
  tempdraw_control_log<-rnorm(length(GilbertControl),NormGilbert_Control_ID50$estimate[1],exp(NormGilbert_Control_ID50$estimate[2]))
  tempdraw_control_log[tempdraw_control_log<log10(LODGilbertcID50)]=log10(LODGilbertcID50/2)
  tempdraw=10^tempdraw_log
  tempdraw_control=10^tempdraw_control_log
  tempsort<-split(tempdraw,cut(tempdraw,FoldDilutions2))
  tempfreq<-sapply(tempsort,FUN=length)/length(tempdraw)
  tempsort_control<-split(tempdraw_control,cut(tempdraw_control,FoldDilutions2))
  tempfreq_control<-sapply(tempsort_control,FUN=length)/length(tempdraw_control)
  tempEstimates=(1-EfficacyModerna)*(tempfreq/tempfreq_control)
  TableBootstrapEstimates_Plot$Estimate[TableBootstrapEstimates_Plot$index==i]=tempEstimates
}

GilbertRawDataEfficacy$Lower=1-aggregate(Estimate~NeutLevel,data=TableBootstrapEstimates_Plot,FUN=function(x){quantile(x,c(0.975))})[,2]
GilbertRawDataEfficacy$Upper=1-aggregate(Estimate~NeutLevel,data=TableBootstrapEstimates_Plot,FUN=function(x){quantile(x,c(0.025))})[,2]

TableBootstrapEstimates_NaN0_Plot=TableBootstrapEstimates_Plot
TableBootstrapEstimates_NaN0_Plot$Estimate[is.nan(TableBootstrapEstimates_NaN0_Plot$Estimate)]=1
TableBootstrapEstimates_NaN0_Plot$Estimate[is.nan(TableBootstrapEstimates_NaN0_Plot$Estimate)]=1

GilbertRawDataEfficacy$Lower_NaN0=1-aggregate(Estimate~NeutLevel,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.975))})[,2]
GilbertRawDataEfficacy$Upper_NaN0=1-aggregate(Estimate~NeutLevel,data=TableBootstrapEstimates_NaN0_Plot,FUN=function(x){quantile(x,c(0.025))})[,2]
GilbertRawDataEfficacy$CombinedLower=pmin(GilbertRawDataEfficacy$Lower_NaN0,GilbertRawDataEfficacy$Lower)
GilbertRawDataEfficacy$CombinedUpper=pmax(GilbertRawDataEfficacy$Upper_NaN0,GilbertRawDataEfficacy$Upper)


yticks=100*c(0,0.2,0.4,0.6,0.8,1)

Figure3B<-ggplot(data=EstimatedEfficacy,aes(x=(10^(NeutValuelog10New))))+
  geom_errorbar(data=GilbertRawDataEfficacy[GilbertRawDataEfficacy$NumberIndividuals>0,],aes(x=Titre,ymin=100*CombinedLower,ymax=100*CombinedUpper,alpha=NumberIndividuals),width=0.1) +
  geom_line(aes(y=100*(ProbRemainUninfected_var)),linetype=1,color="red") +
  geom_point(data=GilbertRawDataEfficacy[GilbertRawDataEfficacy$NumberIndividuals>0,],aes(x=Titre,y=100*Efficacy,alpha=NumberIndividuals),shape=21,fill="dodgerblue") +
  geom_line(data=PopulationEfficacyEstimates[PopulationEfficacyEstimates$Study=="Gilbert cID50" & PopulationEfficacyEstimates$Pop_Indv=="Individual",],aes(x=10^log10NeutLevel,y=100*Efficacy)) +
  geom_ribbon(aes(ymin=100*(Lower),ymax=100*Upper),color=NA,fill="red",alpha=0.15) +
  scale_x_log10(breaks=xticks,labels=xticks) +
  scale_y_continuous(breaks=yticks) +
  coord_cartesian(ylim=c(0,105),xlim=c(0.0156/1.2,1.2*max(xticks)),expand=FALSE) +
  labs(y="Individual protection (%)",
       x="Neutralisation\n(fold of convalescence)",
       alpha="Number individuals") +
  scale_alpha(breaks=c(1,10,30,100,200)) +
  # facet_wrap(~Assay,ncol=2) +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # panel.spacing.y = unit(-0.5, "lines"),
        # panel.spacing.x = unit(-1, "lines"),
        axis.text.x=element_text(angle=45,hjust=1),
        # legend.spacing.y = unit(0, 'cm'),
        # legend.key = element_rect(size = 1),
        # legend.key.size = unit(0.45, "cm"),
        # legend.margin = margin(t=0,b=0.7,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
  ) 


ggsave(plot=Figure3B,"Figure3B.pdf",height=4,width=5.5)

