#### Estimate population efficacy using Feng and Gilbert Models

## Make figure S1#######

###Estimating population efficacy for a group of vaccinated individuals
# from Feng and Gilbert efficacy curves

EstimatePopulationEfficacyBreakthroughModel<-function(mu_titre,sig_titre,p){
  Output<-NULL
  for (i in 1:length(mu_titre)){
    if (is.na(mu_titre[i])) {
      Output[i]=NA
    } else {
      Step=sig_titre*0.001
      IntegralVector=seq(mu_titre[i]-5*sig_titre,mu_titre[i]+5*sig_titre,by=Step)
      Output[i]=sum(Step*BreakthroughEfficacyFunction(IntegralVector[1:(length(IntegralVector)-1)],p[1],p[2],p[3])*dnorm(IntegralVector[1:(length(IntegralVector)-1)],mu_titre[i],sig_titre))
    }
  }
  Output
}




PopulationEfficacyEstimates<-data.frame("Study"=c(rep("Khoury",2*length(NeutValuelog10New)),rep("Gilbert cID50",2*length(NeutValuelog10New)),rep("Gilbert cID80",2*length(NeutValuelog10New))),
                                        "Pop_Indv"=rep(c("Population","Individual","Population","Individual","Population","Individual"),each=length(NeutValuelog10New)),
                                        "log10NeutLevel"=rep(NeutValuelog10New,6),
                                        "Efficacy"=NA,
                                        "Efficacy_LB"=NA,
                                        "Efficacy_UB"=NA)


PopulationEfficacyEstimates$Efficacy[PopulationEfficacyEstimates$Study=="Khoury" & PopulationEfficacyEstimates$Pop_Indv=="Population"]=LogisticModel_PercentUninfected(NeutValuelog10New,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])
PopulationEfficacyEstimates$Efficacy[PopulationEfficacyEstimates$Study=="Gilbert cID50" & PopulationEfficacyEstimates$Pop_Indv=="Population"]=EstimatePopulationEfficacyBreakthroughModel(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"]+log10(ModernaMean_IC50),NeutModernaControl_sd,exp(ApproximateEfficacyFunctionFit$estimate))
PopulationEfficacyEstimates$Efficacy[PopulationEfficacyEstimates$Study=="Gilbert cID80" & PopulationEfficacyEstimates$Pop_Indv=="Population"]=EstimatePopulationEfficacyBreakthroughModel(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"]+log10(ModernaMean_IC80),NeutModernaControlID80_sd,exp(ApproximateEfficacyFunctionFit_cID80$estimate))
PopulationEfficacyEstimates$Efficacy[PopulationEfficacyEstimates$Study=="Khoury" & PopulationEfficacyEstimates$Pop_Indv=="Individual"]=ProbRemainUninfected(NeutValuelog10New,tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])
PopulationEfficacyEstimates$Efficacy[PopulationEfficacyEstimates$Study=="Gilbert cID50" & PopulationEfficacyEstimates$Pop_Indv=="Individual"]=BreakthroughEfficacyFunction(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"]+log10(ModernaMean_IC50),exp(ApproximateEfficacyFunctionFit$estimate[1]),exp(ApproximateEfficacyFunctionFit$estimate[2]),exp(ApproximateEfficacyFunctionFit$estimate[3]))
PopulationEfficacyEstimates$Efficacy[PopulationEfficacyEstimates$Study=="Gilbert cID80" & PopulationEfficacyEstimates$Pop_Indv=="Individual"]=BreakthroughEfficacyFunction(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"]+log10(ModernaMean_IC80),exp(ApproximateEfficacyFunctionFit_cID80$estimate[1]),exp(ApproximateEfficacyFunctionFit_cID80$estimate[2]),exp(ApproximateEfficacyFunctionFit_cID80$estimate[3]))
PopulationEfficacyEstimates$Normalise="Moderna"


PopulationEfficacyEstimates_temp<-data.frame("Study"=c(rep("Khoury",2*length(NeutValuelog10New)),rep("Feng Pseudo",2*length(NeutValuelog10New)),rep("Feng Live",2*length(NeutValuelog10New))),
                                             "Pop_Indv"=rep(c("Population","Individual","Population","Individual","Population","Individual"),each=length(NeutValuelog10New)),
                                             "log10NeutLevel"=rep(NeutValuelog10New,6),
                                             "Efficacy"=NA,
                                             "Efficacy_LB"=NA,
                                             "Efficacy_UB"=NA)

###The protection curve for AZ was plotted (and extacted) on the IU scale. So we must normalise to this.
PopulationEfficacyEstimates_temp$Efficacy[PopulationEfficacyEstimates_temp$Study=="Khoury" & PopulationEfficacyEstimates_temp$Pop_Indv=="Population"]=LogisticModel_PercentUninfected(NeutValuelog10New,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])
PopulationEfficacyEstimates_temp$Efficacy[PopulationEfficacyEstimates_temp$Study=="Feng Pseudo" & PopulationEfficacyEstimates_temp$Pop_Indv=="Population"]=EstimatePopulationEfficacyBreakthroughModel(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"]+log10(FengMean_Pseudo_IU),SDFeng_Pseudo,exp(ApproximateFeng_Pseudo$estimate))
PopulationEfficacyEstimates_temp$Efficacy[PopulationEfficacyEstimates_temp$Study=="Feng Live" & PopulationEfficacyEstimates_temp$Pop_Indv=="Population"]=EstimatePopulationEfficacyBreakthroughModel(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"]+log10(FengMean_Live),SDFeng_Live,exp(ApproximateFeng_Live$estimate))
PopulationEfficacyEstimates_temp$Efficacy[PopulationEfficacyEstimates_temp$Study=="Khoury" & PopulationEfficacyEstimates_temp$Pop_Indv=="Individual"]=ProbRemainUninfected(NeutValuelog10New,tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])
PopulationEfficacyEstimates_temp$Efficacy[PopulationEfficacyEstimates_temp$Study=="Feng Pseudo" & PopulationEfficacyEstimates_temp$Pop_Indv=="Individual"]=BreakthroughEfficacyFunction(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"]+log10(FengMean_Pseudo_IU),exp(ApproximateFeng_Pseudo$estimate[1]),exp(ApproximateFeng_Pseudo$estimate[2]),exp(ApproximateFeng_Pseudo$estimate[3]))
PopulationEfficacyEstimates_temp$Efficacy[PopulationEfficacyEstimates_temp$Study=="Feng Live" & PopulationEfficacyEstimates_temp$Pop_Indv=="Individual"]=BreakthroughEfficacyFunction(NeutValuelog10New-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"]+log10(FengMean_Live),exp(ApproximateFeng_Live$estimate[1]),exp(ApproximateFeng_Live$estimate[2]),exp(ApproximateFeng_Live$estimate[3]))
PopulationEfficacyEstimates_temp$Normalise="AstraZeneca"


PopulationEfficacyEstimates=rbind(PopulationEfficacyEstimates,PopulationEfficacyEstimates_temp)





####Table for plotting
PlotVaccine=SummaryTable_Efficacy_NeutRatio_SD_SEM
PlotVaccine[nrow(PlotVaccine)+1,]<-NA
PlotVaccine$Study[nrow(PlotVaccine)]<-"Covaxin"
PlotVaccine$TechnicalName[nrow(PlotVaccine)]<-"Covaxin"
CovaxinNeutMean<-127.6773881
CovaxinNeutConv<-161.1587484
PlotVaccine$NeutRatio_Reported[nrow(PlotVaccine)]<-log10(CovaxinNeutMean/CovaxinNeutConv)
PlotVaccine$Efficacy[nrow(PlotVaccine)]<-0.806
PlotVaccine$Pop_Indv="Population"

###Adding in AstraZeneca data on efficacy of different dosing intervals from Metaanalysis.
FengByDoseInterval=data.frame("Study"="AstraZeneca",
                              "Interval"=c("<6","6-8","9-11",">12"),
                              "NeutRatioAbove4w"=c(1,1.43,1.55,1.85),
                              "Efficacy"=c(55.1,59.7,72.2,80))



# ShadingDataRegion80$log10NeutLevel=NA
# ShadingDataRegion80$Efficacy=NA
# ShadingDataRegion70$log10NeutLevel=NA
# ShadingDataRegion70$Efficacy=NA
# ShadingDataRegion60$log10NeutLevel=NA
# ShadingDataRegion60$Efficacy=NA
# ShadingDataRegion50$log10NeutLevel=NA
# ShadingDataRegion50$Efficacy=NA

xticks=c(0.0625,0.125,0.25,0.5,1,2,4,8,16)
PopulationEfficacyEstimates$Normalise<-factor(PopulationEfficacyEstimates$Normalise,levels=c("Moderna","AstraZeneca"))


FacetLabel=c("mRNA-1273 (Gilbert et al.)","ChAdOx1 (Feng et al.)")
names(FacetLabel)=c("Moderna","AstraZeneca")

colorlist=c("mediumpurple","darkgoldenrod","chocolate4","deepskyblue3","red","darkorange2")
linetypelist=c("solid","solid","solid","solid","longdash")
##The 12w spacing group had 80% efficacy (pulled twelce week out of the meta analysis paper)
FigureS1<-ggplot(data=PopulationEfficacyEstimates[PopulationEfficacyEstimates$Pop_Indv=="Population",],aes(x=10^(log10NeutLevel),y=100*Efficacy)) +
  geom_rect(data=ExtremesOf80_ForPlotting,mapping=aes(xmin=xmin,xmax=xmax,ymin=0,ymax=100),color=NA,alpha=0.1,inherit.aes = FALSE) +
  geom_segment(data=ExtremesOf80_ForPlotting,
               aes(x=xmin,xend=xmax,y=20,yend=20),
               arrow = arrow(ends = "both", angle = 45, length = unit(.2,"cm")),
               color="gray50",
               inherit.aes = FALSE) +
  geom_text(data=data.frame("Normalise"=factor(c("AstraZeneca","Moderna"),levels=c("Moderna","AstraZeneca")),"Adjustment"=c(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"],SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])),
            aes(x=10^Adjustment,y=25),
            label="10-90th\npercentiles",
            color="gray50",
            size=3,
            lineheight=0.5,
            inherit.aes = FALSE) +
  geom_vline(data=data.frame("Normalise"=factor(c("AstraZeneca","Moderna"),levels=c("Moderna","AstraZeneca")),"Adjustment"=c(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"],SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])),aes(xintercept=10^Adjustment),col="grey",size=0.1) +
  geom_line(aes(color=Study,linetype=Study)) +
  geom_point(data=PlotVaccine,aes(x=10^(NeutRatio_Reported),y=100*Efficacy,shape=TechnicalName,size=TechnicalName),fill="black") +
  geom_point(data=FengByDoseInterval,aes(x=(NeutRatioAbove4w)*(10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"]),y=Efficacy),shape=21,col="blue") +
  geom_segment(data=ExtremesOf80_ForPlotting,
               aes(x=xmin,xend=xmax,y=20,yend=20),
               arrow = arrow(ends = "both", angle = 45, length = unit(.2,"cm")),
               color="gray85",
               inherit.aes = FALSE) +
  scale_x_log10(limits=c(0.0625,20),breaks=xticks,labels=xticks) +
  scale_y_continuous(breaks=c(0,20,40,60,80,100)) +
  coord_cartesian(ylim=c(0,105),expand=FALSE) +
  scale_shape_manual(values=list(23,22,21,3,4,8,-as.hexmode("002A"),25,24),
                     guide=guide_legend(order=1,
                                        override.aes = list(fill=NA,size=c(rep(1.5,6),4,rep(1.5,2))))) +
  scale_size_manual(values=c(rep(1.5,6),4,rep(1.5,2)),guide="none") +
  scale_color_manual(values=colorlist) +
  scale_linetype_manual(values=linetypelist,guide="none") +
  facet_rep_wrap(~Normalise,
                 ncol=2,
                 labeller=labeller(Normalise=FacetLabel)
  ) +
  labs(x="Neutralisation titre\n(fold of convalescence)",
       y="Efficacy (%)",
       shape="Vaccine/Serum",
       color="Model") +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        axis.text.x = element_text(angle = 45,hjust=0.5,vjust=0.75),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        # panel.spacing.y = unit(-0.5, "lines"),
        # panel.spacing.x = unit(-1, "lines"),
        legend.spacing.y = unit(0, 'cm'),
        # legend.key = element_rect(size = 1),
        # legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(t=1,b=-0.5,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
  )


ggsave(plot=FigureS1,filename="FigureS1.pdf",height=4,width=8.5)


