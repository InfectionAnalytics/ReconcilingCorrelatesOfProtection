###Figure 2

##Plotting Figure 2


#Import data
FengGilbertCurves<-read.csv(paste0(currentdirectory,"/Data/FengGilbertCurves_Long.csv"),fileEncoding="UTF-8-BOM")


######Adjusted Curves (to fold of convalescent scale) for Published Models 
# by dividing by matching reported GMT for each vaccine and using Khoury et al. 
# fold-convalescent estimate for each vaccine
FengGilbertCurves$NeutAdj[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Pseudo"]<-(10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"])*FengGilbertCurves$Neut[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Pseudo"]/(FengMean_Pseudo_IU)
FengGilbertCurves$NeutAdj[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Live"]<-(10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"])*FengGilbertCurves$Neut[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Live"]/FengMean_Live
FengGilbertCurves$NeutAdj[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID50"]<-(10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])*FengGilbertCurves$Neut[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID50"]/ModernaMean_IC50
FengGilbertCurves$NeutAdj[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID80"]<-(10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])*FengGilbertCurves$Neut[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID80"]/ModernaMean_IC80

FengGilbertCurves$Neutlog10=log10(FengGilbertCurves$Neut)



#Plotting the scaled models from each study Figure 2
tempGilbertCIsLB<-FengGilbertCurves[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Item=="LB",]
tempGilbertCIsLB<-tempGilbertCIsLB[order(tempGilbertCIsLB$NeutAdj),]
tempGilbertCIsUB<-FengGilbertCurves[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Item=="UB",]
tempGilbertCIsUB<-tempGilbertCIsUB[order(tempGilbertCIsUB$NeutAdj,decreasing = TRUE),]
tempGilbertCIs<-rbind(tempGilbertCIsLB,tempGilbertCIsUB)
extrarow<-tempGilbertCIsUB[tempGilbertCIsUB$Assay=="cID50" & tempGilbertCIsUB$NeutAdj==min(tempGilbertCIsUB$NeutAdj[tempGilbertCIsUB$Assay=="cID50"] ),]
extrarow<-rbind(extrarow,tempGilbertCIsUB[tempGilbertCIsUB$Assay=="cID80" & tempGilbertCIsUB$NeutAdj==min(tempGilbertCIsUB$NeutAdj[tempGilbertCIsUB$Assay=="cID80"]),])
extrarow$Efficacy<-0
tempGilbertCIsLB<-rbind(extrarow,tempGilbertCIsLB)
tempGilbertCIs<-rbind(tempGilbertCIsLB,tempGilbertCIsUB)



### Shading 90th percentile of data.
###These IQR are extracted from figures and tables in papers.
yminvector=c(10,20,10,20)
ymaxvector=c(15,25,15,25)

ShadingDataRegion80=data.frame(Normalise=factor(rep(c("AstraZeneca","Moderna"),each=2),levels=c("Moderna","AstraZeneca")),
                               Study=c("Feng Live","Feng Pseudo","Gilbert cID50","Gilbert cID80"),
                               xmin=c((10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"])*(10^qnorm(0.1,log10(FengMean_Pseudo_IU),SDFeng_Pseudo))/FengMean_Pseudo_IU,
                                      (10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"])*(10^qnorm(0.1,log10(FengMean_Live),SDFeng_Live))/FengMean_Live,
                                      (10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])*(10^qnorm(0.1,log10(ModernaMean_IC50),NeutModernaControl_sd))/ModernaMean_IC50,
                                      (10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])*(10^qnorm(0.1,log10(ModernaMean_IC80),NeutModernaControlID80_sd))/ModernaMean_IC80),
                               xmax=c((10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"])*(10^qnorm(0.9,log10(FengMean_Pseudo_IU),SDFeng_Pseudo))/FengMean_Pseudo_IU,
                                      (10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="AstraZeneca"])*(10^qnorm(0.9,log10(FengMean_Live),SDFeng_Live))/FengMean_Live,
                                      (10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])*(10^qnorm(0.9,log10(ModernaMean_IC50),NeutModernaControl_sd))/ModernaMean_IC50,
                                      (10^SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"])*(10^qnorm(0.9,log10(ModernaMean_IC80),NeutModernaControlID80_sd))/ModernaMean_IC80),
                               ymin=yminvector,
                               ymax=ymaxvector
)
MinOfAssay_ShadingDataRegion80=aggregate(xmin~Normalise,data=ShadingDataRegion80,FUN=min)
MaxOfAssay_ShadingDataRegion80=aggregate(xmax~Normalise,data=ShadingDataRegion80,FUN=max)
ExtremesOf80_ForPlotting=join(MinOfAssay_ShadingDataRegion80,MaxOfAssay_ShadingDataRegion80)



###Plotting
colorlist=c("chocolate4","deepskyblue3","darkorange2")
yticks=c(40,50,60,70,80,90,100)
xticks=c(0.0078,0.0156,0.0313,0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
Figure2A<-ggplot(data=EstimatedEfficacy,aes(x=(10^(NeutValuelog10New)))) +
  geom_vline(xintercept=as.numeric(ExtremesOf80_ForPlotting[ExtremesOf80_ForPlotting$Normalise=="Moderna",c("xmin","xmax")]),color="gray85") +
  annotate(geom="segment",x=ExtremesOf80_ForPlotting$xmin[ExtremesOf80_ForPlotting$Normalise=="Moderna"],
           xend=ExtremesOf80_ForPlotting$xmax[ExtremesOf80_ForPlotting$Normalise=="Moderna"],
           y=50,
           yend=50,
           arrow = arrow(ends = "both", angle = 45, length = unit(.2,"cm")),
           color="gray85") +
  annotate(geom="text",
           x=exp(mean(log(as.numeric(ExtremesOf80_ForPlotting[ExtremesOf80_ForPlotting$Normalise=="Moderna",c("xmin","xmax")])))),
           y=52,
           label = "10-90th\npercentile") +
  geom_polygon(data=tempGilbertCIs,aes(x=NeutAdj,y=(Efficacy),fill=Assay),color=NA,alpha=0.15) +
  geom_line(data=FengGilbertCurves[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Item=="mean",],aes(x=NeutAdj,y=Efficacy,color=Assay)) +
  geom_ribbon(aes(ymin=100*(Lower),ymax=100*Upper),color=NA,fill="red",alpha=0.15) +
  geom_line(aes(y=100*(ProbRemainUninfected_var)),linetype=1,color="red") +
  scale_x_log10(breaks=xticks,labels=xticks) +
  scale_y_continuous(breaks=yticks) +
  scale_color_manual(values=colorlist) +
  scale_fill_manual(values=colorlist) +
  coord_cartesian(ylim=c(40,105),xlim=c(0.0625,1.2*32),expand=FALSE) +
  labs(y="Individual protection (%)",
       x="Neutralisation\n(fold of convalescence)") +
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


ggsave(plot=Figure2A,filename="Figure2A.pdf",height=3.4,width=4)



tempFengCIsLB<-FengGilbertCurves[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Item=="LB",]
tempFengCIsLB<-tempFengCIsLB[order(tempFengCIsLB$NeutAdj),]
tempFengCIsUB<-FengGilbertCurves[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Item=="UB",]
tempFengCIsUB<-tempFengCIsUB[order(tempFengCIsUB$NeutAdj,decreasing = TRUE),]
extrarow<-tempFengCIsUB[tempFengCIsUB$Assay=="Live" & tempFengCIsUB$NeutAdj==min(tempFengCIsUB$NeutAdj[tempFengCIsUB$Assay=="Live"] ),]
extrarow<-rbind(extrarow,tempFengCIsUB[tempFengCIsUB$Assay=="Pseudo" & tempFengCIsUB$NeutAdj==min(tempFengCIsUB$NeutAdj[tempFengCIsUB$Assay=="Pseudo"]),])
extrarow$Efficacy<-40
tempFengCIsLB<-rbind(extrarow,tempFengCIsLB)
tempFengCIs<-rbind(tempFengCIsLB,tempFengCIsUB)



colorlist=c("mediumpurple","darkgoldenrod")

xticks=c(0.0078,0.0156,0.0313,0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
Figure2B<-ggplot(data=EstimatedEfficacy,aes(x=(10^(NeutValuelog10New)))) +
  geom_vline(xintercept=as.numeric(ExtremesOf80_ForPlotting[ExtremesOf80_ForPlotting$Normalise=="AstraZeneca",c("xmin","xmax")]),color="gray85") +
  annotate(geom="segment",x=ExtremesOf80_ForPlotting$xmin[ExtremesOf80_ForPlotting$Normalise=="AstraZeneca"],
           xend=ExtremesOf80_ForPlotting$xmax[ExtremesOf80_ForPlotting$Normalise=="AstraZeneca"],
           y=50,
           yend=50,
           arrow = arrow(ends = "both", angle = 45, length = unit(.2,"cm")),
           color="gray85") +
  annotate(geom="text",
           x=exp(mean(log(as.numeric(ExtremesOf80_ForPlotting[ExtremesOf80_ForPlotting$Normalise=="AstraZeneca",c("xmin","xmax")])))),
           y=52,
           label = "10-90th\npercentile") +
  geom_polygon(data=tempFengCIs,aes(x=NeutAdj,y=(Efficacy),fill=Assay),color=NA,alpha=0.15) +
  geom_line(data=FengGilbertCurves[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Item=="mean",],aes(x=NeutAdj,y=Efficacy,color=Assay)) +
  geom_ribbon(aes(ymin=100*(Lower),ymax=100*Upper),color=NA,fill="red",alpha=0.15) +
  geom_line(aes(y=100*(ProbRemainUninfected_var)),linetype=1,color="red") +
  scale_x_log10(breaks=xticks,labels=xticks) +
  scale_y_continuous(breaks=yticks) +
  scale_color_manual(values=colorlist) +
  scale_fill_manual(values=colorlist) +
  coord_cartesian(ylim=c(40,105),xlim=c(0.0625,1.2*32),expand=FALSE) +
  labs(y="Individual protection (%)",
       x="Neutralisation\n(fold of convalescence)") +
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


ggsave(plot=Figure2B,filename="Figure2B.pdf",height=3.4,width=4.3)


