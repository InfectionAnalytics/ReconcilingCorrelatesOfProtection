###Estimating the non-inferiority margins for different vaccines.


ReferenceVaccines = c("Pfizer", "AstraZeneca", "Moderna")

NonInferiorityMargin = data.frame("Reference"=SummaryTable_Efficacy_NeutRatio_SD_SEM$Study[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study %in% ReferenceVaccines],
                                  "NeutRatio_Reported"=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study %in% ReferenceVaccines])


NeutVector=seq(log10(0.125),log10(8),by=0.1)
ModelMean<-LogisticModel_PercentUninfected(NeutVector,
                                SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],
                                tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])
TableForReferenceFigure=data.frame("Reference"=rep(ReferenceVaccines,each=length(NeutVector)),
                                   "Neut"=rep(NeutVector,length(ReferenceVaccines)),
                                   "Efficacy"=rep(ModelMean,length(ReferenceVaccines)),
                                   "LowerBound"=NA)


#################################################### Adding 95% Confidence Intervals #########################################################

###Using Bootstrap to estimate CIs for the ribbon
N = 50000
lowerquant <- 0.05
ModelParamtemp = rmvnorm(N, mean = c(
  tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate, 2)
), sigma = Cov)
SDrandom = rnorm(N, mean = SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1], sd =
                   SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1])


EstimateLowerBound = function(Neut, Reference) {
  # SE_difference=sqrt((SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]/sqrt(StudySizeTemp))^2+(SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]/sqrt(StudySizeTemp))^2)
  SE_Ref = SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM[as.character(SummaryTable_Efficacy_NeutRatio_SD_SEM$Study) ==
                                                        Reference]
  LowerBound = NULL
  for (i in 1:length(Neut)) {
    tempNeut = rnorm(N, mean = Neut[i], sd = ((SE_Ref)))
    tempEvaluateFunction <-
      LogisticModel_PercentUninfected(tempNeut, SDrandom, ModelParamtemp[, 1], ModelParamtemp[, 2])
    
    LowerBound[i] = quantile(tempEvaluateFunction, lowerquant)
  }
  LowerBound
}



# cores = detectCores()
# coresuse = cores[1] - 1 #not to overload your computer
# cl <- makeCluster(coresuse)
# registerDoParallel(cl)


start_time <- Sys.time()

output <- foreach(i = 1:nrow(TableForReferenceFigure),
                      .combine = rbind) %do% {
                        Reference = as.character(TableForReferenceFigure$Reference[i])
                        Neut=TableForReferenceFigure$Neut[i]
                        
                        OutputVect <- EstimateLowerBound(Neut, Reference)
                        
                        tempMatrix = OutputVect
                        
                        tempMatrix
                        
                      }
end_time <- Sys.time()
end_time - start_time

#stop cluster
# stopCluster(cl)

TableForReferenceFigure$LowerBound=output
TableForReferenceFigure$TechnicalName_Ref[TableForReferenceFigure$Reference=="Pfizer"]="BNT162b2"
TableForReferenceFigure$TechnicalName_Ref[TableForReferenceFigure$Reference=="Moderna"]="mRNA-1273"
TableForReferenceFigure$TechnicalName_Ref[TableForReferenceFigure$Reference=="AstraZeneca"]="ChAdOx1 nCoV-19"




######

LevelGiveLB80 = function(NeutVector, Reference) {
  LowerBound <- EstimateLowerBound(NeutVector, Reference)
  indexListmin = which(abs(LowerBound - 0.8) < 0.01 &
                         (LowerBound - 0.8) <= 0)
  indexListmax = which(abs(LowerBound - 0.8) < 0.01 &
                         (LowerBound - 0.8) >= 0)
  minindex = min(indexListmin)
  maxindex = max(indexListmax)
  interpy = LowerBound[minindex:maxindex] - 0.8
  interpx = NeutVector[minindex:maxindex]
  tempfit <- lm(interpy ~ interpx)
  Neutfor80LB <- -tempfit$coef[1] / tempfit$coef[2]
  Neutfor80LB
}


##Parralel for speed.
# cores = detectCores()
# coresuse = min(nrow(NonInferiorityMargin), cores[1] - 1) 
# cl <- makeCluster(coresuse)
# registerDoParallel(cl)
# 
# 
# start_time <- Sys.time()



estimateNeutfor80 = log10(1)


fullMatrix <- foreach(i = 1:nrow(NonInferiorityMargin),
                      .combine = rbind) %do% {
                        Reference = as.character(NonInferiorityMargin$Reference[i])
                        NeutVector = seq((estimateNeutfor80), log10(2) + (estimateNeutfor80), by =
                                           0.01)
                        
                        OutputVect <- LevelGiveLB80(NeutVector, Reference)
                        
                        tempMatrix = OutputVect
                        
                        tempMatrix
                        
                      }
# end_time <- Sys.time()
# end_time - start_time

#stop cluster
# stopCluster(cl)

NonInferiorityMargin$MarginForAbove80 = fullMatrix[, 1]

NonInferiorityMargin$NonInf=NonInferiorityMargin$NeutRatio_Reported>NonInferiorityMargin$MarginForAbove80

NonInferiorityMargin$Margin=10^(NonInferiorityMargin$MarginForAbove80-NonInferiorityMargin$NeutRatio_Reported)

checkn=10
NonInferiorityMargin$ConfirmLowerBound=NA
for (i in 1:nrow(NonInferiorityMargin)){
  tempresult<-NULL
  # cores = detectCores()
  # coresuse = cores[1] - 1 #not to overload your computer
  # cl <- makeCluster(coresuse)
  # registerDoParallel(cl)
  outputtemp<-foreach(j=1:checkn,.combine=rbind) %do% {
    tempresult<-EstimateLowerBound(NonInferiorityMargin$MarginForAbove80[i],NonInferiorityMargin$Reference[i])
    tempresult
  }
  # stopCluster(cl)
  NonInferiorityMargin$ConfirmLowerBound[i]<-mean(outputtemp)
}



xticks=2^seq(log2(0.125),log2(8),by=1)
yticks=c(0,20,40,60,80,100)

PointsTable<-SummaryTable_Efficacy_NeutRatio_SD_SEM[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study %in% ReferenceVaccines,]
PointsTable$TechnicalName_Ref<-PointsTable$TechnicalName
NonInferiorityMargin$TechnicalName_Ref<-c("BNT162b2","mRNA-1273","ChAdOx1 nCoV-19")

tempTable=NonInferiorityMargin[,c("Reference","MarginForAbove80","TechnicalName_Ref","ConfirmLowerBound")]
colnames(tempTable)[2]="Neut"
colnames(tempTable)[4]="LowerBound"
tempTable$Efficacy=LogisticModel_PercentUninfected(tempTable$Neut,
                                                   SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],
                                                   tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                                   tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])


tempTablev2<-tempTable
tempTablev2$LowerBound[tempTablev2$LowerBound<0.8]=0.8000
tempTablev2$LowerBound[tempTablev2$LowerBound>=0.8]=0.79999
tempTablev2_combine=rbind(tempTable,tempTablev2)

TableForReferenceFigure=rbind(TableForReferenceFigure,tempTablev2_combine)
TableForReferenceFigure=TableForReferenceFigure[order(TableForReferenceFigure$Reference,TableForReferenceFigure$Neut),]

TableForReferenceFigure$TechnicalName_Ref<-factor(TableForReferenceFigure$TechnicalName_Ref,levels=c("BNT162b2","mRNA-1273","ChAdOx1 nCoV-19"))
PointsTable$TechnicalName_Ref<-factor(PointsTable$TechnicalName_Ref,levels=c("BNT162b2","mRNA-1273","ChAdOx1 nCoV-19"))
NonInferiorityMargin$TechnicalName_Ref<-factor(NonInferiorityMargin$TechnicalName_Ref,levels=c("BNT162b2","mRNA-1273","ChAdOx1 nCoV-19"))


colorlist=c("chocolate","mediumorchid4","orangered1")

FigureS2<-ggplot(data=TableForReferenceFigure,aes(x=(10^(Neut))))+
  geom_errorbar(data=SummaryTable_Efficacy_NeutRatio_SD_SEM,aes(x=10^NeutRatio_Reported,ymin=Lower,ymax=Upper),alpha=0.1,width=0.05) +
  geom_errorbar(data=PointsTable,aes(x=10^NeutRatio_Reported,ymin=Lower,ymax=Upper),width=0.05) +
  geom_errorbarh(data=SummaryTable_Efficacy_NeutRatio_SD_SEM,aes(y=100*Efficacy,xmin=RatioReported_LB,xmax=RatioReported_UB),height=2,inherit.aes=FALSE,alpha=0.1) +
  geom_errorbarh(data=PointsTable,aes(y=100*Efficacy,xmin=RatioReported_LB,xmax=RatioReported_UB,color=TechnicalName_Ref),height=2,inherit.aes=FALSE) +
  geom_line(aes(y=100*(Efficacy)),linetype=1,color="red") +
  geom_point(data=SummaryTable_Efficacy_NeutRatio_SD_SEM,aes(x=10^NeutRatio_Reported,y=100*Efficacy),alpha=0.1) +
  geom_point(data=PointsTable,aes(x=10^NeutRatio_Reported,y=100*Efficacy)) +
  # geom_segment(data=NonInferiorityMargin,aes(x=0.125,xend=10^MarginForAbove80,y=80,yend=80),linetype=2) +
  geom_segment(data=NonInferiorityMargin,aes(x=10^MarginForAbove80,xend=10^MarginForAbove80,y=0,yend=100,color=TechnicalName_Ref),linetype=2) +
  geom_ribbon(aes(ymin=100*(LowerBound),ymax=100,fill=LowerBound>=0.80),color=NA,alpha=0.1) +
  scale_x_log10(breaks=xticks,labels=xticks) +
  scale_y_continuous(breaks=yticks) +
  coord_cartesian(ylim=c(20,105),xlim=c(0.1,14),expand=FALSE) +
  labs(y="Efficacy (%)",
       x="Neutralisation\n(fold of convalescence)") +
  #      alpha="Number individuals") +
  scale_alpha_manual(values=c(0.2,1)) +
  scale_fill_manual(values=c("red","dodgerblue"),guide="none") +
  scale_color_manual(values=colorlist,guide="none") +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray98"),
        panel.grid.minor = element_line(colour = "gray98"),
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
  ) +
  facet_rep_wrap(~TechnicalName_Ref,ncol=3)


ggsave(plot=FigureS2,filename="FigureS2.pdf",height=3,width=7)

