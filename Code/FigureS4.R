
####Approximating the VE models from the Feng and Gilbert 
#Papers (used for estimating population immunity in figure S1)

###Also make - Figure S4


#### Functional form of Gilbert and Feng Models
BreakthroughEfficacyFunction=function(s,a,b,c){
  CVE=1-(c*(1-exp(-a*exp(-b*s))))
  
  CVE
}
###Fitting Breakthrough efficacy function to extracted model
SSE=function(log10neut,Efficacy,p){
  
  sum((log(Efficacy)-log(BreakthroughEfficacyFunction(log10neut,p[1],p[2],p[3])))^2)
}





###Fitting extracted VE curves from Feng and Gilbert 
# with functional Form
ApproximateEfficacyFunctionFit<-nlm(function(p){SSE(log10(FengGilbertCurves$Neut[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID50" & FengGilbertCurves$Item=="mean"]),FengGilbertCurves$Efficacy[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID50" & FengGilbertCurves$Item=="mean"]/100,exp(c(p)))},c(log(0.001),log(5),log(1/0.01003131)))
ApproximateEfficacyFunctionFit_cID80<-nlm(function(p){SSE(log10(FengGilbertCurves$Neut[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID80" & FengGilbertCurves$Item=="mean"]),FengGilbertCurves$Efficacy[FengGilbertCurves$Study=="Gilbert" & FengGilbertCurves$Assay=="cID80" & FengGilbertCurves$Item=="mean"]/100,exp(c(p)))},c(log(0.001),log(5),log(1/0.01003131)))
ApproximateFeng_Pseudo<-nlm(function(p){SSE(log10(FengGilbertCurves$Neut[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Pseudo" & FengGilbertCurves$Item=="mean"]),FengGilbertCurves$Efficacy[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Pseudo" & FengGilbertCurves$Item=="mean"]/100,exp(c(p)))},c(log(0.001),log(1),log(1/0.01003131)))
ApproximateFeng_Live<-nlm(function(p){SSE(log10(FengGilbertCurves$Neut[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Live" & FengGilbertCurves$Item=="mean"]),FengGilbertCurves$Efficacy[FengGilbertCurves$Study=="Feng" & FengGilbertCurves$Assay=="Live" & FengGilbertCurves$Item=="mean"]/100,exp(c(p)))},c(log(0.001),log(1),log(1/0.01003131)))

###Evaluating fitted models for visualisation
neutvector=seq(log10(1),log10(100000),by=0.01)
IndividualRiskFittedModel<-data.frame("Study"=rep(c("Gilbert","Feng"),each=2*length(neutvector)),
                                      "Assay"=rep(c("cID50","cID80","Pseudo","Live"),each=length(neutvector)),
                                      "Neut"=rep(neutvector,4),
                                      "Efficacy"=NA)
IndividualRiskFittedModel$Efficacy[IndividualRiskFittedModel$Study=="Gilbert" & IndividualRiskFittedModel$Assay=="cID50"]=BreakthroughEfficacyFunction(IndividualRiskFittedModel$Neut[IndividualRiskFittedModel$Study=="Gilbert" & IndividualRiskFittedModel$Assay=="cID50"],
                                                                                                                                                       exp(ApproximateEfficacyFunctionFit$estimate[1]),
                                                                                                                                                       exp(ApproximateEfficacyFunctionFit$estimate[2]),
                                                                                                                                                       exp(ApproximateEfficacyFunctionFit$estimate[3]))
IndividualRiskFittedModel$Efficacy[IndividualRiskFittedModel$Study=="Gilbert" & IndividualRiskFittedModel$Assay=="cID80"]=BreakthroughEfficacyFunction(IndividualRiskFittedModel$Neut[IndividualRiskFittedModel$Study=="Gilbert" & IndividualRiskFittedModel$Assay=="cID80"],
                                                                                                                                                       exp(ApproximateEfficacyFunctionFit_cID80$estimate[1]),
                                                                                                                                                       exp(ApproximateEfficacyFunctionFit_cID80$estimate[2]),
                                                                                                                                                       exp(ApproximateEfficacyFunctionFit_cID80$estimate[3]))
IndividualRiskFittedModel$Efficacy[IndividualRiskFittedModel$Study=="Feng" & IndividualRiskFittedModel$Assay=="Pseudo"]=BreakthroughEfficacyFunction(IndividualRiskFittedModel$Neut[IndividualRiskFittedModel$Study=="Feng" & IndividualRiskFittedModel$Assay=="Pseudo"],
                                                                                                                                                     exp(ApproximateFeng_Pseudo$estimate[1]),
                                                                                                                                                     exp(ApproximateFeng_Pseudo$estimate[2]),
                                                                                                                                                     exp(ApproximateFeng_Pseudo$estimate[3]))
IndividualRiskFittedModel$Efficacy[IndividualRiskFittedModel$Study=="Feng" & IndividualRiskFittedModel$Assay=="Live"]=BreakthroughEfficacyFunction(IndividualRiskFittedModel$Neut[IndividualRiskFittedModel$Study=="Feng" & IndividualRiskFittedModel$Assay=="Live"],
                                                                                                                                                   exp(ApproximateFeng_Live$estimate[1]),
                                                                                                                                                   exp(ApproximateFeng_Live$estimate[2]),
                                                                                                                                                   exp(ApproximateFeng_Live$estimate[3]))




####Plot of functional form of model
FigureS4<-ggplot(data=FengGilbertCurves[FengGilbertCurves$Item=="mean",],aes(x=Neut,y=Efficacy,color=Assay)) +
  geom_point() +
  geom_line(data=IndividualRiskFittedModel,aes(x=10^Neut,y=100*Efficacy,color=Assay)) +
  scale_x_log10()+
  coord_cartesian(ylim=c(40,100),xlim=c(1,10^4)) +
  facet_rep_wrap(~Study,ncol=2) +
  labs(x="Neutralisation titre (as reported in study)",
       y="Efficacy (%)",
       color="Assay") +
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
        # legend.spacing.y = unit(0, 'cm'),
        # legend.key = element_rect(size = 1),
        # legend.key.size = unit(0.45, "cm"),
        # legend.margin = margin(t=0,b=0.7,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
  )


ggsave(plot=FigureS4,filename = "FigureS4.pdf",height=4,width=7)
