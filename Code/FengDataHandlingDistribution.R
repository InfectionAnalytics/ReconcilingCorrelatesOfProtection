### Importing and Handling Feng data

##Estimating distribution of data and SD.

FengData<-read.csv(paste0(currentdirectory,"/Data/ExtractedFengData.csv"),fileEncoding="UTF-8-BOM")


####Feng Data Handling

### Adding back excluded zeros from figure in Feng et al.
Assay=c(rep("ID50",5),rep("NF50",5))
Group=c("NAAT Negative","NAAT Positive","Primary","Asymptomatic","Non-Primary","NAAT Negative","NAAT Positive","Primary","Asymptomatic","Non-Primary")
NumberIndividuals=c(828,149,47,86,16,412,110,36,62,12)
FengLOD_Pseudo=40
FengLOD_Live=8.6

for (i in 1:length(Group)) {
  ZeroList<-FengData[1:(NumberIndividuals[i]-length(FengData$NeutLevel[FengData$Assay==Assay[i] & FengData$Group==Group[i]])),]  
  ZeroList$x<-NA
  ZeroList$y<-NA
  ZeroList$Panel<-NA
  ZeroList$Assay=Assay[i]
  ZeroList$Group=Group[i]
  if (Assay[i]=="ID50") {
    ZeroList$NeutLevel<-FengLOD_Pseudo/2  
  } else {
    ZeroList$NeutLevel<-FengLOD_Live/2
  }
  
  FengData<-rbind(FengData,ZeroList)
}

####


##### Feng et al.  Means
#### (The Pseudo virus mean reported was not in IU but figures for this marker were all in IU, 
#### the factor for transformation was given in methods)
FengMean_Pseudo=158 #From paper/supplement
FengMean_Live=184 #From paper/supplement
FengConversionFactorPseudoIU=0.1428 #From paper/supplement
FengMean_Pseudo_IU=FengMean_Pseudo*FengConversionFactorPseudoIU


### Estimating Distribution of Feng Data
SDFeng_Live_temp<-nlm(function(p){-sum(log(dnorm(log10(FengData$NeutLevel[FengData$NeutLevel>FengLOD_Live & FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"]),p[1],exp(p[2]))))-length(FengData$NeutLevel[FengData$NeutLevel<=FengLOD_Live & FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"])*log(pnorm(log10(FengLOD_Live),p[1],exp(p[2])))},c(mean(log10(FengData$NeutLevel[FengData$NeutLevel>0 & FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"])),log(sd(log10(FengData$NeutLevel[FengData$NeutLevel>0 & FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"])))))
SDFeng_Live=exp(SDFeng_Live_temp$estimate[2])
MeanFeng_Live=SDFeng_Live_temp$estimate[1]
SDFeng_Pseudo_temp<-nlm(function(p){-sum(log(dnorm(log10(FengData$NeutLevel[FengData$NeutLevel>FengLOD_Pseudo & FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"]),p[1],exp(p[2]))))-length(FengData$NeutLevel[FengData$NeutLevel<=FengLOD_Pseudo & FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"])*log(pnorm(log10(FengLOD_Pseudo),p[1],exp(p[2])))},c(mean(log10(FengData$NeutLevel[FengData$NeutLevel>0 & FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"])),log(sd(log10(FengData$NeutLevel[FengData$NeutLevel>0 & FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"])))))
SDFeng_Pseudo=exp(SDFeng_Pseudo_temp$estimate[2])
MeanFeng_Pseudo=SDFeng_Pseudo_temp$estimate[1]

