### Figure S5

# Generate the qqplot showing normality of extracted data



###### Normality plots
pdf("FigureS5.pdf",height=4,width=10)
par(mfrow=c(1,3))
qqplotdata<-qqnorm(log10(GilbertControl),main="Gilbert et al.")
qqline(log10(GilbertControl))
abline(v=head(sort(qqplotdata$x),sum(GilbertControl<=LOQGilbertcID50*1.1))[sum(GilbertControl<=LOQGilbertcID50*1.1)],col="red")

qqplotdata<-qqnorm(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"]),main="Feng et al. (Pseudo)")
qqline(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"]))
abline(v=head(sort(qqplotdata$x),sum(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"]<=FengLOD_Pseudo))[sum(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="ID50"]<=FengLOD_Pseudo)],col="red")

qqplotdata<-qqnorm(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"]),main="Feng et al. (Live)")
qqline(log10(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"]))
abline(v=head(sort(qqplotdata$x),sum(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"]<=FengLOD_Live))[sum(FengData$NeutLevel[FengData$Group=="NAAT Negative" & FengData$Assay=="NF50"]<=FengLOD_Live)],col="red")
dev.off()

