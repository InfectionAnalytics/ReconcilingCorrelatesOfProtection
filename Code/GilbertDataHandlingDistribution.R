### Importing and Handling Gilbert data

##Estimating distribution of data and SD.


##Importing
GilbertcID50Raw<-read.csv(paste0(currentdirectory,"/Data/GilbertcID50_extracted.csv"),fileEncoding="UTF-8-BOM")


###Gilbert Data Handling

###Gilbert Study Means
ModernaMean_IC50=247 #From paper
ModernaMean_IC80=478 #From paper

##### Estimating Distribution of Gilbert Data
ModernaMean_IC50_upper=264 #From paper
ModernaMean_IC50_lower=231 #From paper
ModernaMean_IC80_upper=508 #From paper
ModernaMean_IC80_lower=450 #From paper
GilbertNumberIndividiual=1005 #From paper


###Gilbert - Estimating Means and SD with with and without censoring.
NeutModernaControl_sd_v1=((log10(ModernaMean_IC50_upper)-log10(ModernaMean_IC50))/1.96)*sqrt(GilbertNumberIndividiual)
NeutModernaControl_sd_v2=((log10(ModernaMean_IC50)-log10(ModernaMean_IC50_lower))/1.96)*sqrt(GilbertNumberIndividiual)
NeutModernaControlID80_sd_v1=((log10(ModernaMean_IC80_upper)-log10(ModernaMean_IC80))/1.96)*sqrt(GilbertNumberIndividiual)
NeutModernaControlID80_sd_v2=((log10(ModernaMean_IC80)-log10(ModernaMean_IC80_lower))/1.96)*sqrt(GilbertNumberIndividiual)

NeutModernaControlID80_sd=mean(NeutModernaControlID80_sd_v1,NeutModernaControlID80_sd_v2)
NeutModernaControl_sd=mean(NeutModernaControl_sd_v1,NeutModernaControl_sd_v2)


GilbertCases=GilbertcID50Raw$cases
GilbertCases=GilbertCases[!is.na(GilbertCases)]
GilbertControl=GilbertcID50Raw$Control

###LOD and LOQ from paper
LODGilbertcID50=2.42 #From paper
LOQGilbertcID50=4.477 #From paper
NormGilbert_Control_ID50<-nlm(function(p){-sum(log(dnorm(log10(GilbertControl[GilbertControl>LODGilbertcID50]),p[1],exp(p[2]))))-sum(GilbertControl<=LODGilbertcID50)*log(pnorm(log10(LODGilbertcID50),p[1],exp(p[2])))},c(mean(log10(GilbertControl)),log(sd(log10(GilbertControl)))))
NormGilbert_Case_ID50<-nlm(function(p){-sum(log(dnorm(log10(GilbertCases[GilbertCases>LODGilbertcID50]),p[1],exp(p[2]))))-sum(GilbertCases<=LODGilbertcID50)*log(pnorm(log10(LODGilbertcID50),p[1],exp(p[2])))},c(mean(log10(GilbertCases)),log(sd(log10(GilbertCases)))))




