#Plot efficacy figures (move per variant)
library(stringr)
library(mvtnorm)
library(lemon)
library(ggplot2)
library(magic)
library(plyr)
library(Hmisc)
library(splines)
library(foreach)
library(doParallel)


#### Set directory to same directory as the r-script
currentdirectory<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentdirectory)

#Import datasets and correlate curves
source(paste0(currentdirectory,"/Code/RunKhouryModel.R")) #Run Khoury et al. Model and estimate CIs
source(paste0(currentdirectory,"/Code/FengDataHandlingDistribution.R")) #Get Means and SD of Neut data distribution
source(paste0(currentdirectory,"/Code/GilbertDataHandlingDistribution.R")) #Get Means and SD of Neut data distribution
source(paste0(currentdirectory,"/Code/Figure2.R")) ### Visualise Khoury,Feng Gilbert curves normalised
source(paste0(currentdirectory,"/Code/FigureS3.R")) ##Gilbert vs Feng plot
source(paste0(currentdirectory,"/Code/FigureS4.R")) ##Approximate the Feng and Gilbert Efficacy Model functions in order to use to estimate population efficacy (in Figure S1)
source(paste0(currentdirectory,"/Code/FigureS1.R")) ##Estimate a vaccinated cohorts overall protection from each model
source(paste0(currentdirectory,"/Code/Figure3.R"))  ##Generate point estimates from each Breakthrough study of the %efficacy and compare to models
source(paste0(currentdirectory,"/Code/FigureS5.R")) ##Normality plots
source(paste0(currentdirectory,"/Code/FigureS2.R")) ##Estimate Non-inferiority margins




##Uncertainity in neut titres for individual
sdofmeasure=0.41235756
PredictBand=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2]+qnorm(c(0.025,0.975),0,sdofmeasure)
DistributionOfAverageOfDuplicates=qnorm(c(0.025,0.5,0.975),0,0.5*sqrt(2*(sdofmeasure^2)))
PredictBand_duplicate=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2]+DistributionOfAverageOfDuplicates
EfficacyBand_duplicate=ProbRemainUninfected(PredictBand_duplicate,
                                            tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                            tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])

