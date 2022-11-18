This repositry contains:
Data Folder:
All the data files containing extracted and preivousy published data from the vaccine comparison model (Khoury et al. Nature Medicine 2021), and the breakthrough infeciton studies (Bergwerk et al. NEJM, Feng et al. Nature Medicine and Gilbert et al. Science). 

  1. Bergwerk_FalseData.csv - A file containing the format of the data provided by Bergwerk et al. from their study published in NEJM - and as specified in the original publication the data can be requested from the authors of the study. This file has had the real data replaced by false placeholder data in order to illustrate the analysis pipeline.
  
  2. ExtractedFengData.csv - file containing the extracted Feng et al. Nature Medicine Data (extracted with high accuracy using Adobe Illustrator).
  
  3. FengGilbertCurves_Long.csv - file containing the manually extracted Feng et al. and Gilbert et al. models of VE with neutralisation (Extracted using Webplotdigitizer)
  
  4. GilbertcID50_extracted.csv - file containing the manually extracted data from Gilbert et al.  (Extracted using Webplotdigitizer)
  
  5. SummaryTable_Efficacy_NeutRatio_SD_SEM.csv - file containing the neutralisation titres from vaccinated and convalescent groups of Phase I/II trials and matching VE from Phase III trials of the same vaccines (as reported in Khoury et al.)
  

Code: Folder containing all R scripts used to run analysis and plot figures.
  
  1. RunKhouryModel.R - Script runs the model from Khoury et al. and generates the model curves and CIs for plotting later.
  
  2. FengDataHandlingDistribution.R - Imports the Feng et al. data and generates key statistics (means, SD of distribution).
  
  3. GilbertDataHandlingDistribution.R - Imports the Gilbert et al. data and generates key statistics (means, SD of distribution).
  
  4. Figure2.R - Script imports the extracted Gilbert and Feng model data from "FengGilbertCurves_Long.csv", and plots this data with confidence intervals compared with the Khoury et al. Model (above).
  
  5. FigureS3.R - generates figure S3
  
  6. FigureS4.R - fits the extracted Feng et al. and Gilbert et al. models (from "FengGilbertCurves_Long.csv") with the function describing the form of these curves and generates figure S4.
  
  7. FigureS1.R - Uses the Gilbert and Feng models (fitted in Figure S4) and the distribution in neutralisation data (estimated in 2 and 3 above) to estimate the overall vaccine efficacy observed in a population, and compare this to reported vaccine efficacy for different vaccines and the Khoury et al. model (generates Figure S1).
  
  8. Figure3.R - Uses much of the above analysis to estimate from, the breakthrough infection studies, the proportion of people with breakthrough infecitons in each 2-fold neutralisation bucket (along with overall vaccine efficacy) to give a "raw" estimate of the efficacy at each neutralisation titer for comparison with the data.
  
  9. FigureS5.R - generates the QQ-plots showing the normality of neutralisation data from all avaliable sources in this analysis.
  
  10. FigureS2.R - uses the inverse Khoury et al. model to estimate the non-inferiority margins to deliver a specified vaccine efficacy. (generates Table S2 and Figure S2)


To run code:
Generating the analysis can be done by running the scripte: "ReconcilingCorrelates_Master_20221118.R"
  - This will run all analysis and generate all figures.
  


