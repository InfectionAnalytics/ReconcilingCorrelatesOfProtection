

###Figure S3

##### Gilbert vs Feng Curves
FigureS3<-ggplot(data=EstimatedEfficacy,aes(x=(10^(NeutValuelog10New)))) +
  geom_polygon(data=tempGilbertCIs[tempGilbertCIs$Assay=="cID50",],aes(x=Neut,y=(Efficacy),fill=Assay),color=NA,alpha=0.15) +
  geom_polygon(data=tempFengCIs[tempFengCIs$Assay=="Pseudo",],aes(x=Neut,y=(Efficacy),fill=Assay),color=NA,alpha=0.15) +
  geom_line(data=FengGilbertCurves[FengGilbertCurves$Study  %in% c("Gilbert","Feng") & FengGilbertCurves$Assay %in% c("Pseudo","cID50") & FengGilbertCurves$Item=="mean",],aes(x=Neut,y=Efficacy,color=Assay)) +
  scale_x_log10() +
  scale_y_continuous(breaks=yticks) +
  scale_color_manual(values=colorlist) +
  scale_fill_manual(values=colorlist) +
  coord_cartesian(ylim=c(40,105),expand=FALSE) +
  labs(y="Vaccine Efficacy (%)",
       x="Neutralisation\n(IU/mL)") +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(angle=45,hjust=1),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
  ) 


ggsave(plot=FigureS3,filename="FigureS3.pdf",height=3.4,width=4)
