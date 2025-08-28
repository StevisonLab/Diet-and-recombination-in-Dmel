#script to analyze testis images

#load packages
library(stringr)
library(lme4)
library(car)
library(emmeans)
library(multcompView)
library(multcomp)
library(ggplot2)
library(ggpubr)
library(agricolae)
library(doBy)
library(readr)
library(ggthemes)

# Read in Data ------------------------------------------------------------


#set working directory
#script location

#data location
##Needs 



#put both raw data and meta data into rawdata folder

#read in the meta data
meta <- read_csv("rawdata/Testis_meta_data.csv")
#read in the raw data
raw<-read_csv("rawdata/Testis_raw_9_13.csv")


#before merging
#need to make a column in the meta data that matches the column in raw (seems to be a concatenation of vial number, age and filename without the file extension)
#function to looks at are "paste" for concatenating and the package "stringr" for find and replace to get rid of file extension
meta$Picture_Code=paste(meta$Vial_Code,meta$Age,meta$Picture_Number,sep="_")
#con<-meta$Picture_Code
meta$Picture_Code=str_remove(meta$Picture_Code,".tif")
#merged datasets using function "merge"
Data_Merge = (merge(meta,raw,by="Picture_Code",all=F))

#make model parameters factors
Data_Merge$Treatment=as.factor(Data_Merge$Treatment)
Data_Merge$Age.x=as.factor(Data_Merge$Age.x)
Data_Merge$Strain=as.factor(Data_Merge$Strain)
Data_Merge$Length=Data_Merge$`Length mm`

#Summary of the data for sample sizes
summaryBy(Length~ Strain+Treatment+Age.x, data = Data_Merge)
#Data_Merge$Length=Data_Merge$`Length mm`
d<-summaryBy(Length~ Strain + Treatment, data = Data_Merge, FUN = length)

# Statistical Analysis ----------------------------------------------------


#run anova on data
#model: linear model (lm) of length ~ Treatment * Day * Strain
#function: 'aov'
#example: aov(yield ~ block + N * P + K, data = npk)
#complete this module of mini course: https://stevisonlab.github.io/R-Mini-Course/docs/module4.html
#model_length=aov(Data_Merge$Length.mm~Data_Merge$Treatment+Data_Merge$Age+Data_Merge$Strain)
#or
model_length=aov(`Length mm`~Treatment*Age.x*Strain, data=Data_Merge)
anova_table=summary(model_length)
anova_table
cld2<-emmeans(model_length, ~Treatment+Strain)
cld3 <- cld(cld2, Letters = letters, alpha = 0.05)
d$CLD=cld3$.group


##Boxplot with CLD by treatment and strain, Supplemental Figure 5
col_pal2=c("#fdd0a2","#fd8d3c","#a63603","#c6dbef","#6baed6","#08519c")
Data_Merge$treat=ifelse(Data_Merge$Treatment=="0.5x",0.5,ifelse(Data_Merge$Treatment=="1x",1,2))
ggplot(data=Data_Merge, aes(x=Treatment, y=Length, fill=Treatment))+
  ylab("Testes Length (mm)") +
  xlab("Caloric Density") + theme_base() + facet_wrap(~Strain) +
  geom_boxplot() + scale_fill_manual(values = col_pal2[1:3]) +
  geom_text(data = d, aes(y=1.5, x=Treatment, label=CLD), position=position_dodge2(width = 0.75), color = "black", fontface = "bold")
ggsave("images/Testes_42.png", height=5)  


ggplot(data=Data_Merge, aes(x=Treatment, y=Length, fill=Treatment))+
  ylab("Testes Length (mm)") +
  xlab("Caloric Density") + theme_base() + facet_wrap(~Strain) +
  geom_boxplot() + scale_fill_manual(values = col_pal2[4:6]) +
  geom_text(data = d, aes(y=1.5, x=Treatment, label=CLD), position=position_dodge2(width = 0.75), color = "black", fontface = "bold")

ggsave("images/Testes_217.png", height=5)  

