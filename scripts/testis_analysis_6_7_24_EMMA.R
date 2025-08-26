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

#boxplot for Treatment
boxplot(Data_Merge$`Length mm`~Data_Merge$Treatment,xlab="Caloric Density", ylab="Testes length (mm)", col= c("purple","blue","green"))
legend("topright", c("0.5x", "1x", "2x"),
       fill = c("purple", "blue", "green"))+
  mtext(c("A","B","B"),3,0.25,at=c(1,2,3),col="blue")+
  text(1,1.75,paste("Treatment: p=",1.01e-06,sep=""),col="blue")
  
#add posthoc letters

#add p-value from model


#boxplot for Age
boxplot(Data_Merge$`Length mm`~Data_Merge$Age.x,xlab="Age (Days)", ylab="Testes length (mm)", col= c("blue","orange","red"))
#add p-value from model
text(1,1.75,paste("Age: p=",0.042696,sep=""),col="blue")
#mtext(c("A","B","B"),3,0.25,at=c(1,2,3),col="blue")
legend("topright", c("0D", "2D", "5D"),
       fill = c("blue", "orange","red"))

#boxplot for Strain
boxplot(Data_Merge$`Length mm`~Data_Merge$Strain,xlab="Strain", ylab="Testes length (mm)", col= c("blue", "orange"))
#add posthoc letters
#mtext(c("A","B"),3,0.25,at=c(1,2),col="blue")
#add p-value from model
text(1,1.75,paste("Strain: p<",2e-16,sep=""),col="blue")
legend("topright", c("42", "217"),
       fill = c("blue", "orange"))

#boxplot for Age by Treatment
boxplot(Data_Merge$`Length mm`~Data_Merge$Treatment+Data_Merge$Age.x+Data_Merge$Strain,xlab="Caloric Density by Age (Days)", ylab="Testes length (mm)",sep="*", col = c("purple", "purple","purple", "blue","blue","blue","green","green","green"))

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

#add posthoc letters
#mtext(c("A","A","A","B","BC","C"),3,0.25,at=c(1,2,3,4,5,6),col="red")
#add p-value from model
text(1.9,1.75,paste("Treatment: p=",1.01e-06,sep=""),col="black")
text(1.9,1.65,paste("Age: p=",0.042696,sep=""),col="black")
text(1.9,1.55,paste("Treatment*Age: p=",0.000615,sep=""),col="black")


#boxplot for Strain by Treatment
boxplot(Data_Merge$`Length mm`~Data_Merge$Treatment+Data_Merge$Strain,xlab="Caloric Density by Strain", ylab="Testes length (mm)",sep="*", col=c("blue","blue","orange","orange"))
#add posthoc letters
#mtext(c("A","A","A","B","BC","C"),3,0.25,at=c(1,2,3,4,5,6),col="green")
#add p-value from model for Treatment
text(1.9,1.75,paste("Treatment: p=",1.01e-06,sep=""),col="black")
text(1.9,1.65,paste("Strain: p<",2e-16,sep=""),col="black")
text(1.9,1.55,paste("Treatment*Strain p=",0.013515,sep=""),col="black")

#boxplot for Age by Strain
boxplot(Data_Merge$`Length mm`~Data_Merge$Strain+Data_Merge$Age.x,xlab="Strain by Age (Days)", ylab="Testis length (mm)",sep="*", col= c("purple","purple","blue","blue","green","green"))
#add posthoc letters
#mtext(c("A","A","A","B","BC","C"),3,0.25,at=c(1,2,3,4,5,6),col="black")
#add p-value from model
text(1.9,1.75,paste("Age: p=",0.042696,sep=""),col="black")
text(1.9,1.65,paste("Strain: p<",2e-16,sep=""),col="black")
text(1.9,1.55,paste("Age*Strain: p=",1.50e-05,sep=""),col="black")

#boxplot for Treatment by Age by Strain
boxplot(Data_Merge$`Length mm`~Data_Merge$Treatment+Data_Merge$Age.x+Data_Merge$Strain,xlab="Caloric Density by Age (Days) by Strain", ylab="Testis length (mm)",sep="*")
#add posthoc letters
#mtext(c("A","A","A","B","BC","C"),3,0.25,at=c(1,2,3,4,5,6),col="orange")
#add p-value from model
text(1.9,1.75,paste("Treatment: p=",1.01e-06,sep=""),col="orange")
text(1.9,1.65,paste("Age: p=",0.042696,sep=""),col="orange")
text(1.9,1.55,paste("Strain: p<",2e-16,sep=""),col="orange")
text(2.1,1.45,paste("Treatment*Age*Strain: p=",0.012196,sep=""),col="orange")

