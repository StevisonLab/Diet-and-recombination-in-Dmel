##RQ Nutrition Project Spring 2023
library(readr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())
library(doBy)
library(dplyr)
library(tidyr)
library(scales)
library(agricolae)
library(broom.mixed)
library(lme4)
library(car)
library(emmeans)
library(multcompView)
library(multcomp)

resp <- read_csv("rawdata/RQ_Spring2023.csv", na = "NA")

resp$Weight = resp$`Pre-Weight` * 1000
resp$Stock= as.factor(resp$Stock)
resp$Treatment= as.factor(resp$Treatment)

##Convert raw data to ml
resp$Area_O2cor = (resp$TVO2)/100
resp$Area_CO2cor = (resp$TVCO2)/100

##Calculate O2 and CO2 per hour per grams and RQ
resp$Vol_Correction = resp$`Chm Vol`+ 0.048 + 0.03 - (10.1*resp$`Pre-Weight`)
resp$O2pHr = (resp$Area_O2cor*((resp$Vol_Correction/resp$`Inj Vol`)/resp$Time))
resp$CO2pHr = (resp$Area_CO2cor*((resp$Vol_Correction/resp$`Inj Vol`)/resp$Time))
resp$O2pHrpGr = (resp$O2pHr/resp$`Pre-Weight`)
resp$CO2pHrpGr = (resp$CO2pHr/resp$`Pre-Weight`)
resp$RQ = resp$CO2pHrpGr/resp$O2pHrpGr

#Check variables trends and distributions before and after conversions
hist(as.numeric(resp$TVO2))
hist(as.numeric(resp$TVCO2))
hist(as.numeric(resp$Area_O2cor))
hist(as.numeric(resp$Area_CO2cor))
hist((resp$O2pHr))
hist((resp$CO2pHr))
hist(resp$O2pHrpGr)
hist(resp$CO2pHrpGr)
hist(as.numeric(resp$RQ))

plot(resp$TVO2~resp$`Pre-Weight`)
plot(resp$Area_O2cor~resp$`Pre-Weight`)
plot(resp$O2pHr~resp$`Pre-Weight`)
plot(resp$O2pHrpGr~resp$`Pre-Weight`)
plot(resp$TVCO2~resp$`Pre-Weight`)
plot(resp$Area_CO2cor~resp$`Pre-Weight`)
plot(resp$CO2pHr~resp$`Pre-Weight`)
plot(resp$CO2pHrpGr~resp$`Pre-Weight`)
hist(resp$RQ)
##Plot CO2pHr and O2pHr in ul against initial weight
plot(resp$`Pre-Weight`, resp$CO2pHr*1000)
plot(resp$`Pre-Weight`, resp$O2pHr*1000)

#Check correlation
cor(resp$Area_O2cor, resp$`Pre-Weight`,  method = "pearson", use = "complete.obs")

## Convert grams to mg and ml to ul
resp$Weight = resp$`Pre-Weight` * 1000
resp$CO2pHr2 = resp$CO2pHr *1000
resp$O2pHr2 = resp$O2pHr *1000
#pdf("RQ_Summary_Summer_2022.pdf")

##sub-setting data with RQ < 1.5
RQW<-subset(resp, resp$RQ>2, na.rm=TRUE, select = c(Date, Round, `Syringe No`, Species, Sex, Vial, Treatment, Weight, O2pHr2, CO2pHr2, RQ, Weight, Area_O2cor, Area_CO2cor))
write_csv(RQW, "output/CleanUp.csv")
resp<-subset(resp, resp$RQ<2, na.rm=TRUE, select = c(Date, Round, Stock, Treatment, Species, Sex, Vial, Weight, O2pHr2, CO2pHr2, RQ, Area_O2cor, Area_CO2cor))

#Function to add N to the boxplot of RQ, O2, CO2, and Weight 
#RQ
stat_box_data <- function(y, upper_limit = 1.62) {
  return( 
    data.frame(
      y = 1.55,
      label = paste('n =', length(y))
    )
  )
}
#O2
stat_box_data1 <- function(y, upper_limit = 5) {
  return( 
    data.frame(
      y = 0.8 * upper_limit,
      label = paste('n =', length(y))
    )
  )
}
#CO2
stat_box_data2 <- function(y, upper_limit = 300) {
  return( 
    data.frame(
      y = 0.8 * upper_limit,
      label = paste('n =', length(y))
    )
  )
}
#Weight
stat_box_data3 <- function(y, upper_limit = max(resp$Weight) * 0.1) {
  return( 
    data.frame(
      y =1.6,
      label = paste('n =', length(y))
    )
  )
}
pdf("images/RQ_Updated_10_4.pdf")

RQst=ggplot(resp, aes(x=Sex, y=RQ, fill=Stock)) + 
  ylab(bquote(Respiratory~Quotient~(RQ))) +
  scale_fill_manual(breaks = c("42", "217"), values=c("blue","orange")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_text(data = RQS.summarized, aes(y = 1.5, x = Sex, label = group), position = position_dodge(width = 0.75)) +
  geom_boxplot()
RQst

RQBoth1=ggplot(resp, aes(x=Sex, y=CO2pHr2, fill=Stock)) + 
  ylab(bquote(ul ~ CO[2] ~ h^-1 ~ produced)) +
  scale_fill_manual(breaks = c("42", "217"), values=c("blue","orange")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.data = stat_box_data2, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_text(data = CO2S.summarized, aes(y = 220, x = Sex, label = group), position = position_dodge(width = 0.75)) +
  geom_boxplot() 
RQBoth1

RQBoth2=ggplot(resp, aes(x=Sex, y=O2pHr2, fill=Stock)) + 
  ylab(bquote(ul ~ O[2] ~ h^-1 ~ consumed)) + scale_fill_manual(breaks = c("42", "217"), values=c("blue","orange")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_text(data = O2S.summarized, aes(y = 290, x = Sex, label = group), position = position_dodge(width = 0.75)) +
  geom_boxplot() 
RQBoth2

RQBoth2=ggplot(resp, aes(x=Sex, y=Weight, fill=Stock)) + 
  ylim(0.2,1.8)+
  ylab(bquote(Body ~ mass ~ (mg))) + scale_fill_manual(breaks = c("42", "217"), values=c("blue","orange")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.data = stat_box_data3, geom = "text", fun = max,
               position = position_dodge(width = 0.75)) +
  geom_boxplot() 
RQBoth2

##By treatment:sex
RQst=ggplot(resp, aes(x=Treatment, y=RQ, fill=Sex)) + 
  ylab(bquote(bold(.(Respiratory~Quotient~(RQ))))) +
  scale_fill_manual(breaks = c("F", "M"), values=c("salmon","green4")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic() +
  theme(
  axis.text.y = element_text(size=8, face= "bold", colour= "black", family = "sans" ),
  axis.text.x = element_text(size=10, face="bold", colour = "black", family = "sans" )) +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_text(data = RQS.summarized, aes(y = 1.5, x = Treatment, label = group, fontface = "bold",  family = "sans"), position = position_dodge(width = 0.75)) +
  geom_boxplot(size=0.8)+theme( strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"), strip.text.x = element_text(size = 10, color = "black", face = "bold"))
RQst

RQBoth1=ggplot(resp, aes(x=Treatment, y=CO2pHr2, fill=Sex)) + 
  ylab(bquote(bold(.(ul ~ CO[2] ~ h^-1)))) +
  scale_fill_manual(breaks = c("F", "M"), values=c("salmon","green4")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size=8, face= "bold", colour= "black", family = "sans" ),
    axis.text.x = element_text(size=10, face="bold", colour = "black", family = "sans" )) +
  geom_text(data = CO2S.summarized, aes(y = 220, x = Treatment, label = group, fontface = "bold",  family = "sans"), position = position_dodge(width = 0.75)) +
  geom_boxplot(size=0.8)+theme( strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"), strip.text.x = element_text(size = 10, color = "black", face = "bold"))
RQBoth1

RQBoth2=ggplot(resp, aes(x=Treatment, y=O2pHr2, fill=Sex)) + 
  ylab(bquote(bold(.(ul ~ O[2] ~ h^-1)))) +
  xlab(bquote(bold(.("Treatment")))) +
  scale_fill_manual(breaks = c("F", "M"), values=c("salmon","green4")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size=8, face= "bold", colour= "black", family = "sans" ),
    axis.text.x = element_text(size=10, face="bold", colour = "black", family = "sans" )) +
  geom_text(data = O2S.summarized, aes(y = 290, x = Treatment, label = group, fontface = "bold",  family = "sans"), position = position_dodge(width = 0.75)) +
  geom_boxplot(size=0.8)+theme( strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"), strip.text.x = element_text(size = 10, color = "black", face = "bold"))
RQBoth2

#Subset by stock 
resp42<-subset(resp, resp$Stock=="42", na.rm=TRUE, select = c(Species, Sex, Vial, Treatment, O2pHr2, CO2pHr2, RQ, Weight, Area_O2cor, Area_CO2cor))     
resp217<-subset(resp, resp$Stock=="217", na.rm=TRUE, select = c(Species, Sex, Vial, Treatment, O2pHr2, CO2pHr2, RQ, Weight, Area_O2cor, Area_CO2cor))

##Summary of the data
fun <- function(x){
  c(m=mean(x), v=var(x), n=length(x), min=min(x), max=max(x))
}
summary_by(cbind(O2pHr2,CO2pHr2, Weight, RQ) ~ Treatment, data=resp42, FUN=fun)
summary_by(cbind(O2pHr2,CO2pHr2, Weight, RQ) ~ Treatment, data=resp217, FUN=fun)

SumAll <- summary_by(cbind(O2pHr2,CO2pHr2, RQ) ~ Treatment + Sex + Stock, data=resp, FUN=fun)
write.csv(SumAll, "output/SumAll_RQ.csv")

#RQ for 42
RQTreatment42=ggplot(resp42, aes(x=Sex, y=RQ, fill=Treatment)) + 
  ggtitle("RQ 42") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Respiratory Quotient (RQ)") +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_boxplot() 
RQTreatment42

#O2 for 42
O2_42=ggplot(resp42, aes(x=Sex, y=O2pHr2, fill=Treatment)) + 
  ylab(bquote(ul ~ O[2] ~ h^-1)) +
  ggtitle("Oxygen Consumed 42") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot()
O2_42

#CO2 for 42
CO2_42=ggplot(resp42, aes(x=Sex, y=CO2pHr2, fill=Treatment)) + 
  ggtitle("Carbon Dioxide produced 42") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(bquote(ul ~ CO[2] ~ h^-1)) +
  geom_boxplot()
CO2_42

#42 Weight at day 10
W_42=ggplot(resp42, aes(x=Sex, y= Weight, fill=Treatment)) +
  ggtitle("Body mass 42") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(bquote(Body ~ mass ~ (mg))) +
  geom_boxplot() 
W_42

#RQ for 217
RQTreatment217=ggplot(resp217, aes(x=Sex, y=RQ, fill=Treatment)) + 
  ylab("Respiratory Quotient (RQ)") +
  ggtitle("RQ 217") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_boxplot() 
RQTreatment217

#O2 for 217
O2_217=ggplot(resp217, aes(x=Sex, y=O2pHr2, fill=Treatment)) + 
  ylab(bquote(ul ~ O[2] ~ h^-1)) +
  ggtitle("Oxygen Consumed 217") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_boxplot() 
O2_217 

#CO2 for 217
CO2_217=ggplot(resp217, aes(x=Sex, y=CO2pHr2, fill=Treatment)) + 
  ggtitle("Carbon Dioxide produced 217") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(bquote(ul ~ CO[2] ~ h^-1)) +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_boxplot() 
CO2_217

#217 Weight
W_217=ggplot(resp217, aes(x=Sex, y=Weight, fill=Treatment)) + 
  ggtitle("Body mass 217") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(bquote(Body ~ mass ~ (mg))) +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_boxplot() 
W_217

dev.off()

####Linear mixed models pooling all data
res.aov <- lmer(RQ ~ (1|Vial) + Treatment * Sex * Stock, data = resp)
as.data.frame(summary(res.aov)$coefficients)
RQ1 <- Anova(res.aov)
RQ1
RQpd <- emmeans(res.aov, ~ Sex * Treatment)
summary(RQpd)
cldRQ <-cld(RQpd, Letters = letters, alpha = 0.05)
plot(cldRQ)
RQS=summaryBy(RQ~Sex+Treatment,data=resp,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
RQS.summarized=merge(RQS,cldRQ, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
unique(RQS.summarized$.group)
RQS.summarized$group=ifelse(RQS.summarized$.group==" ab  ","ab", ifelse(RQS.summarized$.group==" abc ","abc", ifelse(RQS.summarized$.group==" a   ","a", ifelse(RQS.summarized$.group=="   cd","cd", ifelse(RQS.summarized$.group=="  bcd","bcd", ifelse(RQS.summarized$.group=="    d","d",RQS.summarized$.group))))))

res.aov2 <- lmer(O2pHr2 ~ (1|Vial) + Treatment * Sex * Stock, resp)
summary(res.aov2)
as.data.frame(summary(res.aov2)$coefficients)
O2 <- Anova(res.aov2)
O2
O2_1 <- emmeans(res.aov2, ~Sex*Treatment)
plot(O2_1)
cldO2 <-cld(O2_1, Letters = letters, alpha = 0.05)
O2S=summaryBy(O2pHr2~Treatment+Sex,data=resp,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
O2S.summarized=merge(O2S,cldO2, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
unique(O2S.summarized$.group)
O2S.summarized$group=ifelse(O2S.summarized$.group=="  bc","bc", ifelse(O2S.summarized$.group==" a  ","a", ifelse(O2S.summarized$.group==" abc","abc",  ifelse(O2S.summarized$.group==" ab ","ab", ifelse(O2S.summarized$.group=="   c","c", O2S.summarized$.group)))))

res.aov4 <- lmer( CO2pHr2 ~ (1|Vial) + Treatment * Sex * Stock, data = resp)
summary(res.aov4)
CO2 <- Anova(res.aov4)
CO2
CO2_1 <- emmeans(res.aov4, ~Sex*Treatment)
plot(CO2_1)
cldCO2 <-cld(CO2_1, Letters = letters, alpha = 0.05)
CO2S=summaryBy(CO2pHr2~Treatment+Sex,data=resp,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
CO2S.summarized=merge(CO2S,cldCO2, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
unique(CO2S.summarized$.group)
CO2S.summarized$group=ifelse(CO2S.summarized$.group=="  b","b", ifelse(CO2S.summarized$.group==" a","a", ifelse(CO2S.summarized$.group==" ab","ab", CO2S.summarized$.group)))


##Lmer by stocks
##NOT SURE IF THIS IS INFORMATIVE BUT 42 AS A SEPARATE STOCK SHOWS SIGNIFICANT 
##DIFFERENCES DUE TO TREATMENT:SEX WHILE 217 DUE TO SEX ONLY. HOWEVER THE
##POSTHOC DOES NOT SHOW ANY SPECIFIC DIFFERENCES (ALL GROUPED AS "A")
res.aov42 <- lmer(RQ ~ (1|Vial) +Treatment * Sex, data = resp42)
Anova(res.aov42)
RQ42 <- emmeans(res.aov42, ~Sex*Treatment)
plot(RQ42)
cldRQ42 <-cld(RQ42, Letters = letters, alpha = 0.05)

res.aov217 <- lmer(RQ ~ (1|Vial) +Treatment * Sex, data = resp217)
Anova(res.aov217)
RQ217 <- emmeans(res.aov217, ~Sex*Treatment)
cldRQ217 <-cld(RQ217, Letters = letters, alpha = 0.05)

res.aov2_42 <- lmer(O2pHr2 ~ (1|Vial) + Treatment * Sex, resp42)
summary(res.aov2_42)
emmeans(res.aov2_42, list(pairwise ~ Treatment*Sex), adjust = "tukey")

res.aov2_217 <- lmer(O2pHr2 ~ (1|Vial) + Treatment * Sex, resp217)
summary(res.aov2_217)
emmeans(res.aov2_217, list(pairwise ~ Treatment*Sex), adjust = "Sidak")

res.aov3_42 <- lmer( CO2pHr2 ~ (1|Vial) + Treatment * Sex, data = resp42)
summary(res.aov3_42)
emmeans(res.aov3_42, list(pairwise ~ Treatment*Sex), adjust = "Sidak")

res.aov3_217 <- lmer( CO2pHr2 ~ (1|Vial) + Treatment * Sex, data = resp217)
summary(res.aov3_217)
emmeans(res.aov3_217, list(pairwise ~ Treatment*Sex), adjust = "Sidak")

fun <- function(x){
  c(m=mean(x), v=var(x), n=length(x), min=min(x), max=max(x))
}
Sum <- summary_by(cbind(O2pHr2,CO2pHr2, Weight, RQ) ~ Treatment + Stock +Sex, data=resp, FUN=fun)
