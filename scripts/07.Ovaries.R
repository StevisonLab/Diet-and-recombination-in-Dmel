##Oogenesis
library(readr)
library(dplyr)
library(doBy)
library(gridExtra)
library(multcomp)

#Load results from model characterization for round2 nutrition
moms <- read_csv("rawdata/P2_Crosses.csv")
F1 <- read_csv("rawdata/P2_F1Collections.csv", 
              col_types = cols(`Virgin Vial Number` = col_number(), 
              F1v = col_number()))

F1$`Virgin Vial Number`=as.numeric(F1$`Virgin Vial Number`)
F1$vial_num=F1$`Virgin Vial Number`
f1 <- F1 %>% distinct(F1v, .keep_all = TRUE)

key <- merge(f1,moms, by.x ="vial_num", by.y = "CrossV")

key1 <- read_csv("rawdata/DAPI.csv")

##Load models of counts
Model1 <-read_csv("rawdata/OvarioleN2.csv")
Model2 <-read_csv("rawdata/OvariolesN1.csv", col_types = cols(Age = col_character()))

##Fill with 0 missing values
Model1[is.na(Model1)] <- 0
Model2[is.na(Model2)] <- 0

##Merge with key
N1Key <- merge(Model1, key, by.x= "Vial", by.y = "F1v", all.x = TRUE)
N2Key <- merge(Model2, key1, by.x= "Vial", by.y = "Vial", all.x = TRUE)

## Create unique replica ID
N1Key$Ovary_cor = paste(N1Key$Replicate, N1Key$Ovary, sep="_")
N2Key$Ovary_cor = paste(N2Key$Replicate, N2Key$Ovary, sep="_")

##Compare two rounds 
N1_D2 <- N1Key %>% filter(Age == "2D")
N2Key$Age=N2Key$Age.x
D2 <- bind_rows(N1_D2, N2Key)
D2$Treatment = as.factor(D2$Treatment)
D2$Round = as.factor(D2$Round)
D2$Stage = as.factor(D2$Stage)
D2$`Male Stock` = as.factor(D2$`Male Stock`)
D2A <- aov(Number_Oocyte ~  Round, data=D2)
summary(D2A)

##Because is not significantly different both rounds will be analyzed as one
N1Key$Stock=N1Key$`Male Stock`
N2Key$vial_num = N2Key$`Mom Vial`
N2Key$Age = N2Key$Age.x
N1 <- N1Key %>% 
  select(1:10, 11,13,21,22)
N2 <- N2Key %>% 
  select(1,2,4:10, 12,13, 16:18)

##Merge 2 rounds
Ovas <- bind_rows(N1, N2)

Ovas$Number_Oocyte=as.numeric(Ovas$Number_Oocyte)

##Sum per stage category per picture

fun <- function(x){
  c(m=mean(x), v=var(x), n=sum(x), min=min(x), max=max(x))}
fun1 <- function(x){
  c(m=mean(x), v=var(x), n=max(x), min=min(x), max=max(x))}

Ovas$oocyte_cor <- Ovas$Number_Oocyte/Ovas$Ovariole
SOocytes <- aggregate(oocyte_cor~Treatment+Age+Stock, FUN=sum,data=Ovas)
hist(Ovas$oocyte_cor)
##Lots of 0 in the data, a negative binomial seems to be a good fit 
Ovas$Treatment=as.factor(Ovas$Treatment)
Ovas$Stock=as.factor(Ovas$Stock)
Ovas$Age=as.factor(Ovas$Age)
Ovas$Stage=as.factor(Ovas$Stage)
Ova <- na.omit(Ovas)
Ovariole_summary <- summary_by(Ovariole ~ Stock + Age + Treatment, data=Ova, FUN=fun)
Ovary_summary <- summary_by(Ovary ~ Stock + Age + Treatment, data=Ova, FUN=fun1)
Oocyte_summary <- summary_by(oocyte_cor ~ Stock + Age + Treatment, data=Ova, FUN=fun)

combined_summary <- Ovariole_summary %>%
  full_join(Ovary_summary, by = c("Stock", "Age", "Treatment")) %>%
  full_join(Oocyte_summary, by = c("Stock", "Age", "Treatment"))

write.csv(combined_summary,file="output/Ovary_summary.csv",row.names = F,quote = F)

#Dispersion test of the data 

dispersion_test <- function(x) 
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  
  cat("Dispersion test of count data:\n",
      length(x), " data points.\n",
      "Mean: ",mean(x),"\n",
      "Variance: ",var(x),"\n",
      "Probability of being drawn from Poisson distribution: ", 
      round(res, 3),"\n", sep = "")
  
  invisible(res)
}
dispersion_test(Ova$oocyte_cor)

##Data underdispersed, with a lot of 0's and skewed to the left
#zero-inflated models are recommended -- Poisson and NB

library(pscl)
library(glmmTMB)
library(MASS)

#Zero-inflated negative binomial
zinb_model <- glmmTMB(oocyte_cor ~ Age * Treatment * Stock * Stage,
                      zi = ~ oocyte_cor,
                      family = nbinom2,
                      data = Ova)
summary(zinb_model)

# Fit the ZIP model with fixed and random effects
zip_model <- glmmTMB(oocyte_cor ~ Age * Treatment * Stock * Stage, 
                     zi = ~ 1,  # You can add predictors for zero inflation here
                     family = poisson, 
                     data = Ova)

summary(zip_model)

##Both models produced NA's in several parameters
#a more simple NB could be a better fit

#fit_Oocyte <- glmer.nb(oocyte_cor ~ (1|Replicate)+(1|Ovary_cor) + (1|Replicate:Ovary_cor)+Age * Treatment * Stock * Stage, data = Ova)
fit_Oocyte=readRDS("output/Oocyte_Model.rds")
#saveRDS(fit_Oocyte,file="output/Oocyte_Model.rds")
summary(fit_Oocyte)
anova_fitOoocyte <- Anova(fit_Oocyte)
anova_fitOoocyte

lsm <- pairs(emmeans(fit_Oocyte, ~ Stock * Treatment * Stage))
summary(fit_Oocyte)$coefficients
summary(lsm, type = "response")

##Odds ratios
fit_contrast_Ooc <- emmeans::emmeans(fit_Oocyte, specs = pairwise ~ Treatment, by=c("Stock","Age"), mode="kenward-roger")
fit_contr_Ooc <- contrast(fit_contrast_Ooc, method="pairwise")
pheno_contr_Ooc <- as.data.frame(summary(fit_contr_Ooc))
pheno_contr_Ooc
odds_ratioOC=pheno_contr_Ooc
odds_ratioOC$ooc_or=exp(odds_ratioOC$estimate)
odds_ratioOC$ooc_sig=ifelse(odds_ratioOC$p.value<0.001,"***",ifelse(odds_ratioOC$p.value<0.01,"**",ifelse(odds_ratioOC$p.value<0.05,"*","")))
odds_ratioOC

##Checking for dispersion in the model
PearsonResiduals <- resid(fit_Oocyte, type = "pearson")
# extract number of cases in model
Cases <- nrow(Ova)
# extract number of predictors (plus intercept)
NumberOfPredictors <- length(fixef(fit_Oocyte)) +1
# calculate overdispersion
Overdispersion <- sum(PearsonResiduals^2) / (Cases-NumberOfPredictors)
# inspect overdispersion
Overdispersion

cld <- cld(lsm, Letters = letters)
cld
fit_Oocyte


OvaS<- OvaKey %>% count(Stock, Age, Treatment)
OvaS1<- Ova1Key %>% count(Stock, Age, Treatment)
write.csv(OvaKey, "output/Nutri_test.csv")
OvaS=summaryBy(file~Age+Stock,data=OvaKey,FUN = sum,na.rm=T, stringsAsFactors = FALSE)


# Ensure Stock is a factor with correct levels
Ova$Stock <- factor(Ova$Stock, levels = c("42", "217"))

# Ensure Stage is a factor with the specified order
Ova$Stage <- factor(Ova$Stage, levels = c("Germarium", "S1-7", "S8-10", "S11", "S12-14"))

# Define colors for each combination of Stock and Stage
stock_stage_colors <- c(
  "42.Germarium" = "#fee391", "42.S1-7" = "#fec44f", "42.S8-10" = "#fe9929", "42.S11" = "#d95f0e", "42.S12-14" = "#993404", 
  "217.Germarium" = "#d0d1e6", "217.S1-7" = "#a6bddb", "217.S8-10" = "#74a9cf", "217.S11" = "#2b8cbe", "217.S12-14" = "#045a8d"
)

# Plot data (Supplemental Figure 4)
StagesM <- ggplot(Ova, aes(x=Treatment, y=oocyte_cor, fill=interaction(Stock, Stage, sep = "."))) +
  geom_col() +
  xlab("Stock") +
  ylab("Oocyte per ovariole") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = stock_stage_colors, guide = "none") +
  facet_wrap(~Stock*Age)

# Custom legend
legend_data <- data.frame(
  Stock = factor(c(rep("42", 5), rep("217", 5)), levels = c("42", "217")),
  Stage = factor(rep(c("Germarium", "S1-7", "S8-10", "S11", "S12-14"), 2), 
                 levels = c("Germarium", "S1-7", "S8-10", "S11", "S12-14")),
  Color = c("#fee391", "#fec44f", "#fe9929", "#d95f0e", "#993404", 
            "#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")
)

legend_plot <- ggplot(legend_data, aes(x = 1, y = Stage, fill = Color)) +
  geom_tile() +
  facet_wrap(~Stock) +
  scale_fill_identity() +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
grid.arrange(StagesM, legend_plot, ncol = 2, widths = c(3, 1))

Ooc42 <- subset(Ova, Ova$Stock=="42", na.rm=TRUE, select = c(Treatment, Round, Stock, oocyte_cor, TUNEL_Cell))
Ooc217 <- subset(Ova, Ova$Stock=="217", na.rm=TRUE, select = c(Treatment, Round, Stock, oocyte_cor, TUNEL_Cell))

Oocytes42=ggplot(Ooc42, aes(x=Treatment, y=oocyte_cor, fill=Treatment)) +
  ylab("Number of oocytes") +
  theme_base()+
  ylim(-0.5,8)+
  #theme(strip.background=element_rect(fill="white"))+
  #theme(strip.text=element_text(color="black", face="bold")) +
  scale_fill_manual(breaks = c("0.5", "1", "2"), values=c("#fdd0a2","#fd8d3c","#a63603")) +
  #geom_text(data = W0_10S.summarized, aes(y = 1.75, x = Age, label = group), position = position_dodge(width = 0.75)) +
  ggtitle("Total number of oocytes per ovary per ovariole for line 42") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=12, face = "bold")) + 
  geom_boxplot()
Oocytes42

Oocytes217=ggplot(Ooc217, aes(x=Treatment, y=oocyte_cor, fill=Treatment)) +
  ylab("Number of oocytes") +
  theme_base()+
  ylim(-0.5,8)+
  #theme(strip.background=element_rect(fill="white"))+
  #theme(strip.text=element_text(color="black", face="bold")) +
  scale_fill_manual(breaks = c("0.5", "1", "2"), values=c("#c6dbef","#6baed6","#08519c")) +
  #geom_text(data = W0_10S.summarized, aes(y = 1.75, x = Age, label = group), position = position_dodge(width = 0.75)) +
  ggtitle("Total number of oocytes per ovary per ovariole for line 217") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=12, face = "bold")) + 
  geom_boxplot()
Oocytes217


# TUNEL / APOPTOSIS -------------------------------------------------------------------
 
Ova$TUNEL_Ovariole=ifelse(Ova$TUNEL_Cell>0, "Positive", "Negative") 
Ova$TUNEL_Ovariole=as.factor(Ova$TUNEL_Ovariole)

##We need to compare which model is more appropriate

##We might need to just remove the random effects
TUNELlm <- glm(response ~ Treatment * Age * Stage * Stock, data = Ova)
summary(TUNELlm)
Anova(TUNELlm)
##Letters for posthoc - 
lsmTBn <- emmeans(TUNELlm, ~ Age * Treatment * Stock * Stage)

cldT <- cld(lsmTBn, Letters = letters)
cldT

Ova$Tunel_ova = (Ova$TUNEL_Cell/Ova$Number_Oocyte)*100
##Odds ratios
fit_contrast_TUN <- emmeans::emmeans(TUNELlm, ~ Treatment, by=c("Stock","Age", "Stage"), mode="kenward-roger")
fit_contr_TUN <- contrast(fit_contrast_TUN, method="pairwise")
pheno_contr_TUN <- as.data.frame(summary(fit_contr_TUN))
pheno_contr_TUN

odds_ratioTUN=pheno_contr_TUN
odds_ratioTUN$tun_or=exp(odds_ratioTUN$estimate)
odds_ratioTUN$ooc_sig=ifelse(odds_ratioTUN$p.value<0.001,"***",ifelse(odds_ratioTUN$p.value<0.01,"**",ifelse(odds_ratioTUN$p.value<0.05,"*","")))
odds_ratioTUN

TUNEL_summary <- summary_by(Tunel_ova ~ Stock + Age + Treatment + Stage, data=Ova, FUN=fun)
TUNEL_summary <- as.data.frame(TUNEL_summary)
TUNEL_summary[is.na(TUNEL_summary)] <- 0
merged_TUNEL <- merge(TUNEL_summary, cldT[, c('Age', 'Treatment', 'Stock', 'Stage', '.group')], 
                   by = c('Age', 'Treatment', 'Stock', 'Stage'))
# Ensure Stock is a factor with correct levels
TUNEL_summary$Stock <- factor(TUNEL_summary$Stock, levels = c("42", "217"))

# Ensure Stage is a factor with the specified order
TUNEL_summary$Stage <- factor(TUNEL_summary$Stage, levels = c("Germarium", "S1-7", "S8-10", "S11", "S12-14"))

# Define colors for each combination of Stock and Stage
stock_stage_colors <- c(
  "42.Germarium" = "#fee391", "42.S1-7" = "#fec44f", "42.S8-10" = "#fe9929", "42.S11" = "#d95f0e", "42.S12-14" = "#993404", 
  "217.Germarium" = "#d0d1e6", "217.S1-7" = "#a6bddb", "217.S8-10" = "#74a9cf", "217.S11" = "#2b8cbe", "217.S12-14" = "#045a8d"
)
# Extract significant contrasts (p.value < 0.05)
significant_contrasts <- function(data) {
  data <- data[data$p.value < 0.05, ]
  return(data)
}

significant_lines = significant_contrasts(pheno_contr_TUN)
significance_df <- pheno_contr_TUN %>%
  filter(p.value < 0.05) %>%  # Filter for significant results
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    combination = paste(Stock, Stage, sep = "."),
    color = stock_stage_colors[combination]
  )



###GxE for Oocytes and TUNEL

Ova_summary <- summary_by(oocyte_cor ~ Stock + Treatment + Age, data=Ova, FUN=mean)
Ova_summary$Age_cor=ifelse(Ova_summary$Age=="0D","0",ifelse(Ova_summary$Age=="2D",2,5))
Ova_summary$Age_cor= as.numeric(Ova_summary$Age_cor)


tun_summary <- summary_by(TUNEL_Cell ~ Stage + Treatment + Age, data=Ova, FUN=mean)
tun_summary$Age_cor=ifelse(tun_summary$Age=="0D","0",ifelse(tun_summary$Age=="2D",2,5))
tun_summary$Age_cor= as.numeric(tun_summary$Age_cor)
tun_summary$Stage_cor=ifelse(tun_summary$Stage=="Germarium","1",ifelse(tun_summary$Stage=="S1-7","2",ifelse(tun_summary$Stage=="S8-10","3",ifelse(tun_summary$Stage=="S11","4","5"))))
tun_summary$Stage_cor= as.numeric(tun_summary$Stage_cor)

#Plot mean TUNEL cells per day per treatment per stage
tun_gxe_map <- ggplot(tun_summary, aes(x=Stage_cor, y=TUNEL_Cell.mean,col=Treatment)) +
  facet_wrap(~Age_cor) +
  scale_x_continuous(limits=c(1,5),breaks=c(1,2,3,4,5),labels=c("Germarium","S1-7","S8-10","S11","S12-14")) + 
  theme(plot.title = element_blank()) +
  theme_base()+ 
  geom_line(linewidth = 3, lineend = "round")+ theme(legend.title=element_blank())
tun_gxe_map



Order_Sta=c("Germarium", "S1-7", "S8-10", "S11", "S12-14")

##TUNEL Supplemental Figure 9
odds_figure_TUN=ggplot(aes(y=tun_or,x=factor(Stage,levels = Order_Sta),group=contrast,col=contrast,shape=contrast),data=odds_ratioTUN)+
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+theme(axis.text.x = element_text(angle = 25))+
  #scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+
  geom_line()+geom_errorbar(aes(ymin=tun_or-SE,ymax=tun_or+SE))+
  facet_grid(Stock ~ Age)
odds_figure_TUN



#Figure 4B without Stages
fit_contrast_TUN <- emmeans::emmeans(TUNELlm, ~ Treatment, by=c("Stock"), mode="kenward-roger")
fit_contr_TUN <- contrast(fit_contrast_TUN, method="pairwise")
pheno_contr_TUN <- as.data.frame(summary(fit_contr_TUN))
pheno_contr_TUN

odds_ratioTUN=pheno_contr_TUN
odds_ratioTUN$tun_or=exp(odds_ratioTUN$estimate)
odds_ratioTUN$ooc_sig=ifelse(odds_ratioTUN$p.value<0.001,"***",ifelse(odds_ratioTUN$p.value<0.01,"**",ifelse(odds_ratioTUN$p.value<0.05,"*","")))
odds_ratioTUN


step1 <- odds_ratioTUN %>%
  mutate(
    Stock = dplyr::recode(Stock,
                          "42" = "DGRP_42",
                          "217" = "DGRP_217"
    ))

step2 <- step1 %>%
  filter(contrast == "Treatment0.5 - Treatment2")

odds_ratioTUN_summary <- step2 %>%
  mutate(
    color_category = case_when(
      Stock == "DGRP_42" & tun_or >= 1 ~ "DGRP_42_high",
      Stock == "DGRP_42" & tun_or < 1 ~ "DGRP_42_low", 
      Stock == "DGRP_217" & tun_or >= 1 ~ "DGRP_217_high",
      Stock == "DGRP_217" & tun_or < 1 ~ "DGRP_217_low",
      TRUE ~ "unknown"  # Catch any unexpected cases
    )
  )

color_map <- c(
  "DGRP_42_low" = "#a63603",      # Original dark color for OR >= 1
  "DGRP_42_high" = "#fdd0a2",       # Light color for OR < 1
  "DGRP_217_low" = "#08519c",     # Original dark color for OR >= 1
  "DGRP_217_high" = "#c6dbef"       # Light color for OR < 1
)

odds_figure_TUN_stock <- ggplot(
  odds_ratioTUN_summary,
  aes(x = as.factor(Stock), y = tun_or, fill = color_category)
) + 
  # Create bars that start from y = 1
  geom_rect(aes(xmin = as.numeric(as.factor(Stock)) - 0.3,
                xmax = as.numeric(as.factor(Stock)) + 0.3,
                ymin = 1,
                ymax = tun_or), position = position_dodge()) +
#  geom_col(width = 0.6, position = position_dodge()) +
#  geom_errorbar(aes(ymin = tun_or - SE, ymax = tun_or + SE, color = color_category),
#                width = 0.2, size = 1, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  ylab("Odds Ratio from TUNEL assay") +
  xlab("") + ylim(c(0.96,1.04))+
  theme_base() +
  theme( legend.position = "none")

odds_figure_TUN_stock

# Figure 4B with Stages
fit_contrast_TUN <- emmeans::emmeans(TUNELlm, ~ Treatment, by=c("Stock","Stage"), mode="kenward-roger")
fit_contr_TUN <- contrast(fit_contrast_TUN, method="pairwise")
pheno_contr_TUN <- as.data.frame(summary(fit_contr_TUN))
pheno_contr_TUN

odds_ratioTUN=pheno_contr_TUN
odds_ratioTUN$tun_or=exp(odds_ratioTUN$estimate)
odds_ratioTUN$ooc_sig=ifelse(odds_ratioTUN$p.value<0.001,"***",ifelse(odds_ratioTUN$p.value<0.01,"**",ifelse(odds_ratioTUN$p.value<0.05,"*","")))
odds_ratioTUN


step1 <- odds_ratioTUN %>%
  mutate(
    Stock = dplyr::recode(Stock,
                          "42" = "DGRP_42",
                          "217" = "DGRP_217"
    ))

step2 <- step1 %>%
  filter(contrast == "Treatment0.5 - Treatment2")

odds_ratioTUN_summary <- step2 %>%
  mutate(
    color_category = case_when(
      Stock == "DGRP_42" & tun_or >= 1 ~ "DGRP_42_high",
      Stock == "DGRP_42" & tun_or < 1 ~ "DGRP_42_low", 
      Stock == "DGRP_217" & tun_or >= 1 ~ "DGRP_217_high",
      Stock == "DGRP_217" & tun_or < 1 ~ "DGRP_217_low",
      TRUE ~ "unknown"  # Catch any unexpected cases
    )
  )

color_map <- c(
  "DGRP_42_low" = "#a63603",      # Original dark color for OR >= 1
  "DGRP_42_high" = "#fdd0a2",       # Light color for OR < 1
  "DGRP_217_low" = "#08519c",     # Original dark color for OR >= 1
  "DGRP_217_high" = "#c6dbef"       # Light color for OR < 1
)

odds_figure_TUN_stages <- ggplot(
  odds_ratioTUN_summary,
  aes(x = as.factor(Stage), y = tun_or, fill = color_category,group=Stock)
) + #facet_wrap(~Stock) +
  # Create bars that start from y = 1
  geom_rect(aes(xmin = as.numeric(as.factor(Stage)) - 0.3,
                xmax = as.numeric(as.factor(Stage)) + 0.3,
                ymin = 1,
                ymax = tun_or), position = position_dodge()) +
  #  geom_col(width = 0.6, position = position_dodge()) +
#  geom_errorbar(aes(ymin = tun_or - SE, ymax = tun_or + SE, color = color_category),
 #               width = 0.2, size = 1, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  ylab("Odds Ratio from TUNEL assay") +
  xlab("Stage of Oogenesis") +
  theme_base() +
  theme( legend.position = "none")

odds_figure_TUN_stages
ggsave("images/Figure4B.png",plot=odds_figure_TUN_stages, height=7)


T42 <- subset(Ova, Ova$Stock=="42", na.rm=TRUE, select = c(Treatment, Round, Stock, oocyte_cor, Tunel_ova))
T217 <- subset(Ova, Ova$Stock=="217", na.rm=TRUE, select = c(Treatment, Round, Stock, oocyte_cor, Tunel_ova))
Ova[is.nan(Ova)] <- 0
TP42=ggplot(T42, aes(x=Treatment, y=Tunel_ova, fill=Treatment)) +
  ylab("TUNEL positive oocytes (%)") +
  theme_base()+
  ylim(-0.5,70)+
  #theme(strip.background=element_rect(fill="white"))+
  #theme(strip.text=element_text(color="black", face="bold")) +
  scale_fill_manual(breaks = c("0.5", "1", "2"), values=c("#fdd0a2","#fd8d3c","#a63603")) +
  #geom_text(data = W0_10S.summarized, aes(y = 1.75, x = Age, label = group), position = position_dodge(width = 0.75)) +
  ggtitle("Percentage of TUNEL positive oocytes per ovary per ovariole for line 42") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=12, face = "bold")) + 
  geom_boxplot()
TP42


TP217=ggplot(T217, aes(x=Treatment, y=Tunel_ova, fill=Treatment)) +
  ylab("TUNEL positive oocytes (%)") +
  theme_base()+
  ylim(-0.5,70)+
  #theme(strip.background=element_rect(fill="white"))+
  #theme(strip.text=element_text(color="black", face="bold")) +
  scale_fill_manual(breaks = c("0.5", "1", "2"), values=c("#c6dbef","#6baed6","#08519c")) +
  #geom_text(data = W0_10S.summarized, aes(y = 1.75, x = Age, label = group), position = position_dodge(width = 0.75)) +
  ggtitle("Percentage of TUNEL positive oocytes per ovary per ovariole for line 217") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=12, face = "bold")) + 
  geom_boxplot()
TP217

