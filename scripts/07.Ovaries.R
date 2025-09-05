##Oogenesis
library(readr)
library(dplyr)
library(doBy)
library(gridExtra)
library(multcomp)
library(DHARMa)
library(performance)

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

##Compare two rounds (note: Round 2 is ONLY Day 2)
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
Ovary_summary <- summary_by(Ovary ~ Stock + Age + Treatment, data=Ova, FUN=fun)
Oocyte_summary <- summary_by(oocyte_cor ~ Stock + Age + Treatment, data=Ova, FUN=fun)

Total_summary = summary_by(oocyte_cor + Number_Oocyte + Ovariole + Ovary~Ovary_cor+Stock + Age + Treatment,data=Ova, FUN=c(sum,length,min,max,mean))
Total_summary$oocyte_per_ovariole=Total_summary$Number_Oocyte.sum/Total_summary$Ovariole.max


combined_summary2=summary_by(oocyte_per_ovariole+Number_Oocyte.sum+Ovary.sum + Ovariole.sum + Ovariole.mean ~Stock + Age + Treatment,data=Total_summary,FUN = c(sum,mean,length))


combined_summary2=combined_summary2[,c(1:3,5,16,13,10)]
write.csv(combined_summary2,file="output/Ovary_summary.csv",row.names = F,quote = F)

#make sure parameters are coded as factors correctly
Ova$Treatment=as.factor(Ova$Treatment)
Ova$Stock=as.factor(Ova$Stock)
Ova$Age=as.factor(Ova$Age)
Ova$Stage=as.factor(Ova$Stage)

#Load libraries for overdispersion models
library(pscl)
library(glmmTMB)
library(MASS)

#Ideally we would use oocyte per ovariole, but it complicates the model. Let's see if it contains more information or variation than the simpler Oocyte count variable?
plot(Ova$Number_Oocyte~Ova$oocyte_cor)
cor.test(Ova$oocyte_cor,Ova$Number_Oocyte,method="spearman")
#rho=0.96 suggesting these two metrics are highly correlated. Since it is easier to model raw counts rather than counts divided by number of ovarioles, let's take that route!


#Let's try a simple poisson first:
fit_m1 <- glmer(Number_Oocyte ~ Age + Treatment + Stock + Stage + 
                (1|Replicate/Ovary_cor),
                    family = poisson,
                    data = Ova)
#Error when running this model suggests it is inappropriate!
Anova(fit_m1)

# Check for overdispersion
check_overdispersion(fit_m1)
#Overdispersion detected.

#Add observation-level random effect for extra overdispersion
Ova$obs_id <- 1:nrow(Ova)

#Significant dispersion in the data, so switching to Negative Binomial
fit_nb <- glmer.nb(Number_Oocyte ~ Age * Treatment * Stock + 
                            (1|Replicate/Ovary_cor) + (1|obs_id),
                          data = Ova)
#Several errors when running this model too
# Model Table:
Anova(fit_nb)

# Model diagnostics
simulateResiduals(fit_nb, plot = TRUE)
#Still some overdispersion
# Compare to Poisson (just to confirm NB is better)
AIC(fit_m1, fit_nb)  # Lower AIC is better
#df      AIC
#fit_m1 12 8705.989
#fit_nb 22 6747.759
#NB approach is better than simple poisson, but we can do better!

# Let's try zero-inflation approach
#First, check which fixed effect variable accounts for the variation in zeros:
ggplot(Ova, aes(x = Number_Oocyte)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~ Stage)
#Highly variable across stages with Germarium having the most zeros!
#We can incorporate this into the model!

#Fit zero-inflation model
fit_zinb_no_interactions <- glmmTMB(Number_Oocyte ~ Treatment + Stock + Age + Stage +
                      (1|Replicate/Ovary_cor) + (1|obs_id), 
                    ziformula = ~ Stage  ,  # Stage predicts structural zeros
                    family = nbinom2,
                    data = Ova,control = glmmTMBControl(optimizer = nlminb))
#summary(fit_zinb_no_interactions)

simulateResiduals(fit_zinb_no_interactions, plot = TRUE)
#This looks SO much better! No overdispersion. Now we can add in some interactions!

#Adding in interactions; this was done iteratively to avoid overfitting:
fit_zinb_m2 <- glmmTMB(Number_Oocyte ~  Stage + Treatment * Stock * Age  + Stage*Treatment +
                         (1|Replicate/Ovary_cor) + (1|obs_id),
                                    ziformula = ~ Stage  ,  # Stage predicts structural zeros
                                    family = nbinom2,
                                    data = Ova,control = glmmTMBControl(optimizer = nlminb))
#fit_zinb_m2$fit$convergence #should be 0
#summary(fit_zinb_m2)

# Check if this fixes the residual patterns
simulateResiduals(fit_zinb_m2, plot = TRUE)
#Still looks great!

#Check AIC against NB model
AIC(fit_offset_nb,fit_zinb_no_interactions,fit_zinb_m2)
#df      AIC
#fit_offset_nb            22 6710.539
#fit_zinb_no_interactions 19 5511.502
#fit_zinb_m2              39 5497.570
#Model with interactions is the best fit to the data!

#Model Table:
#Anova(fit_zinb_no_interactions)
Anova(fit_zinb_m2)
cld2<-emmeans(fit_zinb_m2, ~Treatment+Age+Stock)
cld3 <- cld(cld2, Letters = letters, alpha = 0.05)

m2_contrast <- emmeans::emmeans(fit_zinb_m2, specs = pairwise ~ Treatment, by=c("Stock","Age"), mode="kenward-roger")
m2_cont <- contrast(m2_contrast, method="pairwise")

m2_posthoc <- as.data.frame(summary(m2_cont))
m2_posthoc$full_sig=ifelse(m2_posthoc$p.value<0.001,"***",ifelse(m2_posthoc$p.value<0.01,"**",ifelse(m2_posthoc$p.value<0.05,"*","")))
m2_posthoc

#Write model tables:
write.csv(Anova(fit_zinb_m2), "output/Oocyte_count_anova_table.csv")
write.csv(m2_posthoc, "output/Oocyte-count_posthoc_table.csv",row.names = F)

# Oocyte Count Plots -------------------------------------------------------------------

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
StagesM <- ggplot(Ova, aes(x=Treatment, y=Number_Oocyte, fill=interaction(Stock, Stage, sep = "."))) +
  geom_col() +
  xlab("Stock") +
  ylab("Number of Oocytes") +
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
t=grid.arrange(StagesM, legend_plot, ncol = 2, widths = c(3, 1))
ggsave("images/Oocyte_count.png",t)

Ooc42 <- subset(Ova, Ova$Stock=="42", na.rm=TRUE, select = c(Treatment, Round, Stock, oocyte_cor, TUNEL_Cell,Number_Oocyte))
Ooc217 <- subset(Ova, Ova$Stock=="217", na.rm=TRUE, select = c(Treatment, Round, Stock, oocyte_cor, TUNEL_Cell,Number_Oocyte))

Oocytes42=ggplot(Ooc42, aes(x=Treatment, y=Number_Oocyte, fill=Treatment)) +
  ylab("Number of oocytes") +
  theme_base()+
  ylim(-0.5,8)+
  #theme(strip.background=element_rect(fill="white"))+
  #theme(strip.text=element_text(color="black", face="bold")) +
  scale_fill_manual(breaks = c("0.5", "1", "2"), values=c("#fdd0a2","#fd8d3c","#a63603")) +
  #geom_text(data = W0_10S.summarized, aes(y = 1.75, x = Age, label = group), position = position_dodge(width = 0.75)) +
  ggtitle("Total number of oocytes for line 42") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=12, face = "bold")) + 
  geom_boxplot()
Oocytes42

Oocytes217=ggplot(Ooc217, aes(x=Treatment, y=Number_Oocyte, fill=Treatment)) +
  ylab("Number of oocytes") +
  theme_base()+
  ylim(-0.5,8)+
  #theme(strip.background=element_rect(fill="white"))+
  #theme(strip.text=element_text(color="black", face="bold")) +
  scale_fill_manual(breaks = c("0.5", "1", "2"), values=c("#c6dbef","#6baed6","#08519c")) +
  #geom_text(data = W0_10S.summarized, aes(y = 1.75, x = Age, label = group), position = position_dodge(width = 0.75)) +
  ggtitle("Total number of oocytes for line 217") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=12, face = "bold")) + 
  geom_boxplot()
Oocytes217


# TUNEL / APOPTOSIS -------------------------------------------------------------------
 
Ova$TUNEL_Ovariole=ifelse(Ova$TUNEL_Cell>0, "Positive", "Negative") 
Ova$TUNEL_Ovariole=as.factor(Ova$TUNEL_Ovariole)

#get number of oocytes that are TUNEL negative
Ova$TUNEL_negative=Ova$Number_Oocyte-Ova$TUNEL_Cell
Ova$Tunel_ova = (Ova$TUNEL_Cell/Ova$Number_Oocyte)*100

#collapse stages into early and late to have fewer combinations with no data
Ova$Stage_grouped <- ifelse(Ova$Stage %in% c("Germarium", "S1-7"), "Early","Late")

# Get mean, median, etc. for each factor combination
TUNEL_summary_stats <- Ova %>%
  dplyr::group_by(Treatment, Stock,Stage_grouped) %>%
  dplyr::summarise(
    mean_tunel = mean(Tunel_ova, na.rm = TRUE),
    num_oocyte = sum(Number_Oocyte, na.rm = TRUE),
    TUN_Pos = sum(TUNEL_Cell, na.rm = TRUE),
    TUN_Neg = sum(TUNEL_negative, na.rm = TRUE),
    .groups = 'drop'
  )
#save summary
write.csv(TUNEL_summary_stats,file="output/TUNEL_summary.csv",row.names = F,quote = F)

ggplot(Ova, aes(x = Treatment, y = Tunel_ova, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(Stage_grouped~Stock) +
  labs(y = "TUNEL Positive (%)", x = "Treatment")

##Here we removed the nested random effect due to overfitting; Age is also not included in any interaction terms since there is some unevenness in the data leading to misleading interactions
TUNELlm <- glmer(cbind(Ova$TUNEL_Cell,Ova$TUNEL_negative) ~ (1|Replicate) + Treatment * Stock * Stage_grouped + Age, data = Ova, family=binomial(link="logit"),control = glmerControl(optimizer = "bobyqa"))
#summary(TUNELlm)
Anova(TUNELlm)

#Write model tables:
write.csv(Anova(TUNELlm), "output/TUNEL_anova_table.csv")

##Letters for posthoc - 
lsmTBn <- emmeans(TUNELlm, ~ Treatment * Stock * Stage_grouped * Age)

cldT <- cld(lsmTBn, Letters = letters)
#cldT

##Odds ratios (leaving Age out since no interaction terms included in model)
fit_contrast_TUN <- emmeans(TUNELlm, specs = pairwise ~ Treatment, by=c("Stock","Stage_grouped"))
fit_contr_TUN <- contrast(fit_contrast_TUN, method="pairwise", type = "response")
pheno_contr_TUN <- as.data.frame(summary(fit_contr_TUN))
pheno_contr_TUN$ooc_sig=ifelse(pheno_contr_TUN$p.value<0.001,"***",ifelse(pheno_contr_TUN$p.value<0.01,"**",ifelse(pheno_contr_TUN$p.value<0.05,"*","")))
#odds_ratioTUN

write.csv(pheno_contr_TUN, "output/TUNEL_posthoc_table.csv",row.names = F)

#Summarize and Plot
TUNEL_summary <- summary_by(Tunel_ova ~ Stock + Treatment + Stage_grouped + Age, data=Ova, FUN=fun)
TUNEL_summary <- as.data.frame(TUNEL_summary)
TUNEL_summary[is.na(TUNEL_summary)] <- 0
merged_TUNEL <- merge(TUNEL_summary, cldT[, c('Treatment', 'Stock', 'Stage_grouped','Age', '.group')], 
                   by = c('Treatment', 'Stock', 'Stage_grouped','Age'))
# Ensure Stock is a factor with correct levels
TUNEL_summary$Stock <- factor(TUNEL_summary$Stock, levels = c("42", "217"))

# Ensure Stage is a factor with the specified order
TUNEL_summary$Stage <- factor(TUNEL_summary$Stage, levels = c("Early", "Late"))

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
    combination = paste(Stock, Stage_grouped, sep = "."),
    color = stock_stage_colors[combination]
  )



###GxE for Oocytes and TUNEL

Ova_summary <- summary_by(oocyte_cor ~ Stock + Treatment + Age, data=Ova, FUN=mean)
Ova_summary$Age_cor=ifelse(Ova_summary$Age=="0D","0",ifelse(Ova_summary$Age=="2D",2,5))
Ova_summary$Age_cor= as.numeric(Ova_summary$Age_cor)


tun_summary <- summary_by(Tunel_ova ~ Stage_grouped + Treatment + Age + Stock, data=Ova, FUN=mean,na.rm)
tun_summary$Age_cor=ifelse(tun_summary$Age=="0D","0",ifelse(tun_summary$Age=="2D",2,5))
tun_summary$Age_cor= as.numeric(tun_summary$Age_cor)
#tun_summary$Stage_cor=ifelse(tun_summary$Stage=="Germarium","1",ifelse(tun_summary$Stage=="S1-7","2",ifelse(tun_summary$Stage=="S8-10","3",ifelse(tun_summary$Stage=="S11","4","5"))))
tun_summary$Stage_cor=ifelse(tun_summary$Stage_grouped=="Early","1","2")
tun_summary$Stage_cor= as.numeric(tun_summary$Stage_cor)

#Plot mean TUNEL cells per day per treatment per stage
tun_gxe_map <- ggplot(tun_summary, aes(x=Stage_cor, y=Tunel_ova.mean,col=Treatment)) +
  facet_grid(Stock~Age_cor) + ylab("Percent Oocytes TUNEL Positive")+
  scale_x_continuous(limits=c(1,2),breaks=c(1,2),labels=c("Early", "Late")) + 
  theme(plot.title = element_blank()) + xlab("Stage of Oogenesis") + 
  theme_base() + 
  geom_line(linewidth = 3, lineend = "round")+ theme(legend.title=element_blank())
tun_gxe_map



Order_Sta=c("Germarium", "S1-7", "S8-10", "S11", "S12-14")
Order_Sta=c("Early", "Late")

##TUNEL Supplemental Figure 9
odds_figure_TUN=ggplot(aes(y=odds.ratio,x=factor(Stage_grouped),group=contrast,col=contrast,shape=contrast),data=pheno_contr_TUN)+
  geom_point(size=3)+ylab("TUNEL Odds Ratio")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+theme(axis.text.x = element_text(angle = 25))+
  #scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+
  geom_line() + geom_errorbar(aes(ymin=odds.ratio-SE,ymax=odds.ratio+SE))+
  facet_grid(Stock~Age)+xlab("Stage of Oogenesis")
odds_figure_TUN
ggsave("images/FigureS9.png",plot=odds_figure_TUN, height=7,width=14)



#Figure 4B without Stages
##Odds ratios (leaving Age out since no interaction terms included in model)
fit_contrast_TUN <- emmeans(TUNELlm, specs = pairwise ~ Treatment, by=c("Stock"))
fit_contr_TUN <- contrast(fit_contrast_TUN, method="pairwise", type = "response")
pheno_contr_TUN <- as.data.frame(summary(fit_contr_TUN))
pheno_contr_TUN$ooc_sig=ifelse(pheno_contr_TUN$p.value<0.001,"***",ifelse(pheno_contr_TUN$p.value<0.01,"**",ifelse(pheno_contr_TUN$p.value<0.05,"*","")))

odds_ratioTUN=pheno_contr_TUN

step1 <- odds_ratioTUN %>%
  mutate(
    Stock = dplyr::recode(Stock,
                          "42" = "DGRP_42",
                          "217" = "DGRP_217"
    ))

step2 <- step1 %>%
  filter(contrast == "Treatment0.5 / Treatment2")

odds_ratioTUN_summary <- step2 %>%
  mutate(
    color_category = case_when(
      Stock == "DGRP_42" & odds.ratio >= 1 ~ "DGRP_42_high",
      Stock == "DGRP_42" & odds.ratio < 1 ~ "DGRP_42_low", 
      Stock == "DGRP_217" & odds.ratio >= 1 ~ "DGRP_217_high",
      Stock == "DGRP_217" & odds.ratio < 1 ~ "DGRP_217_low",
      TRUE ~ "unknown"  # Catch any unexpected cases
    )
  )

color_map <- c(
  "DGRP_42_low" = "#a63603",      # Original dark color for OR >= 1
  "DGRP_42_high" = "#fdd0a2",       # Light color for OR < 1
  "DGRP_217_low" = "#08519c",     # Original dark color for OR >= 1
  "DGRP_217_high" = "#c6dbef"       # Light color for OR < 1
)

odds_figure_TUN_stock <- ggplot(odds_ratioTUN_summary,aes(x = as.factor(Stock), y = odds.ratio, fill = color_category)
) + 
  # Create bars that start from y = 1
  geom_rect(aes(xmin = as.numeric(as.factor(Stock)) - 0.3,
                xmax = as.numeric(as.factor(Stock)) + 0.3,
                ymin = 1,
                ymax = odds.ratio), position = position_dodge()) +
#  geom_col(width = 0.6, position = position_dodge()) +
  geom_errorbar(aes(ymin = odds.ratio - SE, ymax = odds.ratio + SE, color = color_category),
                width = 0.2, size = 1, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = color_map) +
  scale_color_manual(values = color_map) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  ylab("TUNEL Odds Ratio") +
  xlab("") + ylim(c(0,3))+
  theme_base() +
  theme( legend.position = "none")

odds_figure_TUN_stock
ggsave("images/Figure4B.png",plot=odds_figure_TUN_stock, height=7)



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


#difference between 217 and 42 for 0.5x for paper Discussion:
fit_contrast_TUN <- emmeans(TUNELlm, specs = pairwise ~ Stock, by=c("Treatment","Stage_grouped"))
fit_contr_TUN <- contrast(fit_contrast_TUN, method="pairwise", type = "response")
pheno_contr_TUN <- as.data.frame(summary(fit_contr_TUN))
pheno_contr_TUN$ooc_sig=ifelse(pheno_contr_TUN$p.value<0.001,"***",ifelse(pheno_contr_TUN$p.value<0.01,"**",ifelse(pheno_contr_TUN$p.value<0.05,"*","")))
