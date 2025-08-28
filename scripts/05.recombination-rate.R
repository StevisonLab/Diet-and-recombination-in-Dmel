# Recombination rate analysis

kosambi_correction <- function(input) {
  (100/4)*(log((1+(2*input))/(1-(2*input))))
}

recomb <- by_vialday

recomb$nco_inds <- rowSums(recomb[names(recomb) %in% haps_nco], )
recomb$sco_inds <- rowSums(recomb[names(recomb) %in% haps_sco], )
recomb$dco_inds <- rowSums(recomb[names(recomb) %in% haps_dco], )
recomb$tco_inds <- rowSums(recomb[names(recomb) %in% haps_tco], )
recomb$co_inds <- recomb$sco_inds + 2*recomb$dco_inds + 3*recomb$tco_inds
recomb$num_M <- rowSums(recomb[names(recomb) %in% haplotypes], )

recomb$recomb_rate <- recomb$co_inds / recomb$num_M
#recomb$kosambi <- kosambi_correction(recomb$recomb_rate)

recomb <- na.omit(recomb)

## Rough recombination rate by interval for validation purposes
recomb$ycv_count <-rowSums(recomb[names(recomb) %in% ycv_haps], )
recomb$cvv_count <-rowSums(recomb[names(recomb) %in% cvv_haps], )
recomb$vf_count <-rowSums(recomb[names(recomb) %in% vf_haps], )

recomb$ycv_rr <- recomb$ycv_count / recomb$num_M
recomb$cvv_rr <-recomb$cvv_count / recomb$num_M
recomb$vf_rr <-recomb$vf_count / recomb$num_M

total_ycv_count <- sum(recomb[names(recomb) %in% ycv_haps])
total_cvv_count <- sum(recomb[names(recomb) %in% cvv_haps])
total_vf_count <- sum(recomb[names(recomb) %in% vf_haps])
total_num_M <- sum(recomb["num_M"])

total_ycv_rr <- total_ycv_count / total_num_M
total_cvv_rr <-total_cvv_count / total_num_M
total_vf_rr <-total_vf_count / total_num_M

recomb_intervals <- data.frame(interval = c("y-cv", "cv-v", "v-f"),
                               observed = c(total_ycv_rr, total_cvv_rr, total_vf_rr) * 100,
                               expected = c(13.7, 19.3, 23.7))
recomb_intervals$kosambi=kosambi_correction(recomb_intervals$observed/100)



# Run COI script
source("scripts/05b.COI.R")

## Remove haplotype columns
recomb <- recomb[,!names(recomb) %in% haplotypes]

recomb_rate=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M~MaternalVial+vial_letter+PaternalStock+Treatment,data=recomb,FUN=sum,na.rm=T)
recomb_rate=recomb_rate[recomb_rate$num_M.sum>=5,]
recomb_rate$ycv_rr=kosambi_correction(recomb_rate$ycv_count.sum/recomb_rate$num_M.sum)
recomb_rate$cvv_rr=kosambi_correction(recomb_rate$cvv_count.sum/recomb_rate$num_M.sum)
recomb_rate$vf_rr=kosambi_correction(recomb_rate$vf_count.sum/recomb_rate$num_M.sum)
recomb_rate$recomb_rate=rowSums(recomb_rate[,c("ycv_rr","cvv_rr","vf_rr")],na.rm=TRUE)
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$PaternalStock), mean)
names(recomb_summary)[1:2] <- c("Treatment", "PaternalStock")

recomb_summary$treat=ifelse(recomb_summary$Treatment=="2x",2,ifelse(recomb_summary$Treatment=="1x",1,0.5))




## Stats and figures

## FULL interval recomb model

recomb_rate$nco_inds=recomb_rate$num_M.sum-recomb_rate$co_inds.sum
fit_full <- glmer(cbind(co_inds.sum,nco_inds) ~ (1|MaternalVial) + Treatment * PaternalStock * vial_letter, data=recomb_rate,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_full <- Anova(fit_full, test = "Chisq")
anova_full

write.csv(anova_full, "output/recomb-stats_anova-full.csv")

## Post hoc test for full
fit_contrast_full <- emmeans::emmeans(fit_full, specs = pairwise ~ Treatment, by=c("PaternalStock","vial_letter"), mode="kenward-roger")
fit_contr_full <- contrast(fit_contrast_full, method="pairwise")

pheno_contr_full <- as.data.frame(summary(fit_contr_full))
pheno_contr_full$full_sig=ifelse(pheno_contr_full$p.value<0.001,"***",ifelse(pheno_contr_full$p.value<0.01,"**",ifelse(pheno_contr_full$p.value<0.05,"*","")))
pheno_contr_full

write.csv(pheno_contr_full, "output/recomb-stats_anova-full_posthoc.csv")

odds_ratio=subset(pheno_contr_full,pheno_contr_full$PaternalStock=="42")
odds_ratio$full_or=exp(odds_ratio$estimate)
odds_ratio$full_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

pdf("images/full_odds_42.pdf")

odds_figure_full=ggplot(aes(y=full_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("full odds ratio 42")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=full_or-SE,ymax=full_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$full_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$full_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$full_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$full_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$full_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$full_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$full_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$full_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$full_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[12],color="#7570b3",size=10)
odds_figure_full
dev.off()

#repeat for 217
odds_ratio=subset(pheno_contr_full,pheno_contr_full$PaternalStock=="217")
odds_ratio$full_or=exp(odds_ratio$estimate)
odds_ratio$full_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

pdf("images/full_odds_217.pdf")
odds_figure_full=ggplot(aes(y=full_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("full odds ratio 217")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,2.5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=full_or-SE,ymax=full_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$full_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$full_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$full_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$full_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$full_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$full_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$full_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$full_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$full_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[12],color="#7570b3",size=10)
odds_figure_full
dev.off()

#Boxplot for full
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$PaternalStock,recomb_rate$vial_letter), mean)
names(recomb_summary)[1:3] <- c("Treatment", "PaternalStock","Day")

recomb_summary$treat=ifelse(recomb_summary$Treatment=="2x",2,ifelse(recomb_summary$Treatment=="1x",1,0.5))

pdf("images/diet_gxe_byDay.pdf")
diet_gxe_map <- ggplot(recomb_summary, aes(x=as.numeric(treat), y=recomb_rate*100,col=factor(PaternalStock)))  + facet_wrap(~Day) +
  xlab("Caloric Density") + ylab("Percent Recombination") + ggtitle("GxE diet") + 
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) +  
  theme_base()+ 
  geom_line() + theme(legend.title=element_blank())


diet_gxe_map
dev.off()


D_data <- recomb_rate[recomb_rate$vial_letter=="D",]
D_summary=aggregate(D_data[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(D_data$Treatment, D_data$PaternalStock), mean)
names(D_summary)[1:2] <- c("Treatment", "PaternalStock")

D_summary$treat=ifelse(D_summary$Treatment=="2x",2,ifelse(D_summary$Treatment=="1x",1,0.5))

diet_gxe_map_D <- ggplot(D_summary, aes(x=as.numeric(treat), y=recomb_rate*100,col=factor(PaternalStock)))  +
  xlab("Caloric Density") + ylab("Percent Recombination") + ggtitle("GxE diet") + 
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) +  
  theme_base()+ 
  geom_line() + theme(legend.title=element_blank()) 
#annotate(geom="text", x=2, y=2.5, label=odds_ratio$full_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$full_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$full_sig[12],color="#7570b3",size=10)


diet_gxe_map_D


## y-cv interval recomb model
fit_ycv <- glmer(cbind(ycv_count,nco_inds) ~ (1|MaternalVial) + Treatment * PaternalStock * vial_letter, data=recomb,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_ycv <- Anova(fit_ycv, test = "Chisq")
anova_ycv
write.csv(anova_ycv, "output/recomb-stats_anova-ycv_new.csv")


## Post hoc test for y-cv
fit_contrast_ycv <- emmeans::emmeans(fit_ycv, specs = pairwise ~ Treatment, by=c("PaternalStock","vial_letter"), mode="kenward-roger")
fit_contr_ycv <- contrast(fit_contrast_ycv, method="pairwise")

pheno_contr_ycv <- as.data.frame(summary(fit_contr_ycv))
pheno_contr_ycv

write.csv(pheno_contr_ycv, "output/recomb-stats_anova-ycv_posthoc.csv")


## cv-v interval recomb model
fit_cvv <- glmer(cbind(cvv_count,nco_inds) ~ (1|MaternalVial) + Treatment * PaternalStock * vial_letter, data=recomb,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_cvv <- Anova(fit_cvv, test = "Chisq")
anova_cvv
write.csv(anova_cvv, "output/recomb-stats_anova-cvv.csv")


## Post hoc test for cv-v
fit_contrast_cvv <- emmeans::emmeans(fit_cvv, specs = pairwise ~ Treatment, by=c("PaternalStock","vial_letter"), mode="kenward-roger")
fit_contr_cvv <- contrast(fit_contrast_cvv, method="pairwise")

pheno_contr_cvv <- as.data.frame(summary(fit_contr_cvv))
pheno_contr_cvv

write.csv(pheno_contr_cvv, "output/recomb-stats_anova-cvv_posthoc.csv")


## v-f interval recomb model
fit_vf <- glmer(cbind(vf_count,nco_inds) ~ (1|MaternalVial) + Treatment * PaternalStock * vial_letter, data=recomb,family=binomial(link="logit"),control = glmerControl(
  optimizer ='optimx', optCtrl=list(method='nlminb')))
anova_vf <- Anova(fit_vf, test = "Chisq")
anova_vf

write.csv(anova_vf, "output/recomb-stats_anova-vf.csv")


## Post hoc test for v-f
fit_contrast_vf <- emmeans::emmeans(fit_vf, specs = pairwise ~ Treatment, by=c("PaternalStock","vial_letter"), mode="kenward-roger")
fit_contr_vf <- contrast(fit_contrast_vf, method="pairwise")

pheno_contr_vf <- as.data.frame(summary(fit_contr_vf))
pheno_contr_vf

write.csv(pheno_contr_vf, "output/recomb-stats_anova-vf_posthoc.csv")


# Summary table for paper:
recomb_rate=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M~vial_letter+PaternalStock+Treatment,data=recomb,FUN=sum,na.rm=T)
recomb_rate=recomb_rate[recomb_rate$num_M.sum>=5,]
recomb_rate$ycv_rr=kosambi_correction(recomb_rate$ycv_count.sum/recomb_rate$num_M.sum)
recomb_rate$cvv_rr=kosambi_correction(recomb_rate$cvv_count.sum/recomb_rate$num_M.sum)
recomb_rate$vf_rr=kosambi_correction(recomb_rate$vf_count.sum/recomb_rate$num_M.sum)
recomb_rate$recomb_rate=rowSums(recomb_rate[,c("ycv_rr","cvv_rr","vf_rr")],na.rm=TRUE)
recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$PaternalStock), mean)
names(recomb_summary)[1:2] <- c("Treatment", "PaternalStock")

#recomb_summary <- aggregate(recomb_rate[, c("recomb_rate", "ycv_rr", "cvv_rr", "vf_rr")], by = list(recomb_rate$Treatment, recomb_rate$PaternalStock), sd)
#names(recomb_summary)[1:2] <- c("Treatment", "PaternalStock")

# Differences in recomb rate between diets for strain 42
sum(recomb_summary$ycv_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="0.5x"]+recomb_summary$cvv_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="0.5x"]+recomb_summary$vf_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="0.5x"])
55.23938-51.92768 #low calorie to control
sum(recomb_summary$ycv_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="1x"]+recomb_summary$cvv_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="1x"]+recomb_summary$vf_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="1x"])
51.92768-44.206 #high calorie to control
sum(recomb_summary$ycv_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="2x"]+recomb_summary$cvv_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="2x"]+recomb_summary$vf_rr[recomb_summary$PaternalStock=="42" & recomb_summary$Treatment=="2x"])
55.23938-44.206 #low calorie to high calorie

# Differences in recomb rate between diets for strain 217
sum(recomb_summary$ycv_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="0.5x"]+recomb_summary$cvv_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="0.5x"]+recomb_summary$vf_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="0.5x"])
49.83412-46.28678 #low calorie to control
sum(recomb_summary$ycv_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="1x"]+recomb_summary$cvv_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="1x"]+recomb_summary$vf_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="1x"])
46.28678-45.80224 #high calorie to control
sum(recomb_summary$ycv_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="2x"]+recomb_summary$cvv_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="2x"]+recomb_summary$vf_rr[recomb_summary$PaternalStock=="217" & recomb_summary$Treatment=="2x"])
49.83412-45.80224 #low calorie to high calorie



## ODDS-ratio figures


#extract odds ratio and standard error for 42 ONLY y-cv
odds_ratio=subset(pheno_contr_ycv,pheno_contr_ycv$PaternalStock=="42")
odds_ratio$ycv_or=exp(odds_ratio$estimate)
odds_ratio$ycv_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

#ODDs plot for 42 for y-cv

pdf("images/ycv_odds_42.pdf")

odds_figure_ycv=ggplot(aes(y=ycv_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("ycv odds ratio 42")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=ycv_or-SE,ymax=ycv_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$ycv_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$ycv_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$ycv_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$ycv_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$ycv_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$ycv_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$ycv_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$ycv_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$ycv_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$ycv_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$ycv_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$ycv_sig[12],color="#7570b3",size=10)
odds_figure_ycv
dev.off()


#extract odds ratio and standard error for 217 ONLY y-cv
odds_ratio=subset(pheno_contr_ycv,pheno_contr_ycv$PaternalStock=="217")
odds_ratio$ycv_or=exp(odds_ratio$estimate)
odds_ratio$ycv_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

#ODDs plot for 217 for y-cv

pdf("images/ycv_odds_217.pdf")

odds_figure_ycv=ggplot(aes(y=ycv_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("ycv odds ratio 217")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=ycv_or-SE,ymax=ycv_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$ycv_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$ycv_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$ycv_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$ycv_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$ycv_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$ycv_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$ycv_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$ycv_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$ycv_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$ycv_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$ycv_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$ycv_sig[12],color="#7570b3",size=10)
odds_figure_ycv
dev.off()


#extract odds ratio and standard error for 42 ONLY cv-v
odds_ratio=subset(pheno_contr_cvv,pheno_contr_cvv$PaternalStock=="42")
odds_ratio$cvv_or=exp(odds_ratio$estimate)
odds_ratio$cvv_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

#ODDs plot for 42 for cv-v

pdf("images/cvv_odds_42.pdf")

odds_figure_cvv=ggplot(aes(y=cvv_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("cvv odds ratio 42")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=cvv_or-SE,ymax=cvv_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$cvv_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$cvv_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$cvv_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$cvv_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$cvv_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$cvv_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$cvv_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$cvv_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$cvv_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$cvv_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$cvv_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$cvv_sig[12],color="#7570b3",size=10)
odds_figure_cvv
dev.off()


#extract odds ratio and standard error for 217 ONLY cv-v
odds_ratio=subset(pheno_contr_cvv,pheno_contr_cvv$PaternalStock=="217")
odds_ratio$cvv_or=exp(odds_ratio$estimate)
odds_ratio$cvv_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

#ODDs plot for 217 for cv-v

pdf("images/cvv_odds_217.pdf")

odds_figure_cvv=ggplot(aes(y=cvv_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("cvv odds ratio 217")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=cvv_or-SE,ymax=cvv_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$cvv_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$cvv_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$cvv_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$cvv_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$cvv_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$cvv_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$cvv_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$cvv_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$cvv_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$cvv_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$cvv_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$cvv_sig[12],color="#7570b3",size=10)
odds_figure_cvv
dev.off()


#extract odds ratio and standard error for 42 ONLY v-f
odds_ratio=subset(pheno_contr_vf,pheno_contr_vf$PaternalStock=="42")
odds_ratio$vf_or=exp(odds_ratio$estimate)
odds_ratio$vf_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

#ODDs plot for 42 for v-f

pdf("images/vf_odds_42.pdf")

odds_figure_vf=ggplot(aes(y=vf_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("vf odds ratio 42")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=vf_or-SE,ymax=vf_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$vf_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$vf_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$vf_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$vf_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$vf_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$vf_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$vf_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$vf_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$vf_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$vf_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$vf_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$vf_sig[12],color="#7570b3",size=10)
odds_figure_vf
dev.off()


#extract odds ratio and standard error for 217 ONLY v-f
odds_ratio=subset(pheno_contr_vf,pheno_contr_vf$PaternalStock=="217")
odds_ratio$vf_or=exp(odds_ratio$estimate)
odds_ratio$vf_sig=ifelse(odds_ratio$p.value<0.001,"***",ifelse(odds_ratio$p.value<0.01,"**",ifelse(odds_ratio$p.value<0.05,"*","")))
odds_ratio

#ODDs plot for 217 for v-f

pdf("images/vf_odds_217.pdf")

odds_figure_vf=ggplot(aes(y=vf_or,x=vial_letter,group=contrast,col=contrast,shape=contrast),data=odds_ratio)+ggtitle("vf odds ratio 217")+scale_colour_manual(values=c("#1b9e77","#d95f02","#7570b3"))+ 
  geom_point(size=3)+ylab("Odds Ratios")+theme_base()+geom_hline(yintercept = 1,linetype="dashed",color="grey")+ylim(0.5,5)+theme(axis.text.x = element_text(angle = 45))+
  scale_x_discrete(name="Days post-mating",labels=c("1-2","3-4","5-6","7-8"))+geom_line()+geom_errorbar(aes(ymin=vf_or-SE,ymax=vf_or+SE))+
  annotate(geom="text", x=1, y=2, label=odds_ratio$vf_sig[1],color="#1b9e77",size=10)+annotate(geom="text", x=2, y=2, label=odds_ratio$vf_sig[4],color="#1b9e77",size=10)+annotate(geom="text", x=3, y=2, label=odds_ratio$vf_sig[7],color="#1b9e77",size=10)+annotate(geom="text", x=4, y=2, label=odds_ratio$vf_sig[10],color="#1b9e77",size=10)+
  annotate(geom="text", x=1, y=2.25, label=odds_ratio$vf_sig[2],color="#d95f02",size=10)+annotate(geom="text", x=2, y=2.25, label=odds_ratio$vf_sig[5],color="#d95f02",size=10)+annotate(geom="text", x=3, y=2.25, label=odds_ratio$vf_sig[8],color="#d95f02",size=10)+annotate(geom="text", x=4, y=2.25, label=odds_ratio$vf_sig[11],color="#d95f02",size=10)+
  annotate(geom="text", x=1, y=2.5, label=odds_ratio$vf_sig[3],color="#7570b3",size=10)+annotate(geom="text", x=2, y=2.5, label=odds_ratio$vf_sig[6],color="#7570b3",size=10)+annotate(geom="text", x=3, y=2.5, label=odds_ratio$vf_sig[9],color="#7570b3",size=10)+annotate(geom="text", x=4, y=2.5, label=odds_ratio$vf_sig[12],color="#7570b3",size=10)
odds_figure_vf
dev.off()




