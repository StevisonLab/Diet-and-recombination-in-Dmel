

## Crossover Interference between the y-cv Interval and the cv-v Interval

### Crossover interference calculation

#COI <- recomb
COI=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M+X0100 + X1011+X0010 + X1101+X0110 + X1001~MaternalVial+PaternalStock+Treatment+adj_age,data=recomb,FUN=sum,na.rm=T)
COI$COI=COI$co_inds.sum/COI$num_M.sum
COI$ycv_rr=COI$ycv_count.sum/COI$num_M.sum
COI$cvv_rr=COI$cvv_count.sum/COI$num_M.sum
COI$vf_rr=COI$vf_count.sum/COI$num_M.sum
COI=COI[COI$num_M.sum>=5,]

# expected vs observed
COI$Exp_DCO_ycv_cvv <-  (COI$ycv_rr * COI$cvv_rr)*COI$num_M
COI$Obs_DCO_ycv_cvv <-  COI$X0100 + COI$X1011

COI$Exp_DCO_cvv_vf <-  (COI$cvv_rr * COI$vf_rr)*COI$num_M
COI$Obs_DCO_cvv_vf <-  COI$X0010 + COI$X1101

COI$Exp_DCO_ycv_vf <-  (COI$ycv_rr * COI$vf_rr)*COI$num_M
COI$Obs_DCO_ycv_vf <-  COI$X0110 + COI$X1001


# coefficient of coincidence and interference
COI$COC_ycv_cvv <- COI$Obs_DCO_ycv_cvv / COI$Exp_DCO_ycv_cvv
COI$Interference_ycv_cvv <- 1 - COI$COC_ycv_cvv

COI$COC_cvv_vf <- COI$Obs_DCO_cvv_vf / COI$Exp_DCO_cvv_vf
COI$Interference_cvv_vf <- 1 - COI$COC_cvv_vf

COI$COC_ycv_vf <- COI$Obs_DCO_ycv_vf / COI$Exp_DCO_ycv_vf
COI$Interference_ycv_vf <- 1 - COI$COC_ycv_vf

hist(COI$Interference_ycv_cvv)
hist(COI$Interference_ycv_vf)
hist(COI$Interference_cvv_vf)

### Statistical model

#ycv_cvv
fit_ycv_cvv <- lm(Interference_ycv_cvv ~ Treatment * PaternalStock, data=COI)
#summary(fit_ycv_cvv)

anova_coi_ycv_cvv <- Anova(fit_ycv_cvv, test="Chisq")
anova_coi_ycv_cvv
#write.csv(anova_coi_ycv_cvv,"interference_ycv_cvv_model_table.csv")

#ycv_vf
fit_ycv_vf <- lm(Interference_ycv_vf ~ Treatment * PaternalStock, data=COI)
#summary(fit_ycv_vf)

anova_coi_ycv_vf <- Anova(fit_ycv_vf, test="Chisq")
anova_coi_ycv_vf
#write.csv(anova_coi_ycv_cvv,"interference_ycv_cvv_model_table.csv")

#cvv_vf
fit_cvv_vf <- lm(Interference_cvv_vf ~ Treatment * PaternalStock, data=COI)
#summary(fit_cvv_vf)

anova_coi_cvv_vf <- Anova(fit_cvv_vf, test="Chisq")
anova_coi_cvv_vf
write.csv(anova_coi_ycv_cvv,"interference_ycv_cvv_model_table.csv")

fit_contrast_cvv_vf <- emmeans::emmeans(fit_cvv_vf, "Treatment", by="PaternalStock", mode="kenward-roger")
fit_contr_cvv_vf <- contrast(fit_contrast_cvv_vf, method="trt.vs.ctrl")

pheno_contr_cvv_vf <- as.data.frame(summary(fit_contr_cvv_vf))
pheno_contr_cvv_vf
write.csv(pheno_contr_cvv_vf,"../stats/interference_cvv_vf_posthoc_table.csv")


### Figures
COI$treat=ifelse(COI$Treatment=="2x",2,ifelse(COI$Treatment=="1x",1,0.5))

COI_figure_ycv_cvv <- ggplot(aes(y=Interference_ycv_cvv, x=as.numeric(treat), col=PaternalStock), data=COI) +
  xlab("Caloric Density") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) + 
  ggtitle("Interference between the y-cv and cv-v intervals") +
  theme_base() +
  ylim(min(na.omit(COI$Interference_ycv_cvv))-0.5,
       max(na.omit(COI$Interference_ycv_cvv))+1) +
  geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun = median, geom="line", aes(group=PaternalStock),na.rm=TRUE, linewidth=1) #+
COI_figure_ycv_cvv

COI_figure_ycv_vf <- ggplot(aes(y=Interference_ycv_vf, x=as.numeric(treat), col=PaternalStock), data=COI) +
  xlab("Caloric Density") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) + 
  ggtitle("Interference between the y-cv and v-f intervals") +
  theme_base() +
  ylim(min(na.omit(COI$Interference_ycv_vf))-0.5,
       max(na.omit(COI$Interference_ycv_vf))+1) +
  geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun = median, geom="line", aes(group=PaternalStock), na.rm=TRUE,linewidth=1) #+
COI_figure_ycv_vf

COI_figure_cvv_vf <- ggplot(aes(y=Interference_cvv_vf, x=as.numeric(treat), col=PaternalStock), data=COI) +
  xlab("Caloric Density") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) + 
  ggtitle("Interference between the cv-v and v-f intervals") +
  theme_base() +
  ylim(min(na.omit(COI$Interference_cvv_vf))-0.5,
       max(na.omit(COI$Interference_cvv_vf))+1) +
  geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun = median, geom="line", aes(group=PaternalStock), na.rm=TRUE,linewidth=1) #+
# annotate(geom="text", x=1, y=1.25, label=sig_ycv_cvv[1], size=5) +
# annotate(geom="text", x=2, y=1.25, label=sig_ycv_cvv[2], size=5) +
# annotate(geom="text", x=3, y=1.25, label=sig_ycv_cvv[3], size=5) +
# annotate(geom="text", x=4, y=1.25, label=sig_ycv_cvv[4], size=5)
COI_figure_cvv_vf

#add sample sizes to figure with the code below:
#geom_text(check_overlap = F,hjust = 0, nudge_x = 0.05,angle=45,size=1)+

#dev.off()


### Save figures

ggsave("images/COI-ycv-cvv.png", plot=COI_figure_ycv_cvv, height=5)
ggsave("images/COI-cvv-vf.png", plot=COI_figure_cvv_vf, height=5)
ggsave("images/COI-ycv-vf.png", plot=COI_figure_ycv_vf, height=5)


# Plotting interference against genetic map distance
