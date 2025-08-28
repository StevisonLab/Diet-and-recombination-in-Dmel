

## Crossover Interference between the y-cv Interval and the cv-v Interval

### Crossover interference calculation

#COI <- recomb
COI=summaryBy(co_inds+ycv_count+cvv_count+vf_count+num_M+X0100 + X1011+X0010 + X1101+X0110 + X1001~PaternalStock+Treatment+vial_letter,data=recomb,FUN=sum,na.rm=T)
COI$COI=COI$co_inds.sum/COI$num_M.sum
COI$ycv_rr=COI$ycv_count.sum/COI$num_M.sum
COI$cvv_rr=COI$cvv_count.sum/COI$num_M.sum
COI$vf_rr=COI$vf_count.sum/COI$num_M.sum
COI=COI[COI$num_M.sum>=5,]

# expected vs observed
COI$Exp_DCO_ycv_cvv <-  round((COI$ycv_rr * COI$cvv_rr)*COI$num_M,0)
COI$Obs_DCO_ycv_cvv <-  COI$X0100 + COI$X1011

COI$Exp_DCO_cvv_vf <-  round((COI$cvv_rr * COI$vf_rr)*COI$num_M,0)
COI$Obs_DCO_cvv_vf <-  COI$X0010 + COI$X1101

COI$Exp_DCO_ycv_vf <-  round((COI$ycv_rr * COI$vf_rr)*COI$num_M,0)
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


### Statistical models
COI$PaternalStock=as.factor(COI$PaternalStock)
#ycv_cvv
#remove expected zero
COI_test <- COI %>%
  filter(Exp_DCO_ycv_cvv > 0 )
fit_ycv_cvv <- glm(cbind(Obs_DCO_ycv_cvv) ~ Treatment * PaternalStock + offset(log(Exp_DCO_ycv_cvv)), data=COI_test,family=quasipoisson(link = "log"))
#summary(fit_ycv_cvv)

anova_coi_ycv_cvv <- anova(fit_ycv_cvv, test="F")
anova_coi_ycv_cvv
write.csv(anova_coi_ycv_cvv,"output/interference_ycv_cvv_model_table.csv")

#fit_contrast_ycv_cvv <- emmeans::emmeans(fit_ycv_cvv, "Treatment", by="PaternalStock", mode="kenward-roger")
#fit_contr_ycv_cvv <- contrast(fit_contrast_ycv_cvv, method="trt.vs.ctrl")

#pheno_contr_ycv_cvv <- as.data.frame(summary(fit_contr_ycv_cvv))
#pheno_contr_ycv_cvv
#write.csv(pheno_contr_ycv_cvv,"output/interference_ycv_cvv_posthoc_table.csv")

#ycv_vf
#remove expected zero
COI_test <- COI %>%
  filter(Exp_DCO_ycv_vf > 0 )
fit_ycv_vf <- glm(Obs_DCO_ycv_vf ~ Treatment * PaternalStock + offset(log(Exp_DCO_ycv_vf)), data=COI_test,family=quasipoisson(link = "log"))
#summary(fit_ycv_vf)

anova_coi_ycv_vf <- anova(fit_ycv_vf, test="F")
anova_coi_ycv_vf
write.csv(anova_coi_ycv_cvv,"output/interference_ycv_cvv_model_table.csv")

#fit_contrast_ycv_vf <- emmeans::emmeans(fit_ycv_vf, "Treatment", by="PaternalStock", mode="kenward-roger")
#fit_contr_ycv_vf <- contrast(fit_contrast_ycv_vf, method="trt.vs.ctrl")

#pheno_contr_ycv_vf <- as.data.frame(summary(fit_contr_ycv_vf))
#pheno_contr_ycv_vf
#write.csv(pheno_contr_ycv_vf,"output/interference_ycv_vf_posthoc_table.csv")

#cvv_vf
#remove expected zero
COI_test <- COI %>%
  filter(Exp_DCO_cvv_vf > 0 )
fit_cvv_vf <- glm(Obs_DCO_cvv_vf ~ Treatment * PaternalStock + offset(log(Exp_DCO_cvv_vf)), data=COI_test,family=quasipoisson(link = "log"))
#summary(fit_cvv_vf)

anova_coi_cvv_vf <- anova(fit_cvv_vf, test="F")
anova_coi_cvv_vf
write.csv(anova_coi_ycv_cvv,"output/interference_ycv_cvv_model_table.csv")

fit_contrast_cvv_vf <- emmeans::emmeans(fit_cvv_vf, "Treatment", by="PaternalStock", mode="kenward-roger")
fit_contr_cvv_vf <- contrast(fit_contrast_cvv_vf, method="trt.vs.ctrl")

pheno_contr_cvv_vf <- as.data.frame(summary(fit_contr_cvv_vf))
pheno_contr_cvv_vf
write.csv(pheno_contr_cvv_vf,"output/interference_cvv_vf_posthoc_table.csv")


### Figures
COI$treat=ifelse(COI$Treatment=="2x",2,ifelse(COI$Treatment=="1x",1,0.5))
col_pal2=c("#fdd0a2","#fd8d3c","#a63603","#c6dbef","#6baed6","#08519c")
           
COI_figure_ycv_cvv <- ggplot(aes(y=Interference_ycv_cvv, x=as.numeric(treat), col=PaternalStock), data=COI) + 
  xlab("Caloric Density") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) + 
  ggtitle("Interference of y-cv-v interval") + 
  theme_base() + scale_color_manual(values = c("#fd8d3c","#6baed6"))+
#  ylim(0.6,1.5) + 
#  geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
  stat_summary(fun = median, geom="line", aes(group=PaternalStock),na.rm=TRUE, linewidth=1) +
  theme(plot.title = element_text(size = 8))

COI_figure_ycv_cvv

COI_figure_ycv_vf <- ggplot(aes(y=Interference_ycv_vf, x=as.numeric(treat), col=PaternalStock), data=COI) + 
  xlab("Caloric Density") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) + 
  ggtitle("Interference of y-cv-v-f interval") +
  theme_base() + scale_color_manual(values = c("#fd8d3c","#6baed6"))+
#  ylim(-0.4,1.1) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
#  geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun = median, geom="line", aes(group=PaternalStock), na.rm=TRUE,linewidth=1)+
  theme(plot.title = element_text(size = 8))

COI_figure_ycv_vf

COI_figure_cvv_vf <- ggplot(aes(y=Interference_cvv_vf, x=as.numeric(treat), col=PaternalStock), data=COI) + 
  xlab("Caloric Density") + ylab("Crossover Interference") +
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) + 
  ggtitle("Interference of cv-v-f interval") +
  theme_base() + scale_color_manual(values = c("#fd8d3c","#6baed6"))+
#  ylim(-0.4,1.1) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
#  geom_point(size=1.5,na.rm=TRUE) +
  stat_summary(fun = median, geom="line", aes(group=PaternalStock), na.rm=TRUE,linewidth=1) +
# annotate(geom="text", x=2, y=0.6, label="**", size=7,col="red")+
  theme(plot.title = element_text(size = 8)) 

COI_figure_cvv_vf



### Save figures

ggsave("images/COI-ycv-cvv.png", plot=COI_figure_ycv_cvv, height=5)
ggsave("images/COI-cvv-vf.png", plot=COI_figure_cvv_vf, height=5)
ggsave("images/COI-ycv-vf.png", plot=COI_figure_ycv_vf, height=5)


# Plotting interference against genetic map distance

# Need to reorganize the data frame: Interference measure and Genetic Distance Measure
wide_df=COI[,c("PaternalStock","Treatment","vial_letter","Interference_cvv_vf","Interference_ycv_cvv","Interference_ycv_vf")]
colnames(wide_df)=c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df=melt(wide_df,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Interference")

wide_df2=COI[,c("PaternalStock","Treatment","vial_letter","ycv_count.sum","cvv_count.sum","vf_count.sum","co_inds.sum")]
wide_df2$ycv_cvv=rowSums(wide_df2[,c("ycv_count.sum","cvv_count.sum")],na.rm=TRUE)
wide_df2$cvv_vf=rowSums(wide_df2[,c("cvv_count.sum","vf_count.sum")],na.rm=TRUE)
wide_df2$ycv_vf=wide_df2$co_inds.sum
wide_df2=wide_df2[,c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")]
long_df2=melt(wide_df2,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Num_COs")

wide_df3=COI[,c("PaternalStock","Treatment","vial_letter","num_M.sum")]
wide_df3$sub1=wide_df3$num_M.sum
wide_df3$sub2=wide_df3$num_M.sum
long_df3=melt(wide_df3,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Total")

wide_df5=COI[,c("PaternalStock","Treatment","vial_letter","Exp_DCO_cvv_vf","Exp_DCO_ycv_cvv","Exp_DCO_ycv_vf")]
colnames(wide_df5)=c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df5=melt(wide_df5,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Exp_DCO")

wide_df6=COI[,c("PaternalStock","Treatment","vial_letter","Obs_DCO_cvv_vf","Obs_DCO_ycv_cvv","Obs_DCO_ycv_vf")]
colnames(wide_df6)=c("PaternalStock","Treatment","vial_letter","cvv_vf","ycv_cvv","ycv_vf")
long_df6=melt(wide_df6,id.vars=c("PaternalStock","Treatment","vial_letter"),variable.name = "Interval",value.name = "Obs_DCO")

long_df4=cbind(long_df,long_df2[,5],long_df3[,5],long_df5[,5],long_df6[,5])
colnames(long_df4)=c("PaternalStock","Treatment","vial_letter","Interval","Interference","Num_COs","Total","Exp_DCO","Obs_DCO")
long_df4$Gendist=long_df4$Num_COs/long_df4$Total
long_df4$Gendist=ifelse(long_df4$Num_COs/long_df4$Total<0.5,long_df4$Num_COs/long_df4$Total,0.499+(((long_df4$Num_COs/long_df4$Total)-0.5)/100))
long_df4$kosambi_distances <- (100/4)*(log((1+(2*long_df4$Gendist))/(1-(2*long_df4$Gendist))))
long_df4$PaternalStock=as.factor(long_df4$PaternalStock)

ggplot(data=long_df4,aes(x=kosambi_distances,y=Interference,col=Treatment))+geom_line()+facet_wrap(~PaternalStock)
ggsave("images/COI-vs-RR_byDay.png", height=5)


#COI_final=summaryBy(Num_COs+ Total+ Exp_DCO +Obs_DCO~PaternalStock+Treatment+Interval,data=long_df4,FUN=sum,na.rm=T)

#COI_final$COC <- COI_final$Obs_DCO.sum / COI_final$Exp_DCO.sum
#COI_final$Interference <- 1 - COI_final$COC
#COI_final$Gendist=ifelse(COI_final$Num_COs/COI_final$Total<0.5,COI_final$Num_COs/COI_final$Total,0.499)
#COI_final$kosambi_distances <- (100/4)*(log((1+(2*COI_final$Gendist))/(1-(2*COI_final$Gendist))))

COI_final=summaryBy(Interference+kosambi_distances~PaternalStock+Treatment+Interval,data=long_df4,FUN=mean,na.rm=T)
col_pal2=c("#fdd0a2","#fd8d3c","#a63603","#c6dbef","#6baed6","#08519c")

diet_gxe_int=ggplot(data=COI_final,aes(x=kosambi_distances.mean,y=Interference.mean,col=Treatment))+geom_line()+facet_wrap(~PaternalStock) + xlab("Recombination Rate (cM)") + ylab("Crossover Interference") +   #ggtitle("Interference versus Recombination Rate") +
  theme_base() + scale_color_manual(values = col_pal2[1:3])
  #ylim(0.1,1) + xlim(0.15,0.52)
diet_gxe_int
ggsave("images/COI-vs-RR_42.png", height=5)

diet_gxe_int=ggplot(data=COI_final,aes(x=kosambi_distances.mean,y=Interference.mean,col=Treatment))+geom_line()+facet_wrap(~PaternalStock) + xlab("Recombination Rate (cM)") + ylab("Crossover Interference") +   #ggtitle("Interference versus Recombination Rate") +
  theme_base() + scale_color_manual(values = col_pal2[4:6])
#ylim(0.1,1) + xlim(0.15,0.52)
diet_gxe_int
ggsave("images/COI-vs-RR_217.png", height=5)

#pdf("images/diet_gxe_interference.pdf")
#diet_gxe_int
#dev.off()

#Make recombination Plot to match

#extract total map length 
ycv_vf=subset(long_df4,long_df4$Interval=="ycv_vf")
ycv_vf$treat=ifelse(ycv_vf$Treatment=="2x",2,ifelse(ycv_vf$Treatment=="1x",1,0.5))

#plot
#annotation_df=compare_means(kosambi_distances~treat, data=ycv_vf)
diet_gxe_map <- ggplot(ycv_vf, aes(x=treat, y=kosambi_distances,col=factor(PaternalStock)))  + 
  xlab("Caloric Density") + ylab("Recombination (cM)") + 
  scale_x_continuous(limits=c(0.4,2.1),breaks=c(0.5,1,2),labels=c("0.5x","1x","2x")) +  
  theme_base()+ scale_color_manual(values = c("#fd8d3c","#6baed6"))+
  theme(legend.title=element_blank())+ theme(plot.title = element_text(size = 8)) +
  stat_summary(fun = median, geom="line", aes(group=PaternalStock),na.rm=TRUE, linewidth=1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
  ggtitle("Recombination Rate of y-f interval") 

ggsave("images/RR-total.png", plot=diet_gxe_map, height=5)

pdf("images/diet_gxe_full.pdf")
diet_gxe_map
dev.off()

ggarrange(diet_gxe_map,COI_figure_cvv_vf,ncol=2,nrow=1,legend="bottom",common.legend = TRUE,labels=c("A","B"),font.label = list(size = 16))
ggsave("images/Figure2.png",width=8,height=4)


#Get values for Table in Paper:
COI_summary=summaryBy(Interference~PaternalStock+Treatment+Interval,data=long_df4,FUN=mean,na.rm=T)
COI_summary=summaryBy(Interference~PaternalStock+Treatment+Interval,data=long_df4,FUN=sd,na.rm=T)

