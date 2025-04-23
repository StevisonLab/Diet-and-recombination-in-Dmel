fecund <- by_vialday[, c("Treatment", "PaternalStock", "vial_num", "vial_letter", "num_offspring", "num_moms")]

## Remove rows with less than 10 offspring
#### fecund <- fecund[fecund$num_offspring >= 10,]

## Calculate per vial per letter fecundity
fecund$fecundity <- fecund$num_offspring/fecund$num_moms

fecund_by_vial <- aggregate(fecundity ~ vial_num+vial_letter, data = fecund, sum)
fecund_by_vial <- merge(fecund_by_vial, backcross, by.x = "vial_num", by.y = "CrossID")

fecund_summary <- aggregate(fecundity ~ Treatment + PaternalStock, data = fecund_by_vial, mean)

#make day and vial number a factor
fecund_by_vial$vial_letter=as.factor(fecund_by_vial$vial_letter)
fecund_by_vial$vial_num=as.factor(fecund_by_vial$vial_num)

## Poisson
#fit_fec <- glmer(fecundity ~ (1|MaternalVial)+Treatment*PaternalStock*vial_letter, data = fecund_by_vial, family = "poisson")
#fecundity_stats <- Anova(fit_fec, test = "Chisq")

##Negative Binomial
fit_fec <- glmer.nb(fecundity ~ (1|MaternalVial)+Treatment*PaternalStock*vial_letter, data = fecund_by_vial)
fecundity_stats <- Anova(fit_fec)
fecundity_stats

write.csv(fecundity_stats, "output/fecundity-stats_anova.csv")

## Fecundity figure
fecund_figure <- ggplot(fecund_by_vial, aes(x=Treatment, y=fecundity, col=PaternalStock)) +
  xlab("Caloric Density") + ylab("# progeny per mom") + ggtitle("Fecundity") +
  ylim(0, 80) +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_boxplot() + scale_color_discrete(name = "Treatment")

fecund_figure
ggsave("images/fecundity-figure.png")

## Fecundity figure by day
fecund_figure <- ggplot(fecund_by_vial, aes(x=Treatment, y=fecundity, col=PaternalStock)) + facet_wrap(~vial_letter) +
  xlab("Caloric Density") + ylab("# progeny per mom") + ggtitle("Fecundity") +
  ylim(0, 80) +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_boxplot() + scale_color_discrete(name = "Treatment")

fecund_figure
ggsave("images/fecundity-figure-byDay.png")

mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="42"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="217"])

mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="42" & fecund_by_vial$Treatment=="0.5x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="217" & fecund_by_vial$Treatment=="0.5x"])

mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="42" & fecund_by_vial$Treatment=="1x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="217" & fecund_by_vial$Treatment=="1x"])

mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="42" & fecund_by_vial$Treatment=="2x"])
mean(fecund_by_vial$fecundity[fecund_by_vial$PaternalStock=="217" & fecund_by_vial$Treatment=="2x"])


#get number of vials with zero progeny
length(fecund_by_vial$vial_num[fecund_by_vial$fecundity==0 & fecund_by_vial$PaternalStock=="217"])

length(fecund_by_vial$vial_num[fecund_by_vial$fecundity==0 & fecund_by_vial$PaternalStock=="42"])
