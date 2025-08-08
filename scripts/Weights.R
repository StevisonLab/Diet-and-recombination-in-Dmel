##Weights for nutrition project Spring 2023 Round 2 and 3
library(ggplot2)
library(readr)
library(tidyr)
library(ggsignif)
library(lme4)
library(car)
library(multicompView)
library(multicomp)
library(doBy)
##Read Data
W <- read_csv("rawdata/Weight.csv", col_types = cols(Vial = col_number(), Round = col_number(), Weight = col_number(), Stock = col_character(), Treatment = col_character(), Age = col_number()))
W$Age= as.factor(W$Age)
W$Round= as.factor(W$Round)

## Convert grams to mg in the Weight Data
W$Weights = W$Weight * 1000

## Subset by day and Stock
Day0_10 <- subset(W, W$Age=="0"|W$Age=="10", na.rm=TRUE, select = c(Treatment, Sex, Round, Stock, Weights, Vial, Age))
Day0A <- subset(W, W$Age=="0", na.rm=TRUE, select = c(Treatment, Sex, Round, Stock, Weights, Vial))
AllF <- subset(W, W$Sex=="F", na.rm=TRUE, select = c(Treatment, Sex, Round, Stock, Weights, Vial, Age))
Day0 <- subset(W, W$Age=="0" & W$Sex=="F", na.rm=TRUE, select = c(Treatment, Age, Sex, Round, Stock, Weights, Vial))
Day10 <- subset(W, W$Age=="10" & W$Sex=="F", na.rm=TRUE, select = c(Treatment, Age, Sex, Round, Stock, Weights, Vial))
S42 <- subset(W, W$Sex=="F" & W$Stock=="42" & W$Age=="0", na.rm=TRUE, select = c(Treatment,Sex, Round, Stock, Weights, Vial))
S217 <- subset(W, W$Sex=="F" & W$Stock=="217" & W$Age=="0", na.rm=TRUE, select = c(Treatment,Sex, Round, Stock, Weights, Vial))
A42 <- subset(W, W$Stock=="42", na.rm=TRUE, select = c(Treatment, Age, Sex, Round, Stock, Weights, Vial))
A217 <- subset(W, W$Stock=="217", na.rm=TRUE, select = c(Treatment, Age, Sex, Round, Stock, Weights, Vial))


Day10$Stock= as.factor(Day10$Stock)
Day10$Treatment= as.factor(Day10$Treatment)
A217$Round= as.numeric(A217$Round)
A42$Round= as.numeric(A42$Round)

##Summary of the data
fun <- function(x){
  c(m=mean(x), v=var(x), n=length(x), min=min(x), max=max(x))
}
SummaryW <- summary_by(Weights ~ c(Treatment, Sex, Stock), data=W, FUN=fun)
write_csv(SummaryW, "output/weight_Summary1.csv")


##Model for both sex day 0 and 10
res.aov <- lmer(Weights ~ (1|Vial) + Treatment * Age * Stock * Sex, data = Day0_10)
summary(res.aov)$coefficients
W0_10A <- Anova(res.aov)
W0_10A
W0_10em <- emmeans(res.aov, ~Sex*Treatment*Age*Stock)
plot(W0_10em)
cldW0_10em <-cld(W0_10em, Letters = letters, alpha = 0.05)
W0_10S=summaryBy(Weights~Treatment+Sex+Age+Stock,data=Day0_10,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
W0_10S.summarized=merge(W0_10S,cldW0_10em, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
unique(W0_10S.summarized$.group)
W0_10S.summarized$group=ifelse(W0_10S.summarized$.group=="  bc ","bc", ifelse(W0_10S.summarized$.group=="    d","d", ifelse(W0_10S.summarized$.group=="  bcd","bcd", ifelse(W0_10S.summarized$.group==" a   ","a", ifelse(W0_10S.summarized$.group=="   cd","cd", ifelse(W0_10S.summarized$.group=="  b  ","b", W0_10S.summarized$.group))))))

res.aov2 <- lmer(Weights ~ (1|Vial) + Treatment * Age * Stock, data = AllF)
summary(res.aov2)$coefficients
WFA <- Anova(res.aov2)
WFA
WFem <- emmeans(res.aov2, ~Treatment*Age*Stock)
WFempaired <- pairs(WFem, adjust = "tukey")
write_csv(as.data.frame(WFempaired), "./output/WFemm.csv")
cldWFem <-cld(WFem, Letters = letters, alpha = 0.05)
WFS=summaryBy(Weights~Treatment+Age+Stock,data=AllF,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
WFS.summarized=merge(WFS,cldWFem, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
unique(WFS.summarized$.group)
WFS.summarized$group=ifelse(WFS.summarized$.group==" abcdefgh ","abcdefgh", ifelse(WFS.summarized$.group=="         i","i", ifelse(WFS.summarized$.group==" abcdefghi","abcdefghi", ifelse(WFS.summarized$.group==" a c  f   ","acf", ifelse(WFS.summarized$.group==" abcd fg  ","abcdfg", ifelse(WFS.summarized$.group=="  b de ghi","bdeghi", ifelse(WFS.summarized$.group=="     e  hi","ehi", ifelse(WFS.summarized$.group=="      fghi","fghi", ifelse(WFS.summarized$.group==" ab       ","ab", ifelse(WFS.summarized$.group=="   cdefghi","cdefghi", ifelse(WFS.summarized$.group==" abcde    ","abcde", WFS.summarized$.group)))))))))))


#Function to add N to the boxplot of Weight
stat_box_data3 <- function(y, upper_limit = 2) {
  return( 
    data.frame(
      y = 0.5,
      label = paste("n=",length(y))
    )
  )
}


##Function for SD
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
###---------------------------------##

##All data
SDA0_10 <- data_summary(Day0_10, varname="Weights", 
                      groupnames=c("Age", "Treatment", "Sex", "Stock"))
write_csv(SDA0_10, "output/SummaryWeight0_10.csv")

SDA42 <- data_summary(A42, varname="Weights", 
                     groupnames=c("Age", "Treatment", "Sex"))

SDA217 <- data_summary(A217, varname="Weights", 
                      groupnames=c("Age", "Treatment", "Sex"))

ALLFm <- data_summary(AllF, varname="Weights", 
                      groupnames=c("Age", "Treatment", "Stock"))
write_csv(ALLFm, "output/SummaryFemalesDaily.csv")
pdf("images/Bodymass_Spring23.pdf")
#Plot for day 0 and 10 for both sexes, by treatment and stock
WeightA0_10=ggplot(SDA0_10, aes(x=Age, y=Weights, group=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=Weights-sd, ymax=Weights+sd), width=.1) +
  scale_color_manual(values=c("yellow","green4", "blue")) +
  ggtitle("Body mass (mg) for Day 0 and 10 by treatments and sexes") + theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Body mass (mg)") +
  theme(plot.title = element_text(size=10)) +
  geom_text(data = W0_10S.summarized, aes(y = 1.5, x = Age, label = group), position = position_dodge(width = 0.75)) +
  geom_line() + geom_point() + facet_wrap(vars(Sex,))
WeightA0_10
Weight217_10=ggplot(SDA217, aes(x=Age, y=Weights, group=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=Weights-sd, ymax=Weights+sd), width=.1) +
  scale_color_manual(breaks = c("0.5", "1", "2"), values=c("#FEEDDE","#FD8D3C", "#A63603")) +
  ggtitle("Body mass (mg) for Day 0 and 10 by treatments and sexes") + theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Body mass (mg)") +
  theme(plot.title = element_text(size=10)) +
  geom_line() + geom_point() + facet_wrap(~Sex)
Weight217_10

Weight42_10=ggplot(SDA42, aes(x=Age, y=Weights, group=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=Weights-sd, ymax=Weights+sd), width=.1) +
  scale_color_manual(breaks = c("0.5", "1", "2"), values=c("#EFF3FF","#6BAED6", "#08519C")) +
  ggtitle("Body mass (mg) for Day 0 and 10 by treatments and sexes") + theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Body mass (mg)") +
  theme(plot.title = element_text(size=10)) +
  geom_line() + geom_point() + facet_wrap(~Sex)
Weight42_10






##Boxplots
WeightA0_10=ggplot(Day0_10, aes(x=Age, y=Weights, fill=Treatment)) +
  ylab("Body mass (mg)") +
  theme_bw()+
  theme(strip.background=element_rect(fill="white"))+
  theme(strip.text=element_text(color="black", face="bold")) +
  stat_summary(fun.data = stat_box_data3, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_text(data = W0_10S.summarized, aes(y = 1.75, x = Age, label = group), position = position_dodge(width = 0.75)) +
  ggtitle("Body mass (mg) for Day 0 and 10 by treatments and sexes") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=12, face = "bold")) + 
  geom_boxplot() + facet_grid(vars(Stock), vars(Sex))
WeightA0_10

##Plot for all ages for females by treatment and stock
annotation_df <- data.frame(
  Stock = c("42", "217", "42", "217", "217", "217", "217"),
  Treatment = c("2", "1", "2", "1", "0.5", "1", "1"),
  start = c("2", "0", "2", "6", "0", "0", "0"),
  end = c("8", "6", "10", "10", "10", "10", "4"),
  y = c(1.4, 1.4, 1.45, 1.38, 1.5, 1.47, 1.43),
  label = c("***", " ** ", " ** ", " * ", " * ", " * ", " * ")
)
annotation_df$group <- 1:nrow(annotation_df)


WeightF=ggplot(ALLFm, aes(x=Age, y=Weights, group=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=Weights-sd, ymax=Weights+sd), width=.1) +
  ggtitle("Body mass females by stock, day, and treatment") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=10)) + 
  geom_signif(
    data = annotation_df,
    aes(xmin = annotation_df$start, xmax = annotation_df$end, annotations = annotation_df$label, y_position = annotation_df$y, group = group),
    textsize = 3,
    step_increase = 0.05,
    tip_length = 0.01,
    manual = TRUE
  ) +
  #geom_signif(comparisons = list(c("2", "8"), 
              #                   c("0", "6"), c("2", "10"), c("6", "10"), c("0", "10"), c("0", "4")), map_signif_level = TRUE, margin_top = 0.08,
              #step_increase = 0.05,
              #tip_length = 0.01) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_line(alpha=.6) + geom_point(color="orange",alpha=.6) + facet_wrap(~Stock)
WeightF
AllF$Stock=as.factor(AllF$Stock)
AllF$Treatment=as.factor(AllF$Treatment)

ALLFm <- data_summary(AllFm, varname="Weights", 
                      groupnames=c("Age", "Treatment", "Stock"))
WeightF=ggplot(ALLFm, aes(x=Age, y=Weights, fill = factor(Stock),alpha = Treatment, color=Stock, group=Treatment ), 
    stat = "identity") + 
  geom_signif(
    data = annotation_df,
    aes(xmin = annotation_df$start, xmax = annotation_df$end, annotations = annotation_df$label, y_position = annotation_df$y, group = group),
    textsize = 3,
    step_increase = 0.05,
    tip_length = 0.01,
    manual = TRUE
  ) +
  scale_alpha_discrete(range = c(0.2,1)) +
  scale_color_manual(values = c("orange3","blue3"))+
  geom_line() + geom_point() + facet_wrap(~Stock)


WeightF
data_subset <- Day0_10[Day0_10$Sex %in% "F", ]
data_model <- W0_10S.summarized[W0_10S.summarized$Sex %in% "F", ]
WeightAF=ggplot(data_subset, aes(x=Age, y=Weights, color = Stock,alpha = Treatment, fill=Stock), 
               stat = "identity") + 
  ylab("Body mass (mg)") +
  scale_alpha_discrete(range = c(0.2,1)) +
  scale_color_manual(values = c("orange4","blue4"))+
  scale_fill_manual(values = c("orange3","blue3"))+
  stat_summary(fun.data = stat_box_data3, geom = "text", fun = median,
               position = position_dodge(width = 0.75)) +
  geom_text_repel(data = data_model, aes(y = 1.8, x = Age, label = group), vjust=0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot() + facet_wrap(~Stock)
WeightAF
data_subset1 <- Day0_10[Day0_10$Sex %in% "M", ]
data_model1 <- W0_10S.summarized[W0_10S.summarized$Sex %in% "M", ]
WeightAF=ggplot(data_subset1, aes(x=Age, y=Weights, color = Stock,alpha = Treatment, fill=Stock), 
                stat = "identity") + 
  ylab("Body mass (mg)") +
  scale_alpha_discrete(range = c(0.2,1)) +
  scale_color_manual(values = c("orange4","blue4"))+
  scale_fill_manual(values = c("orange3","blue3"))+
  geom_text_repel(data = data_model1, aes(y = 1, x = Age, label = group), vjust=0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot() + facet_wrap(~Stock)
WeightAF
dev.off()
