

# piRNA View

# Read in Novomagic output
raw_data=read.csv("rawdata/piRNA_get_annot_NovoMagic.csv",header=TRUE,na.strings="--")
raw_data=unique(raw_data)

#get row means
raw_data$S42_0_LC <- rowMeans(
  raw_data[, c("s42_05_0_1_fpkm", "s42_05_0_2_fpkm", "s42_05_0_3_fpkm", "s42_05_0_4_fpkm")],
  na.rm = TRUE
)
raw_data$S42_0_C=rowMeans(
  raw_data[, c("s42_1_0_1_fpkm","s42_1_0_2_fpkm","s42_1_0_3_fpkm","s42_1_0_4_fpkm")],na.rm=TRUE)
raw_data$S42_0_HC=rowMeans(
  raw_data[, c("s42_2_0_1_fpkm","s42_2_0_2_fpkm","s42_2_0_3_fpkm","s42_2_0_4_fpkm")],na.rm=TRUE)
raw_data$S217_0_LC=rowMeans(
  raw_data[, c("s217_05_0_1_fpkm","s217_05_0_2_fpkm","s217_05_0_3_fpkm","s217_05_0_4_fpkm")],na.rm=TRUE)
raw_data$S217_0_C=rowMeans(
  raw_data[, c("s217_1_0_1_fpkm","s217_1_0_2_fpkm","s217_1_0_3_fpkm","s217_1_0_4_fpkm")],na.rm=TRUE)
raw_data$S217_0_HC=rowMeans(
  raw_data[, c("s217_2_0_1_fpkm","s217_2_0_2_fpkm","s217_2_0_3_fpkm","s217_2_0_4_fpkm")],na.rm=TRUE)

raw_data=raw_data[,c(1:10,88:93)]




#reshape to long format
library(tidyr)
library(dplyr)

# Pivot longer: gather all FPKM columns into key-value pairs
long_df <- raw_data %>%
  pivot_longer(
    cols = matches("^s\\d+_"),  # selects all columns starting with 's' and digits
    names_to = "sample",
    values_to = "FPKM"
  ) %>%
  # Separate 'sample' column into strain, diet, rep, etc.
  separate(
    col = sample,
    into = c("strain", "day","diet"),
    sep = "_",
    remove = FALSE
  ) %>%
  dplyr::select(ID, strain, diet, day,FPKM)


#plots

long_df$day=as.factor(long_df$day)
long_df$treat=ifelse(long_df$diet=="C",1.0,ifelse(long_df$diet=="HC",2.0,0.5))
long_df$treat=as.factor(long_df$treat)
long_df=subset(long_df,long_df$diet=="HC" | long_df$diet=="LC")
# Plots:
library(ggplot2)
library(ggthemes)

pale=c("#fdd0a2","#c6dbef")

ggplot(data=long_df, aes(x = treat, y = log2(FPKM + 1),fill=strain)) +
  geom_boxplot() + #facet_wrap(~ID) +
  labs(title = "Expression by Diet and Strain", x = "Diet", y = "FPKM") +
  theme_minimal() +scale_fill_manual(values = pale) 
ggsave("images/piRNA.png")

#Conclusion: S42 has lower overall piRNA expression regardless of diet

ggplot(data=long_df, aes(x = as.numeric(treat), y = log2(FPKM + 1),col=strain)) +
  geom_point() + geom_line() + facet_wrap(~ID) +
  labs(title = "Expression by Diet and Strain", x = "Diet", y = "FPKM") +
  theme_minimal() + scale_color_manual(values = pale) 
ggsave("images/piRNA_byGene.png")
