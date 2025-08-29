##Physiology Figures
library(ggplot2)
library(readr)
library(tidyr)
library(ggsignif)
library(lme4)
library(car)
library(multcompView)
library(multcomp)
library(doBy)
library(dplyr)
library(stringr)

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

RQ_summary <- resp %>%
  group_by(Stock, Treatment) %>%
  summarise(mean_RQ = mean(RQ, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = mean_RQ) %>%
  mutate(fold_change = `0.5` / `2`,
         log2_fold_change = log2(fold_change),
         Trait = "RQ") %>%
  dplyr::select(Stock, Trait, fold_change, log2_fold_change)

##Body mass
W <- read_csv("rawdata/Weight.csv", col_types = cols(Vial = col_number(), Round = col_number(), Weight = col_number(), Stock = col_character(), Treatment = col_character(), Age = col_number()))
W$Age= as.factor(W$Age)
W$Round= as.factor(W$Round)

## Convert grams to mg in the Weight Data
W$Weights = W$Weight * 1000

mass_summary <- W %>%
  group_by(Stock, Treatment) %>%
  summarise(mean_weight = mean(Weights, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = mean_weight) %>%
  mutate(fold_change = `0.5` / `2`,
         log2_fold_change = log2(fold_change),
         Trait = "Weight") %>%
  dplyr::select(Stock, Trait, fold_change, log2_fold_change)

##Oocytes
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

##Because is not significantly different both rounds will be analyzed as one
N1Key$Stock=N1Key$`Male Stock`
N2Key$vial_num = N2Key$`Mom Vial`
N2Key$Age = N2Key$Age.x
N1 <- N1Key %>% 
  dplyr::select(1:10, 11,13,21,22)
N2 <- N2Key %>% 
  dplyr::select(1,2,4:10, 12,13, 16:18)

##Merge 2 rounds
Ovas <- bind_rows(N1, N2)

Ovas$Number_Oocyte=as.numeric(Ovas$Number_Oocyte)

fun <- function(x){
  c(m=mean(x), v=var(x), n=sum(x), min=min(x), max=max(x))}
fun1 <- function(x){
  c(m=mean(x), v=var(x), n=max(x), min=min(x), max=max(x))}

Ovas$oocyte_cor <- Ovas$Number_Oocyte/Ovas$Ovariole
Ovas$Stock = as.factor(Ovas$Stock)
Ovas$TUNEL_Ovariole=ifelse(Ovas$TUNEL_Cell>0, "Positive", "Negative") 
Ovas$TUNEL_Ovariole=as.factor(Ovas$TUNEL_Ovariole)
Ovas$Tunel_ova = (Ovas$TUNEL_Cell/Ovas$Number_Oocyte)*100
Ovas <- Ovas %>%
  filter(!is.na(Stock))

oocyte_summary <- Ovas  %>%
  group_by(Stock, Treatment) %>%
  summarise(mean_oocyte = mean(oocyte_cor, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = mean_oocyte) %>%
  mutate(fold_change = `0.5` / `2`,
         log2_fold_change = log2(fold_change),
         Trait = "oocyte_cor") %>%
  dplyr::select(Stock, Trait, fold_change, log2_fold_change)

tunel_summary <- Ovas %>%
  group_by(Stock, Treatment) %>%
  summarise(mean_tunel = mean(TUNEL_Cell, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = mean_tunel) %>%
  mutate(fold_change = `0.5` / `2`,
         log2_fold_change = log2(fold_change),
         Trait = "TUNEL_cell") %>%
  dplyr::select(Stock, Trait, fold_change, log2_fold_change)

##Testis

meta <- read.csv(file="rawdata/Testis_meta_data.csv", stringsAsFactors = F)
raw<-read.csv(file="rawdata/Testis_raw_9_13.csv")
meta$Picture_Code=paste(meta$Vial_Code,meta$Age,meta$Picture_Number,sep="_")
#con<-meta$Picture_Code
meta$Picture_Code=str_remove(meta$Picture_Code,".tif")
#merged datasets using function "merge"
Data_Merge = (merge(meta,raw,by="Picture_Code",all=F))

#make model parameters factors
Data_Merge$Treatment=as.factor(Data_Merge$Treatment)
Data_Merge$Age=as.factor(Data_Merge$Age)
Data_Merge$Strain=as.factor(Data_Merge$Strain)


testis_summary <- Data_Merge %>%
  rename(Stock = Strain) %>%  # Rename Strain to Stock
  group_by(Stock, Treatment) %>%
  summarise(mean_length = mean(Length.mm, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = mean_length) %>%
  mutate(fold_change = `0.5x` / `2x`,
         log2_fold_change = log2(fold_change),
         Trait = "Length.mm") %>%
  dplyr::select(Stock, Trait, fold_change, log2_fold_change)

##Fecundity
fec <- read.csv("output/fecund_by_vial.csv", stringsAsFactors = F)
fec <- fecund_by_vial
fecund_summary$PaternalStock=as.factor(fecund_summary$PaternalStock)
fecund_summary$Treatment=as.factor(fecund_summary$Treatment)

fec_summary <- fec %>%
  rename(Stock = PaternalStock) %>%  # Rename to Stock
  group_by(Stock, Treatment) %>%
  summarise(mean_fecundity = mean(fecundity, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = mean_fecundity) %>%
  mutate(
    fold_change = `0.5x` / `2x`,
    log2_fold_change = log2(fold_change),
    Trait = "Fecundity"
  ) %>%
  dplyr::select(Stock, Trait, fold_change, log2_fold_change)

fec_stock_avg <- fec %>%
  rename(Stock = PaternalStock) %>%
  filter(!is.na(Stock)) %>%
  group_by(Stock) %>%
  summarise(mean_fecundity = mean(fecundity, na.rm = TRUE), .groups = "drop")

##all data combined
final_results <- bind_rows(
  RQ_summary,
  mass_summary,
  fec_summary,
  testis_summary,
  oocyte_summary
  #tunel_summary
) %>%
  arrange(Stock, Trait)

final_results <- final_results %>%
  mutate(
    fill_color = case_when(
      Stock == 42 & log2_fold_change < 0 ~ "#a63603",  # dark for negative
      Stock == 42 & log2_fold_change >= 0 ~ "#fdd0a2", # light for positive
      Stock == 217 & log2_fold_change < 0 ~ "#08519c",
      Stock == 217 & log2_fold_change >= 0 ~ "#c6dbef",
      TRUE ~ "gray80"
    )
  )
final_results$Trait=as.character(final_results$Trait)
final_results <- final_results %>%
  mutate(
    Trait = as.character(Trait),  # ensure it's character
    Trait = dplyr::recode(
      Trait,
      "Weight" = "Body mass",
      #"TUNEL_Cell" = "DNA degradation",
      "oocyte_cor" = "Number of oocytes",
      "Length.mm" = "Testes Length"
    )
  )
final_results$Trait <- factor(final_results$Trait, levels = c("Body mass", "RQ", "Testes Length", "Number of oocytes", "Fecundity"))
plot_42 <- final_results %>%
  filter(Stock == 42) %>%
  ggplot(aes(x = log2_fold_change, y = Trait, fill = fill_color)) +
  geom_col(width = 0.6) +
  xlim(-1.07,0.25)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_identity() +
  labs(
    title = "Stock 42: log2 Fold Change by Trait",
    x = "log2 Fold Change ",
    y = "Trait",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "black") 
  )

# Plot for Stock 217
plot_217 <- final_results %>%
  filter(Stock == 217) %>%
  ggplot(aes(x = log2_fold_change, y = Trait, fill = fill_color)) +
  geom_col(width = 0.6) +
  xlim(-1.07,0.25)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_identity() +
  labs(
    title = "Stock 217: log2 Fold Change by Trait",
    x = "log2 Fold Change",
    y = "Trait",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")  
  )

plot_42
plot_217

##Combined without the DNA degradation; Figure 4A
#final_results %>%
#  filter(Stock %in% c(42, 217)) %>%
Figure4A=ggplot(final_results,aes(x = log2_fold_change, 
             y = Trait, 
             fill = fill_color, 
             group = Stock)) +   # <-- ensures side-by-side by stock
  geom_col(width = 0.6, position = position_dodge(width = 0.7)) +
  xlim(-1.07, 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_identity() +  # <-- uses your existing colors
  labs(
    #title = "log2 Fold Change by Trait (Stock 42 vs Stock 217)",
    x = expression("log"[2]*"FC"),
    y = "Trait",
    fill = "Treatment"
  ) +
  theme_base()  #+

#  theme(
#    plot.title = element_text(face = "bold"),
#    axis.title.x = element_text(face = "bold"),
#    axis.title.y = element_text(face = "bold"),
#    axis.text.y = element_text(face = "bold", size = 14),
#    axis.text.x = element_text(face = "bold", size = 14),
#    legend.title = element_text(face = "bold"),
#    axis.text = element_text(color = "black")
#  )

Figure4A
ggsave("images/Figure4A.png",plot=Figure4A, height=5)

