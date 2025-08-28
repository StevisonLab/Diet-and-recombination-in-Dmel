##### Load required libraries
library(pbkrtest)
library(dplyr)
library(reshape2)
library(stringi)
library(stringr)
library(ggplot2)
library(emmeans)
library(lme4)
library(car)
library(ggthemes)
library(doBy)
library(MASS)
library(cowplot)
library(lmerTest)
library(optimx)
library(ggpubr)

##### RUN ANALYSIS
##
### Load custom functions
## Assign stars based on level of significance
eval.significance <- function(input) {
  ifelse(input$p.value < 0.001, "***",
         ifelse(input$p.value < 0.01, "**",
                ifelse(input$p.value < 0.05, "*", "")))
}

### Define haplotype groups
source("scripts/01.haplotype-groups.R")

# Read in cleaned up anonymized phenotype data
data=read.csv(file="rawdata/phenotype-data-anon.csv",stringsAsFactors = F)

# Cleanup backcross form data (F1 cross info)
source("scripts/02.backcross-cleanup.R")

# Merge phenotype data with backcross/treatment data
source("scripts/03.merge-data.R")

#save cleaned merged datasets
write.csv(merged,file="output/phenotyping_data_cleaned.csv",row.names = F,quote = F)
write.csv(by_vialday,file="output/data_cleaned_by_vialday.csv",row.names = F,quote = F)

# Check for happlotype skews
source("scripts/04.haplotype-bias.R")

# Calculate mean recombination by treatment and strain; runs COI script within it!
source("scripts/05.recombination-rate.R")

# Calculate fecundity by vial by day and avg fecundity/mom/day
source("scripts/06.fecundity.R")

# Combine all data into one summary data.frame
summary <- merge(recomb_summary, fecund_summary, by = c("Treatment", "PaternalStock"))

# Ovary Analysis
source("scripts/07.Ovaries.R")

# Analyze Testis Measurements
source("scripts/08.testis_analysis.R")

# Analyze Body Mass
source("scripts/09.Weights.R")

# Analyze RQ
source("scripts/10.RQ.R")

# Analyze RNA Seq
source("scripts/11.RNA_sequencing.R")
source("scripts/12.piRNA.R")
