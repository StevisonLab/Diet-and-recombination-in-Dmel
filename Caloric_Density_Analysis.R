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
source("scripts/haplotype-groups.R")

# Read in cleaned up anonymized phenotype data
data=read.csv(file="rawdata/phenotype-data-anon.csv",stringsAsFactors = F)

# Cleanup backcross form data (F1 cross info)
source("scripts/backcross-cleanup.R")

# Merge phenotype data with backcross/treatment data
source("scripts/merge-data.R")

#save cleaned merged datasets
write.csv(merged,file="output/phenotyping_data_cleaned.csv",row.names = F,quote = F)
write.csv(by_vialday,file="output/data_cleaned_by_vialday.csv",row.names = F,quote = F)

# Check for happlotype skews
source("scripts/haplotype-bias.R")

# Calculate mean recombination by treatment and strain
source("scripts/recombination-rate.R")

# Calculate fecundity by vial by day and avg fecundity/mom/day
source("scripts/fecundity.R")

# Combine all data into one summary data.frame
summary <- merge(recomb_summary, fecund_summary, by = c("Treatment", "PaternalStock"))

# Calculate haplotype bias
source("scripts/haplotype-bias.R")

# Calculate Crossover Interference
source("scripts/COI.R")
