# List haplotypes as combinations of 0s and 1s

#haplotypes <- levels(data$haplotype)
haplotypes <- c("X0000", "X0001", "X0010", "X0011", "X0100", "X0101", "X0110", "X0111",
                "X1000", "X1001", "X1010", "X1011", "X1100", "X1101", "X1110", "X1111")

cogs_haps <- data.frame(haplotype = haplotypes)
cogs_haps <- within(cogs_haps, {
  crossover_type = NA
  crossover_type[haplotype == "X0000"] = "NCO1"
  crossover_type[haplotype == "X1111"] = "NCO2"
  
  crossover_type[haplotype == "X1000"] = "SCO11"
  crossover_type[haplotype == "X0111"] = "SCO12"
  
  crossover_type[haplotype == "X0011"] = "SCO21"
  crossover_type[haplotype == "X1100"] = "SCO22"
  
  crossover_type[haplotype == "X0001"] = "SCO31"
  crossover_type[haplotype == "X1110"] = "SCO32"
  
  crossover_type[haplotype == "X0100"] = "DCO11"
  crossover_type[haplotype == "X1011"] = "DCO12"
  
  crossover_type[haplotype == "X0110"] = "DCO21"
  crossover_type[haplotype == "X1001"] = "DCO22"
  
  crossover_type[haplotype == "X0010"] = "DCO31"
  crossover_type[haplotype == "X1101"] = "DCO32"
  
  crossover_type[haplotype == "X0101"] = "TCO1"
  crossover_type[haplotype == "X1010"] = "TCO2"
})

co_groups <-c("NCO1", "NCO2", "SCO11", "SCO12",
              "SCO21", "SCO22", "SCO31", "SCO32",
              "DCO11", "DCO12", "DCO21", "DCO22",
              "DCO31", "DCO32", "TCO1", "TCO2")

cog_pairs <-c("NCO", "SCO1", "SCO2", "SCO3", "DCO1", "DCO2", "DCO3", "TCO")

haps_nco <- cogs_haps$haplotype[grep("^NCO", cogs_haps$crossover_type)]
haps_sco <- cogs_haps$haplotype[grep("^SCO", cogs_haps$crossover_type)]
haps_dco <- cogs_haps$haplotype[grep("^DCO", cogs_haps$crossover_type)]
haps_tco <- cogs_haps$haplotype[grep("^TCO", cogs_haps$crossover_type)]

haps_sco1 <- cogs_haps$haplotype[grep("^SCO1", cogs_haps$crossover_type)]
haps_sco2 <- cogs_haps$haplotype[grep("^SCO2", cogs_haps$crossover_type)]
haps_sco3 <- cogs_haps$haplotype[grep("^SCO3", cogs_haps$crossover_type)]
haps_dco1 <- cogs_haps$haplotype[grep("^DCO1", cogs_haps$crossover_type)]
haps_dco2 <- cogs_haps$haplotype[grep("^DCO2", cogs_haps$crossover_type)]
haps_dco3 <- cogs_haps$haplotype[grep("^DCO3", cogs_haps$crossover_type)]

cogs_nco <- co_groups[grep("^NCO", co_groups)]
cogs_sco <- co_groups[grep("^SCO", co_groups)]
cogs_dco <- co_groups[grep("^DCO", co_groups)]
cogs_tco <- co_groups[grep("^TCO", co_groups)]

cog_sco1 <- co_groups[grep("^SCO1", co_groups)]
cog_sco2 <- co_groups[grep("^SCO2", co_groups)]
cog_sco3 <- co_groups[grep("^SCO3", co_groups)]
cog_dco1 <- co_groups[grep("^DCO1", co_groups)]
cog_dco2 <- co_groups[grep("^DCO2", co_groups)]
cog_dco3 <- co_groups[grep("^DCO3", co_groups)]

# crossover groups with a crossover in the each interval (y-cv, cv-v, and v-f)
ycv_cog <- c("SCO11", "SCO12",
             "DCO11", "DCO12",
             "DCO21", "DCO22",
             "TCO1", "TCO2")

cvv_cog <- c("SCO21", "SCO22",
             "DCO11", "DCO13",
             "DCO31", "DCO32",
             "TCO1", "TCO2")

vf_cog <- c("SCO31", "SCO32",
            "DCO21", "DCO22",
            "DCO31", "DCO32",
            "TCO1", "TCO2")

# which haplotypes include each possible crossover
ycv_haps <- c("X1000", "X0111",
              "X0100", "X1011",
              "X0110", "X1001",
              "X0101", "X1010")

cvv_haps <- c("X0011", "X1100",
              "X0100", "X1011",
              "X0010", "X1101",
              "X0101", "X1010")

vf_haps <- c("X0001", "X1110",
             "X0110", "X1001",
             "X0010", "X1101",
             "X0101", "X1010")

# haplotypes containing each mutation
haps_y <- stri_extract(haplotypes, regex = "^1...") %>% na.omit
haps_cv <- stri_extract(haplotypes, regex = "^.1..") %>% na.omit
haps_v <- stri_extract(haplotypes, regex = "^..1.") %>% na.omit
haps_f <- stri_extract(haplotypes, regex = "^...1") %>% na.omit

# haplotypes containing each wildtype phenotype
haps_wt_body <- stri_extract(haplotypes, regex = "^0...") %>% na.omit
haps_wt_wing <- stri_extract(haplotypes, regex = "^.0..") %>% na.omit
haps_wt_eyec <- stri_extract(haplotypes, regex = "^..0.") %>% na.omit
haps_wt_hair <- stri_extract(haplotypes, regex = "^...0") %>% na.omit

