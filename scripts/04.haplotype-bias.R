##### Haplotype bias calculations and plotting
##
##    Defines the following functions:
##        calc.skew()
##        calc.skew.m()
##        skew.direction()
##        binom.pval()
##        haplotype.bias.pval()
##        hap.bias.pval.m()
##        calc.biases.data.frame()
##        melt.biases.data.frame()
##        plot.hap.bias()
##
##    Needs:
##        From haplotype-groups.R
##            haps_nco, haps_sco1 ... haps_dco3, haps_tco
##        From merge-data.R
##            by_vialday, by_daymaternalvial
##        
##    Creates:
##        --
##        
##    NOTES: 
##    - Warnings about rows missing values (`geom_point()`) are from rows with
##        no flies with either haplotype in that pair
##    - Warnings about rows missing values (`geom_text()`) are haplotypes that
##        had no entries with significant skew toward them when adding N labels
##        
##        



##### DEFINE FUNCTIONS ---------------------------------------------------------

## Calculate skew for a single pair (i.e., two sums)
##      @param  input       named vector with counts of flies summed by haplotype
##      @param  cog_pair    character vector with names corresponding to names
##                            of two values from input to compare
##      @return             0.5 if no skew,
##                            0 if skew toward crossover group "1" (left)
##                            1 if skew toward crossover group "2" (right)
calc.skew <- function(input, cog_pair) {
  . <- input[cog_pair[2]] / sum(input[cog_pair[1]],input[cog_pair[2]])
  .[,1]
  
  ## Alternate method of calculating, where the smaller is divided by the bigger
  #pmin(input[,cog_pair[1]],input[,cog_pair[2]]) / pmax(input[,cog_pair[1]],input[,cog_pair[2]])
}

## Calculate skew of two columns for each row in a data.frame using mapply
##      @param  input       data.frame with columns containing counts of flies
##      @param  cog_pair    character vector of length 2 (not enforced, anything
##                            longer is ignored) with names corresponding to
##                            names of two count columns from input to compare
##      @return             0.5 if no skew,
##                            0 if skew toward crossover group "1" (left)
##                            1 if skew toward crossover group "2" (right)
calc.skew.m <-  function(input, cog_pair) {
  mapply(function(x, y) x / sum(x + y), input[,cog_pair[2]], input[,cog_pair[1]])
}

## Tell which haplotype the pair is skewed toward for a single pair
##      @param  input       data.frame with columns containing counts of flies
##      @param  cog_pair    character vector of length 2 (not enforced, anything
##                            longer is ignored) with names corresponding to
##                            names of two count columns from input to compare
skew.direction <- function(input, cog_pair) {
  if (input[cog_pair[1]] > input[cog_pair[2]]) { cog_pair[1] } else {cog_pair[2]}
}

## Return the p-value from a binomial test comparing m and t
##      @param    m         the first number to be compared
##      @param    t         the second number to be compared
##      @return             the p-value
binom.pval <- function(m, t) 
{
  binom.test(c(m, t), p = 0.5)$p.value
}

## Calculate a p-value comparing a crossover pair using a binomial test. Can use
## crossover group names or haplotype strings as long as they are the same
## between the two parameters.
##      @param    input     data frame containing sums for each haplotype
##      @param    cog_pair  The names of the columns in the data.frame to 
##                            compare. Only the first 2 elements are used.
hap.bias.pval <- function(input, cog_pair)
{
  binom.pval(input[cog_pair[1]], input[cog_pair[2]])
  #paste(pv)
}

## Calculate p-values using a binomial test on the number of flies of both
## haplotypes in the given crossover pair for each row in a data frame using
## mapply. Can use crossover group names or haplotype strings as long as they
## are the same between the two parameters.
##      @param    input     a data.frame containing the haplotype counts
##      @param    cog_pair  The names of the columns in the data.frame to
##                            compare. Only the first 2 elements are used.
hap.bias.pval.m <- function(input, hap_pair)
{
  mapply(binom.pval, m=input[,hap_pair[1]], t=input[,hap_pair[2]])
}

## Return a copy of the input data frame with 16 new columns corresponding to a
## p-value from a binomial test and a skew for each of the 8 haplotype pairs.
## Depends on the names of columns with counts of flies matching the naming of
## haplotypes from haplotype-groups.R
##      @param    input   data.frame containing the haplotype counts, with a
##                          column for each haplotype
calc.biases.data.frame <- function(input) {
  . <- input
  
  .$nco_pval <- hap.bias.pval.m(., haps_nco)
  .$sco1_pval <- hap.bias.pval.m(., haps_sco1)
  .$sco2_pval <- hap.bias.pval.m(., haps_sco2)
  .$sco3_pval <- hap.bias.pval.m(., haps_sco3)
  .$dco1_pval <- hap.bias.pval.m(., haps_dco1)
  .$dco2_pval <- hap.bias.pval.m(., haps_dco2)
  .$dco3_pval <- hap.bias.pval.m(., haps_dco3)
  .$tco_pval <- hap.bias.pval.m(., haps_tco)
  
  .$nco_skew <- calc.skew.m(., haps_nco)
  .$sco1_skew <- calc.skew.m(., haps_sco1)
  .$sco2_skew <- calc.skew.m(., haps_sco2)
  .$sco3_skew <- calc.skew.m(., haps_sco3)
  .$dco1_skew <- calc.skew.m(., haps_dco1)
  .$dco2_skew <- calc.skew.m(., haps_dco2)
  .$dco3_skew <- calc.skew.m(., haps_dco3)
  .$tco_skew <- calc.skew.m(., haps_tco)
  
  .
}

## Melt the data frame to the format required by plot.hap.bias()
##      @param    input   data.frame of haplotpye counts to melt
##      @param    id.vars character vector of column names to group the data by
melt.biases.data.frame <- function(input, id.vars) {
  ## Melt the data to get crossover_group as the variable and skew as the value
  melted_skew <- melt(input,
                      id.vars = id.vars,
                      measure.vars = measure_skew,
                      variable.name = "crossover_group",
                      value.name = "skew")
  ### Remove the _skew suffix from all entries in the crossover_group column to
  ##  make merging easier
  melted_skew$crossover_group <- gsub("_skew", "", melted_skew$crossover_group)
  
  ### Melt the data to get crossover_group as the variable and pval as the value
  melted_pval <- melt(input,
                      id.vars = id.vars,
                      measure.vars = measure_pval,
                      variable.name = "crossover_group",
                      value.name = "pval")
  ### Remove the _pval suffix from all entries in the crossover_group column to
  ##  make merging easier
  melted_pval$crossover_group <- gsub("_pval", "", melted_pval$crossover_group)
  
  ### Merge the two melted data frames to create one data frame with both pvals
  ##  and skew for each combination of id.vars
  melted <- merge(melted_skew, melted_pval)
  ### Making crossover_group a factor vector
  melted$crossover_group <- factor(melted$crossover_group,
                                   levels=rev(c("nco", "sco1", "sco2", "sco3",
                                                "dco1", "dco2", "dco3", "tco")))
  melted
}

## Create and return a ggplot object plotting haplotype bias
##
##      @param  data      melted data.frame with the names of the pairs of
##                          haplotypes in a factor vector called 
##                          "crossover_group" and the skew values in a numeric 
##                          vector called "skew" which has values between 0 and 
##                          1, with 0.5 being no skew, values below 0.5 being 
##                          skew to the left, and values above 0.5 being skew to
##                          the right
##      @param  left.axis.labels    the haplotypes on the left y-axis
##      @param  right.axis.labels   the haplotypes on the right y-axis
##      @param  n.label.position    -- not currently in use --
##                                  The position for the red numbers depicting N
##                                    for each haplotype. If not included, N
##                                    labels are not added without warning
##                                    - "inner": either side of 0.5 
##                                    - "outer": close to 0 and 1
plot.hap.bias <- function(data,
                          n.label.position = "inner",
                          left.axis.labels = axis_left,
                          right.axis.labels = axis_right) {
  ## All rows that have a skew to the left (as plotted)
  left_skew <- na.omit(data[data[,"skew"] < 0.5,])
  left_skew_n <- tapply(left_skew[,"skew"], left_skew[,"crossover_group"], FUN = length)
  ## left_skew_p <- round(tapply(left_skew[,"skew"], left_skew[,"crossover_group"], FUN = length)/left_skew_n,2)
  ## Add zeros so all rows have an N value (not required -- otherwise leaves a blank)
  #left_skew_n[is.na(left_skew_n)] <- 0
  
  ## All rows that have a skew to the right (as plotted)
  right_skew <- na.omit(data[data[,"skew"] > 0.5,])
  right_skew_n <- tapply(right_skew[,"skew"], right_skew[,"crossover_group"], FUN = length)
  ## right_skew_p <- round(tapply(right_skew[,"skew"], right_skew[,"crossover_group"], FUN = length)/right_skew_n,2)
  ## Add zeros so all rows have an N value (not required -- otherwise leaves a blank)
  #right_skew_n[is.na(right_skew_n)] <- 0
  
  ## Convert to numeric to enable plotting as a continuous variable 
  ## Enables use of secondary axis with opposite haplotypes
  data[,"crossover_group"] <- as.numeric(data[,"crossover_group"])
  
  ## Create and return a ggplot object plotting haplotype bias
  ## NOTE: x- and y-axes are swapped at the end, so x is y and y is x
  . <- ggplot(data = data) +
    ## Line at 0.5 represents no bias
    geom_hline(yintercept = 0.5, color = "red") +
    ## Plot each point with an opacity of 25%
    geom_jitter(aes(x=crossover_group, y=skew),
                alpha = 0.25) +
    ## Using continuous scale to represent categorical crossover pair because
    #     of sec.axis argument, which adds a secondary axis. dup_axis() duplicates
    #     the original (left) axis and the labels argument overwrites the labels
    #     with the opposite haplotypes (right)
    scale_x_continuous(breaks = 1:8,
                       limits = c(0.5, 8.5),
                       expand = c(0, 0),
                       labels = axis_left,
                       sec.axis = dup_axis(labels = axis_right)) +
    ## Change limits and padding on sides to make room for N labels
    scale_y_continuous(limit = c(-0.1, 1.1), expand = c(0, 0)) +
    ## Add labels to show how many points are on each side of the red line
#    annotate("text", label = left_skew_n, x = 1:8, y = -0.05, color = "red", fontface = "bold") +
#    annotate("text", label = right_skew_n, x = 1:8, y = 1.05, color = "red", fontface = "bold") +
    ## Select and alter the theme
    theme_light() +
    theme(axis.title = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) +
    ## Switch the x- and y-axes for easier interpretation
    coord_flip()
  
  # ## Add labels to show how many points are on each side
  # # Label position depends on n.label.position
  # if (n.label.position == "inner") {
  #   ## Either side of 0.5
  #   . <- . +
  #     annotate("text", label = left_skew_n, x = 1:8, y = 0.45, color = "red", fontface = "bold") +
  #     annotate("text", label = right_skew_n, x = 1:8, y = 0.55, color = "red", fontface = "bold")
  # } else if (n.label.position == "outer") {
  #   ## Close to 0 and 1
  #   . <- . +
  #     annotate("text", label = left_skew_n, x = 1:8, y = 0.1, color = "red", fontface = "bold") +
  #     annotate("text", label = right_skew_n, x = 1:8, y = 0.9, color = "red", fontface = "bold")
  # }
  
  ## Otherwise don't add N labels
  
  ## Return the ggplot object
  .
}



########## DEFINE VECTORS USED FOR ALL DATASETS --------------------------------

### Names of columns containing p-values and skews to use for melting -- matches
##  names used by calc.biases.data.frame() for the new columns
measure_pval = c("nco_pval",
                 "sco1_pval",
                 "sco2_pval",
                 "sco3_pval",
                 "dco1_pval",
                 "dco2_pval",
                 "dco3_pval",
                 "tco_pval")

measure_skew = c("nco_skew",
                 "sco1_skew",
                 "sco2_skew",
                 "sco3_skew",
                 "dco1_skew",
                 "dco2_skew",
                 "dco3_skew",
                 "tco_skew")

### Haplotypes on the left y-axis
##  NOTE: each pair has a numeric value for positioning on the plot -- starting 
##        with the nco pair, which has a numeric value of 8.0, and counting down 
##        to the tco pair, which has a numeric value of 1.0
axis_left <- rev(c(haps_nco[1],
                   haps_sco1[1],
                   haps_sco2[1],
                   haps_sco3[1],
                   haps_dco1[1], 
                   haps_dco2[1], 
                   haps_dco3[1], 
                   haps_tco[1]))

### Haplotypes on the right y-axis -- order matched with axis_left
##  NOTE: each pair has a numeric value for positioning on the plot -- starting 
##        with the nco pair, which has a numeric value of 8.0, and counting down 
##        to the tco pair, which has a numeric value of 1.0
axis_right <- rev(c(haps_nco[2],
                    haps_sco1[2],
                    haps_sco2[2],
                    haps_sco3[2],
                    haps_dco1[2],
                    haps_dco2[2],
                    haps_dco3[2],
                    haps_tco[2]))



########## PREPARE DATA FOR PLOTTING -------------------------------------------

##### By vial and day
## Create a copy of by_vialday with 16 new columns, corresponding to a p-value
## from a binomial test and a skew for each of the 8 haplotype pairs
by_vialday_hapbias <- calc.biases.data.frame(by_vialday)

### Melt data to form for plotting
melted_vialday <- melt.biases.data.frame(by_vialday_hapbias, c("vial_num", "vial_letter","Treatment","PaternalStock"))

### Subset entries with significant skew
melted_vialday_sig <- melted_vialday[melted_vialday$pval < 0.05,]


##### By day and maternal vial
## Create a copy of by_daymaternalvial with 16 new columns, corresponding to a
## p-value from a binomial test and a skew for each of the 8 haplotype pairs
by_daymat_hapbias <- calc.biases.data.frame(by_daymaternalvial)

### Melt data to form for plotting
melted_daymat <- melt.biases.data.frame(by_daymat_hapbias, c("MaternalVial", "vial_letter"))

### Subset entries with significant skew
melted_daymat_sig <- melted_daymat[melted_daymat$pval < 0.05,]


####### TEST #######
### By treatment and paternal stock
test_calc <- calc.biases.data.frame(by_treatstrain)
test_melt <- melt.biases.data.frame(test_calc, c("Treatment", "PaternalStock"))
test_plot <- plot.hap.bias(test_melt)
#test_plot


########## PLOT ----------------------------------------------------------------
##
### Set what data frame to use for plotting
##  Enables easy swapping of data to plot to try different things
pd <- melted_vialday
pd_sig <- melted_vialday_sig

# pd <- melted_daymat
# pd_sig <- melted_daymat_sig


### Plot haplotype bias for all days
### Only rows with significant skew
hap_bias_plot_sig <- plot.hap.bias(pd_sig, "inner")
plot_grid(hap_bias_plot_sig)

#### Plot haplotype bias by treatment by strain
## Rows for 42
hap_bias_plot_sig_0.5x <- plot.hap.bias(pd_sig[pd_sig$Treatment == "0.5x" & pd_sig$PaternalStock=="42",])
hap_bias_plot_sig_1x <- plot.hap.bias(pd_sig[pd_sig$Treatment == "1x" & pd_sig$PaternalStock=="42",])
hap_bias_plot_sig_2x <- plot.hap.bias(pd_sig[pd_sig$Treatment == "2x" & pd_sig$PaternalStock=="42",])
plot_grid(hap_bias_plot_sig_0.5x,
          hap_bias_plot_sig_1x, 
          hap_bias_plot_sig_2x, 
          labels = "AUTO")
ggsave("images/hap_bias_42.png")

100*length(pd_sig$vial_num[pd_sig$Treatment == "0.5x" & pd_sig$PaternalStock=="42"])/length(pd$vial_num[pd$Treatment == "0.5x" & pd$PaternalStock=="42" & !is.na(pd$skew)])
100*length(pd_sig$vial_num[pd_sig$Treatment == "1x" & pd_sig$PaternalStock=="42"])/length(pd$vial_num[pd$Treatment == "1x" & pd$PaternalStock=="42" & !is.na(pd$skew)])
100*length(pd_sig$vial_num[pd_sig$Treatment == "2x" & pd_sig$PaternalStock=="42"])/length(pd$vial_num[pd$Treatment == "2x" & pd$PaternalStock=="42" & !is.nan(pd$skew)])

## Rows for 217
hap_bias_plot_sig_0.5x <- plot.hap.bias(pd_sig[pd_sig$Treatment == "0.5x" & pd_sig$PaternalStock=="217",])
hap_bias_plot_sig_1x <- plot.hap.bias(pd_sig[pd_sig$Treatment == "1x" & pd_sig$PaternalStock=="217",])
hap_bias_plot_sig_2x <- plot.hap.bias(pd_sig[pd_sig$Treatment == "2x" & pd_sig$PaternalStock=="217",])
plot_grid(hap_bias_plot_sig_0.5x,
          hap_bias_plot_sig_1x, 
          hap_bias_plot_sig_2x, 
          labels = "AUTO")
ggsave("images/hap_bias_217.png")

100*length(pd_sig$vial_num[pd_sig$Treatment == "0.5x" & pd_sig$PaternalStock=="217"])/length(pd$vial_num[pd$Treatment == "0.5x" & pd$PaternalStock=="217" & !is.na(pd$skew)])
100*length(pd_sig$vial_num[pd_sig$Treatment == "1x" & pd_sig$PaternalStock=="217"])/length(pd$vial_num[pd$Treatment == "1x" & pd$PaternalStock=="217" & !is.na(pd$skew)])
100*length(pd_sig$vial_num[pd_sig$Treatment == "2x" & pd_sig$PaternalStock=="217"])/length(pd$vial_num[pd$Treatment == "2x" & pd$PaternalStock=="217" & !is.na(pd$skew)])


##### HAPLOTYPE BIAS ON SUMMED DATA --------------------------------------------
##
# hapBiasByTwoThings <- function(input, by) {
#   
#   for (i in 1:length(input[,1])) {
#     tmp1 <- input[i, names(input) %in% haplotypes]
#     b1 <- input[i, by[1]]
#     b2 <- input[i, by[2]]
#     summed <- colSums(tmp1[, names(tmp1) %in% haplotypes])
#     
#     ## Calculate haplotype bias
#     hap_bias_pvals <- data.frame(crossover_group = cog_pairs,
#                                  p.value = c(
#                                    hap.bias.pval(summed, haps_nco),
#                                    hap.bias.pval(summed, haps_sco1),
#                                    hap.bias.pval(summed, haps_sco2),
#                                    hap.bias.pval(summed, haps_sco3),
#                                    hap.bias.pval(summed, haps_dco1),
#                                    hap.bias.pval(summed, haps_dco2),
#                                    hap.bias.pval(summed, haps_dco3),
#                                    hap.bias.pval(summed, haps_tco)))
#     
#     hap_bias_pvals$p.value <- as.numeric(hap_bias_pvals$p.value)
#     hap_bias_pvals$significance <- eval.significance(hap_bias_pvals)
#     hap_bias_significant <- hap_bias_pvals[!(hap_bias_pvals$significance == "" |
#                                                is.na(hap_bias_pvals$significance)),]
#     
#     ## Calculate skew to evaluate data
#     skew <- data.frame(crossover_group = cog_pairs, 
#                        skew = c(calc.skew(summed, haps_nco),
#                                 calc.skew(summed, haps_sco1),
#                                 calc.skew(summed, haps_sco2),
#                                 calc.skew(summed, haps_sco3),
#                                 calc.skew(summed, haps_dco1),
#                                 calc.skew(summed, haps_dco2),
#                                 calc.skew(summed, haps_dco3),
#                                 calc.skew(summed, haps_tco)),
#                        skew_direction = c(skew.direction(summed, haps_nco),
#                                           skew.direction(summed, haps_sco1),
#                                           skew.direction(summed, haps_sco2),
#                                           skew.direction(summed, haps_sco3),
#                                           skew.direction(summed, haps_dco1),
#                                           skew.direction(summed, haps_dco2),
#                                           skew.direction(summed, haps_dco3),
#                                           skew.direction(summed, haps_tco)))
#     
#     #merge
#     hap_bias_summary <- merge(hap_bias_pvals, skew, by = "crossover_group")
#     name <- paste("output/hap_bias_summary", b1, b2, "csv", sep=".")
#     write.csv(hap_bias_summary, file = name, row.names = F, quote = F)
#   }
# }
# 
# hapBiasByTwoThings(by_daymaternalvial, c("MaternalVial", "vial_letter"))
# 
# 
# 
# # Crossover group pairs defined in haplotype-groups.R
# # Sum of males with each haplotype done in merge-data.R
# # summed <- colSums(by_treatstrain[,names(by_treatstrain) %in% haplotypes])
# 
# ## Calculate haplotype bias to evaluate data
# hap_bias_pvals <- data.frame(crossover_group = cog_pairs,
#                              p.value = c(
#                                 hap.bias.pval(summed, haps_nco),
#                                 hap.bias.pval(summed, haps_sco1),
#                                 hap.bias.pval(summed, haps_sco2),
#                                 hap.bias.pval(summed, haps_sco3),
#                                 hap.bias.pval(summed, haps_dco1),
#                                 hap.bias.pval(summed, haps_dco2),
#                                 hap.bias.pval(summed, haps_dco3),
#                                 hap.bias.pval(summed, haps_tco)))
# 
# hap_bias_pvals$p.value <- as.numeric(hap_bias_pvals$p.value)
# hap_bias_pvals$significance <- eval.significance(hap_bias_pvals)
# hap_bias_significant <- hap_bias_pvals[!(hap_bias_pvals$significance == "" |
#     is.na(hap_bias_pvals$significance)),]
# 
# ## Calculate skew to evaluate data
# skew <- data.frame(crossover_group = cog_pairs, 
#                    skew = c(calc.skew(summed, haps_nco),
#                             calc.skew(summed, haps_sco1),
#                             calc.skew(summed, haps_sco2),
#                             calc.skew(summed, haps_sco3),
#                             calc.skew(summed, haps_dco1),
#                             calc.skew(summed, haps_dco2),
#                             calc.skew(summed, haps_dco3),
#                             calc.skew(summed, haps_tco)),
#                    skew_direction = c(skew.direction(summed, haps_nco),
#                                       skew.direction(summed, haps_sco1),
#                                       skew.direction(summed, haps_sco2),
#                                       skew.direction(summed, haps_sco3),
#                                       skew.direction(summed, haps_dco1),
#                                       skew.direction(summed, haps_dco2),
#                                       skew.direction(summed, haps_dco3),
#                                       skew.direction(summed, haps_tco)))
# 
# #merge
# hap_bias_summary=merge(hap_bias_pvals,skew,by="crossover_group")
# 
# #examine skew by treatment, day and stock
# for (i in 1:length(by_treatstrain$Treatment)) {
#   tmp1=by_treatstrain[i,c(4:19)]
#   treatment=by_treatstrain[i,1]
#   stock=by_treatstrain[i,2]
#   summed <- colSums(tmp1[,names(tmp1) %in% haplotypes])
#   ## Calculate haplotype bias to evaluate data
#   hap_bias_pvals <- data.frame(crossover_group = cog_pairs,
#                                p.value = c(
#                                  hap.bias.pval(summed, haps_nco),
#                                  hap.bias.pval(summed, haps_sco1),
#                                  hap.bias.pval(summed, haps_sco2),
#                                  hap.bias.pval(summed, haps_sco3),
#                                  hap.bias.pval(summed, haps_dco1),
#                                  hap.bias.pval(summed, haps_dco2),
#                                  hap.bias.pval(summed, haps_dco3),
#                                  hap.bias.pval(summed, haps_tco)))
#   
#   hap_bias_pvals$p.value <- as.numeric(hap_bias_pvals$p.value)
#   hap_bias_pvals$significance <- eval.significance(hap_bias_pvals)
#   hap_bias_significant <- hap_bias_pvals[!(hap_bias_pvals$significance == "" |
#                                              is.na(hap_bias_pvals$significance)),]
#   
#   ## Calculate skew to evaluate data
#   skew <- data.frame(crossover_group = cog_pairs, 
#                      skew = c(calc.skew(summed, haps_nco),
#                               calc.skew(summed, haps_sco1),
#                               calc.skew(summed, haps_sco2),
#                               calc.skew(summed, haps_sco3),
#                               calc.skew(summed, haps_dco1),
#                               calc.skew(summed, haps_dco2),
#                               calc.skew(summed, haps_dco3),
#                               calc.skew(summed, haps_tco)),
#                      skew_direction = c(skew.direction(summed, haps_nco),
#                                         skew.direction(summed, haps_sco1),
#                                         skew.direction(summed, haps_sco2),
#                                         skew.direction(summed, haps_sco3),
#                                         skew.direction(summed, haps_dco1),
#                                         skew.direction(summed, haps_dco2),
#                                         skew.direction(summed, haps_dco3),
#                                         skew.direction(summed, haps_tco)))
#   
#   #merge
#   hap_bias_summary=merge(hap_bias_pvals,skew,by="crossover_group")
#   name=paste("output/hap_bias_summary",treatment,stock,"csv",sep=".")
#   write.csv(hap_bias_summary,file=name,row.names = F,quote = F)
# }
# 
# 
# 
# #might be nice to evaluate this more precisely within each vial and 2 day collection period
# #need to use the by_vialday dataframe
# 
# for (i in 1:length(by_daymaternalvial[,1])) {
#   tmp1 <- by_daymaternalvial[i, names(by_daymaternalvial) %in% haplotypes]
#   b1 <- by_daymaternalvial[i, "MaternalVial"] # maternal vial
#   b2 <- by_daymaternalvial[i, "vial_letter"] # vial_letter
#   summed <- colSums(tmp1[, names(tmp1) %in% haplotypes])
#   ## Calculate haplotype bias to evaluate data
#   hap_bias_pvals <- data.frame(crossover_group = cog_pairs,
#                                p.value = c(
#                                  hap.bias.pval(summed, haps_nco),
#                                  hap.bias.pval(summed, haps_sco1),
#                                  hap.bias.pval(summed, haps_sco2),
#                                  hap.bias.pval(summed, haps_sco3),
#                                  hap.bias.pval(summed, haps_dco1),
#                                  hap.bias.pval(summed, haps_dco2),
#                                  hap.bias.pval(summed, haps_dco3),
#                                  hap.bias.pval(summed, haps_tco)))
#   
#   hap_bias_pvals$p.value <- as.numeric(hap_bias_pvals$p.value)
#   hap_bias_pvals$significance <- eval.significance(hap_bias_pvals)
#   hap_bias_significant <- hap_bias_pvals[!(hap_bias_pvals$significance == "" |
#                                              is.na(hap_bias_pvals$significance)),]
#   
#   ## Calculate skew to evaluate data
#   skew <- data.frame(crossover_group = cog_pairs, 
#                      skew = c(calc.skew(summed, haps_nco),
#                               calc.skew(summed, haps_sco1),
#                               calc.skew(summed, haps_sco2),
#                               calc.skew(summed, haps_sco3),
#                               calc.skew(summed, haps_dco1),
#                               calc.skew(summed, haps_dco2),
#                               calc.skew(summed, haps_dco3),
#                               calc.skew(summed, haps_tco)),
#                      skew_direction = c(skew.direction(summed, haps_nco),
#                                         skew.direction(summed, haps_sco1),
#                                         skew.direction(summed, haps_sco2),
#                                         skew.direction(summed, haps_sco3),
#                                         skew.direction(summed, haps_dco1),
#                                         skew.direction(summed, haps_dco2),
#                                         skew.direction(summed, haps_dco3),
#                                         skew.direction(summed, haps_tco)))
#   
#   #merge
#   hap_bias_summary <- merge(hap_bias_pvals, skew, by="crossover_group")
#   name <- paste("output/hap_bias_summary", b1, b2, "csv", sep=".")
#   write.csv(hap_bias_summary, file=name, row.names = F, quote = F)
# }