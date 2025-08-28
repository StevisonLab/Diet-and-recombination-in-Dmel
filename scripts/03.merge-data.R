##### Merge phenotype data with F1 cross information
##        
##    Uses:
##        From haplotype-groups.R
##            haplotypes
##        From main script:
##            data
##        From backcross-cleanup.R
##            backcross
##
##    Creates:
##        merged
##        by_vial
##        by_vialday
##        by_maternalvialday
##        by_treatstrain
##        summed


# Merge phenotyping data with F1 cross information
merged <- merge(data, backcross, by.x = "vial_num", by.y = "CrossID", all=T)

# remove NAs until all data collected
merged <- na.omit(merged)

#add adjusted age column
merged$adj_age=ifelse(merged$vial_letter=="A",merged$avg_age,ifelse(merged$vial_letter=="B",merged$avg_age+2,ifelse(merged$vial_letter=="C",merged$avg_age+4,ifelse(merged$vial_letter=="D",merged$avg_age+6,"error"))))
merged$adj_age=as.numeric(merged$adj_age)

#check age distribution
ggplot(data=merged,aes(x=num_offspring,y=min_age))+geom_point()+facet_wrap(~PaternalStock*Treatment)

by_vialday <- aggregate(merged[, c(haplotypes, "num_offspring")], by = list(merged$vial_num, merged$vial_letter,merged$adj_age), sum)
names(by_vialday)[1:3] <- c("vial_num", "vial_letter","adj_age")
by_vialday <- merge(by_vialday, backcross, by.x = "vial_num", by.y = "CrossID")


### Combine rows from the same vial
by_vial <- aggregate(merged[, c(haplotypes, "num_offspring","adj_age")],
                     by = list(merged$vial_num), sum)
names(by_vial)[1] <- c("vial_num")
by_vial <- merge(by_vial, backcross, by.x = "vial_num", by.y = "CrossID")


### Combine rows with the same maternal vial and day
#by_daymaternalvial <- aggregate(merged, by = list(merged$MaternalVial, merged$vial_letter), sum)
by_daymaternalvial <- summaryBy(
  merged[names(merged) %in% haplotypes] ~ merged$MaternalVial + merged$vial_letter,
  data = merged, id =  merged$average_age, FUN = sum, keep.names = TRUE)

names(by_daymaternalvial)[1:2] <- c("MaternalVial", "vial_letter")

# Sum phenotype counts for all vials with the same treatment and strain
by_treatstrain <- aggregate(merged[names(merged) %in% c(haplotypes, "num_offspring")],
                            by = list(merged$Treatment, merged$PaternalStock), sum)
by_treatstrain$num_M <- rowSums(by_treatstrain[, names(by_treatstrain) %in% haplotypes])

names(by_treatstrain)[1:2] <- c("Treatment", "PaternalStock")


### Total number of male flies observed with each haplotype
summed <- colSums(by_treatstrain[,names(by_treatstrain) %in% haplotypes])

write.csv(by_treatstrain[,c("Treatment", "PaternalStock", "num_M", "num_offspring")], "output/pheno-samp-size-by-treatment.csv", row.names = FALSE)
