rawbackcross <- read.csv("rawdata/F1Collections.csv", header = TRUE)

backcross <- rawbackcross[, c(1:3, 6:8)]#"ï..Maternal.Vial", "Treatment", "Paternal.Stock", "CrossID", "X..mothers")]
names(backcross) <- c("MaternalVial", "Treatment", "PaternalStock", "CrossID", "Maternal_Age","num_moms")
backcross <- na.omit(backcross)

backcross$MaternalVial <- as.factor(backcross$MaternalVial)
backcross$Treatment <- as.factor(tolower(backcross$Treatment))
backcross$PaternalStock <- as.factor(backcross$PaternalStock) 
backcross$CrossID <- as.factor(backcross$CrossID) 
backcross$num_moms <- as.numeric(backcross$num_moms)
backcross$Maternal_Age <- as.numeric(backcross$Maternal_Age)
backcross <- na.omit(backcross)

#remove vials with duplicated labels
#backcross$CrossID[duplicated(backcross$CrossID)]
dups=c(4,143:146)
backcross=backcross[!backcross$CrossID %in% dups,]

#aggregate vials
num_moms_calc <- aggregate(num_moms ~ MaternalVial + CrossID + Treatment + PaternalStock, data = backcross, sum)

min_age_calc <- aggregate(Maternal_Age ~ MaternalVial + CrossID + Treatment + PaternalStock, data = backcross, min)

max_age_calc <- aggregate(Maternal_Age ~ MaternalVial + CrossID + Treatment + PaternalStock, data = backcross, max)

age_merge1=merge(min_age_calc,max_age_calc,by=c("MaternalVial","CrossID","Treatment","PaternalStock"))

colnames(age_merge1)=c("MaternalVial","CrossID","Treatment","PaternalStock","min_age","max_age")

age_calc=merge(age_merge1,num_moms_calc,by=c("MaternalVial","CrossID","Treatment","PaternalStock"))

age_calc$age_diff=age_calc$max_age - age_calc$min_age

#remove vials with age range >2
age_calc=subset(age_calc,age_calc$age_diff<3)

#add column with the average age of each vial
age_calc$avg_age=(age_calc$max_age + age_calc$min_age)/2

#remove vials where min age is above 7
age_calc=subset(age_calc,age_calc$min_age<=7)

backcross=age_calc

