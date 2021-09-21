rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)
setwd("/mnt/work/yunsung.lee/06.ART_DNAm_BW/")

### Load info
info=as.data.frame(feather::read_feather("/mnt/work/yunsung.lee/data/info.feather"))

### Remove TED samples from non-ART newborns.
info=info[-which(info$ANY_ART==0 & info$TED_flag==1),]

### Add a variable for folate intake
### Ellen's code
#/*Generating  a categorical variable indicating folate intake never, before pregnancy, before Q1, both before pregnancy and before Q1 */
#gen folate_check=0
#replace folate_check=1 if (AA940==1 | AA941==1 | AA942==1 ) & (AA946==. & AA943==. & AA944==. & AA945==.) 
#replace folate_check=2 if  (AA940==. & AA941==. & AA942==. ) & (AA946==1 | AA943==1 | AA944==1 | AA945==1) 
#replace folate_check=3 if (AA940==1 | AA941==1 | AA942==1 ) & (AA946==1 | AA943==1 | AA944==1 | AA945==1) 
#tab1 folate_check

cat("Define a variable for folate intake\n")
info$folate_check=0
info$folate_check[which((info$AA940==1|info$AA941==1|info$AA942==1) & 
                       (is.na(info$AA943)==T &is.na(info$AA944)==T &is.na(info$AA945)==T &is.na(info$AA946)==T))]=1

info$folate_check[which((is.na(info$AA940)==T|is.na(info$AA941)==T|is.na(info$AA942)==T) & 
                       (info$AA943==1|info$AA944==1|info$AA945==1|info$AA946==1))]=2

info$folate_check[which((info$AA940==1|info$AA941==1|info$AA942==1) & 
                       (info$AA943==1|info$AA944==1|info$AA945==1|info$AA946==1))]=3


### Ellen's code
#/* very few with folate only before pregnancy (n=25), will therefore make new variable with folate intake never (or only before pregnancy), before Q1,and both before + Q1 */
#cap drop folate
#gen folate=0
#replace folate=1 if (AA940==. & AA941==. & AA942==. ) & (AA946==1 | AA943==1 | AA944==1 | AA945==1)
#replace folate=2 if (AA940==1 | AA941==1 | AA942==1 ) & (AA946==1 | AA943==1 | AA944==1 | AA945==1) 
#label define folatelb 0"Never or before pregnancy" 1"In pregnancy" 2"Before + in pregnancy"
#label values folate folatelb
#tab1 folate

info$folate=0
info$folate[which(
(is.na(info$AA940)==T&is.na(info$AA941)==T&is.na(info$AA942)==T) &
(info$AA943==1|info$AA944==1|info$AA945==1|info$AA946==1)
)]=1

info$folate[which(
(info$AA940==1|info$AA941==1|info$AA942==1) &
(info$AA943==1|info$AA944==1|info$AA945==1|info$AA946==1)
)]=2

### Generate dummy variables
info$folate1=ifelse(is.na(info$folate)==T,NA,ifelse(info$folate==1,1,0))
info$folate2=ifelse(is.na(info$folate)==T,NA,ifelse(info$folate==2,1,0))

cat("Distribution of folate\n")
print(table(info$folate,useNA = "always"))
cat("Distribution of folate1\n")
print(table(info$folate1,useNA = "always"))
cat("Distribution of folate2\n")
print(table(info$folate2,useNA = "always"))

### Select samples that passed the QC.
info=subset(info,info$c_Sample_ID %in% readRDS("data/c_Sample_ID.RData"))

### Save
feather::write_feather(info,paste("/mnt/work/yunsung.lee/06.ART_DNAm_BW/data/info_w_variables.feather",sep=""))
