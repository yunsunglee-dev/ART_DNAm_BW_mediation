rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)
setwd("/mnt/work/yunsung.lee/06.ART_DNAm_BW")

role=commandArgs(TRUE)[1];cat("role:",role,"\n")
exposure=commandArgs(TRUE)[2];cat("exposure:",exposure,"\n")

### Load phenotypic data
cat("Load phenotypic data\n")
info=as.data.frame(feather::read_feather(paste("data/info_",role,"_",exposure,".feather",sep="")))

### Load DNAm data at the ART-associated CpGs
cat("Load DNAm data at the ART-associated CpGs\n")
expr=as.data.frame(feather::read_feather(paste("data/expr_",role,"_",exposure,".feather",sep="")))

### Load EWAS results.
ewas=as.data.frame(feather::read_feather(paste("result_data/ewas_",role,"_",exposure,".feather",sep="")))

### Select significant CpGs after Bonferroni correction.
#ewas=subset(ewas,ewas$p<0.05/nrow(ewas))

### Select significant CpGs after FDR 0.05
ewas=subset(ewas,ewas$q < 0.05)

cat("ewas dim:",dim(ewas))
head(ewas)

info=data.frame(info,expr[,ewas$Name])
feather::write_feather(info,paste("result_data/info_",role,"_",exposure,".feather",sep=""))
