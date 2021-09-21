rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)
setwd("/mnt/work/yunsung.lee/06.ART_DNAm_BW/")

### Load info
cat("Load info:\n")
info=as.data.frame(feather::read_feather(paste("data/info_w_variables.feather",sep="")))
cat("Print dim(info):",dim(info),"\n")

### Exclude samples whose exposure is NA.
exposure=commandArgs(TRUE)[1];cat("Print exposure:",exposure,"\n")
info=subset(info,is.na(info[,exposure])==F)
table(info[,exposure],useNA = "always")

### Load expr
expr=as.data.frame(feather::read_feather(paste("/mnt/work/yunsung.lee/data/bmat_child_autosomal_BMIQ.feather",sep="")))

### Transform expr to matrix.
rownames(expr)=expr$Sample_ID;expr$Sample_ID=c();expr=as.matrix(expr)

### Set x_Sample_ID
x_Sample_ID="c_Sample_ID"

### Exclude samples from expr
expr=expr[which(rownames(expr)%in% info[,x_Sample_ID]),]

### Align expr according to info
expr=expr[match(info[,x_Sample_ID],rownames(expr)),]
table(info[,x_Sample_ID] == rownames(expr),useNA = "always")

### Transform beta values to M values
expr=log2(expr/(1-expr))
print(head(info[,c(x_Sample_ID,exposure)]))
print(head(expr[,1:4]))

### Save
feather::write_feather(info,paste("data/info_child_",exposure,".feather",sep=""))
feather::write_feather(data.frame(expr),paste("data/expr_child_",exposure,".feather",sep=""))
