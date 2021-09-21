rm(list=ls());gc(TRUE)
setwd("/mnt/work/yunsung.lee/06.ART_DNAm_BW/zscore")
options(stringsAsFactors = F)

role=commandArgs(TRUE)[1];cat("role:",role,"\n")
exposure=commandArgs(TRUE)[2];cat("exposure:",exposure,"\n")
method=commandArgs(TRUE)[3];cat("method:",method,"\n")
alpha=as.numeric(commandArgs(TRUE)[4]);cat("alpha:",alpha,"\n")
outcome=commandArgs(TRUE)[5];cat("outcome:",outcome,"\n")
#role="child";exposure="non_fresh";method="bon";alpha=0.05;outcome="VEKT_s"

### Load info_role_... data.
info=as.data.frame(feather::read_feather(paste("data/info_",role,"_",exposure,".feather",sep="")))
info$c_AMP_Plate=as.factor(info$c_AMP_Plate)

adj_vars=c("MORS_ALDER","mom_edu2","mom_edu3","mom_edu4",
            "parity_b","bmi_before","daily_before_during_preg1","daily_before_during_preg2","daily_before_during_preg3","kjonn_epic_b","folate1","folate2")


### Load ewas of ART and BW summary statistics
ewas=readRDS(paste("result_data/list_",role,"_",exposure,"_",method,alpha,".RData",sep=""))

### Select significant CpGs.
if(method=="bon"){
    ewas=subset(ewas,ewas$p<alpha/nrow(ewas))    
}else if(method=="fdr"){
    ewas=subset(ewas,ewas$q<alpha)    
}

### Calculate the composite effect
### Define a work horse for bootstrapping.
boot=function(b){
    #b=1;one_cpg=target_cpgs[1]
    set.seed(b)
    info_boot=na.omit(info[sample(1:nrow(info),nrow(info),replace = T),c(outcome,exposure,adj_vars,cpgs_interest,"c_AMP_Plate")])
    ### DNAm ~ Fresh/non + covariates (model for a total effect)
    myfor1=as.formula(paste(outcome,"~",exposure,"+",paste(adj_vars,collapse = "+"),sep=""))
    model1=nlme::lme(myfor1, random=~1|c_AMP_Plate, data=info_boot)
    tab=summary(model1)$tTable
    beta_total=tab[rownames(tab)==exposure,1]
    ### DNAm ~ Fresh/non + CpG + covariates (model for a direct effect)
    myfor1=as.formula(paste(outcome,"~",exposure,"+",paste(cpgs_interest,collapse="+"),"+",paste(adj_vars,collapse = "+"),sep=""))
    model1=nlme::lme(myfor1, random=~1|c_AMP_Plate, data=info_boot)
    tab=summary(model1)$tTable
    beta_dir=tab[rownames(tab)==exposure,1]    
    return(c(beta_total,beta_dir))
}

cpgs_interest=ewas$Name
cat("Print target CpGs:",cpgs_interest,"\n")
length(cpgs_interest)
result=parallel::mclapply(1:5000,boot,mc.cores = 3)
result=do.call("rbind",result)

### Save
cat("Save\n")
feather::write_feather(data.frame(result),paste("result_data/sobel_test_combined_",role,"_",exposure,"_",method,alpha,".feather",sep=""))