rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)

role=commandArgs(TRUE)[1];cat("role:",role,"\n")
exposure=commandArgs(TRUE)[2];cat("exposure:",exposure,"\n")
method=commandArgs(TRUE)[3];cat("method:",method,"\n")
alpha=as.numeric(commandArgs(TRUE)[4]);cat("alpha:",alpha,"\n")
#role="child";exposure="fresh_frozen";method="fdr";alpha=0.05

### Load info
info=as.data.frame(feather::read_feather(paste("/mnt/work/yunsung.lee/06.ART_DNAm_BW/result_data/info_",role,"_",exposure,".feather",sep="")))
cat("Sample size:", nrow(info),"\n\n\n")
cat("Distribution of ",exposure,":\n")
print(table(info[,exposure],useNA="always"))

### Focus on the ART-associated CpGs
t1=readRDS(paste("/mnt/work/yunsung.lee/06.ART_DNAm_BW/result_data/list_",role,"_",exposure,"_",method,alpha,".RData",sep=""))
if(method=="bon"){
    t1=subset(t1,t1$p < alpha/nrow(t1))    
}else if(method=="fdr"){
    t1=subset(t1,t1$q < alpha)    
}    
target_cpgs=t1$Name
cat("The number of CpGs associated with ART and BW (no CpG*ART interaction):",length(target_cpgs),"\n\n\n")


### Start mediation analysis 
cat("Start mediation analysis\n")
cat("The number of target CpGs:",length(target_cpgs),"\n")

### Define adjusting variables
adj_vars=c("MORS_ALDER","mom_edu2","mom_edu3","mom_edu4",
            "parity_b","bmi_before","daily_before_during_preg1","daily_before_during_preg2","daily_before_during_preg3","kjonn_epic_b","folate1","folate2")
cat("Print adjusting variables:",adj_vars,"\n")

### Define c_AMP_Plate as a factor variable
info$c_AMP_Plate=as.factor(info$c_AMP_Plate)



### Define a work horse for bootstrapping.
boot=function(b){
    #b=1;one_cpg=target_cpgs[1]
    set.seed(b)
    info_boot=na.omit(info[sample(1:nrow(info),nrow(info),replace = T),c("VEKT",exposure,adj_vars,one_cpg,"c_AMP_Plate")])
    ### DNAm ~ Fresh/non + covariates (model for a total effect)
    myfor1=as.formula(paste("VEKT~",exposure,"+",paste(adj_vars,collapse = "+"),sep=""))
    model1=nlme::lme(myfor1, random=~1|c_AMP_Plate, data=info_boot)
    tab=summary(model1)$tTable
    beta_total=tab[rownames(tab)==exposure,1]
    ### DNAm ~ Fresh/non + CpG + covariates (model for a direct effect)
    myfor1=as.formula(paste("VEKT~",exposure,"+",one_cpg,"+",paste(adj_vars,collapse = "+"),sep=""))
    model1=nlme::lme(myfor1, random=~1|c_AMP_Plate, data=info_boot)
    tab=summary(model1)$tTable
    beta_dir=tab[rownames(tab)==exposure,1]    
    return(c(beta_total,beta_dir))
}

### Run difference models over all the CpGs differentially methylated between non-ART and fresh.
cat("Run difference models over all the CpGs of interest.\n")
result=list()
for(i in 1:length(target_cpgs)){
    one_cpg=target_cpgs[i]
    cat("Start ",one_cpg,"\n")
    result[[i]]=unlist(parallel::mclapply(1:5000,boot,mc.cores = 12))    ##30
    cat("Finish ",one_cpg,"\n\n")
}
names(result)=target_cpgs
result=t(do.call("rbind",result))

### Save
cat("Save\n")
feather::write_feather(data.frame(result),
                       paste("/mnt/work/yunsung.lee/06.ART_DNAm_BW/result_data/sobel_test_",role,"_",exposure,"_",method,alpha,".feather",sep=""))
