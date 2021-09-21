rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)
setwd("/mnt/work/yunsung.lee/06.ART_DNAm_BW/")

role=commandArgs(TRUE)[1];cat("role:",role,"\n")
exposure=commandArgs(TRUE)[2];cat("exposure:",exposure,"\n")
method=commandArgs(TRUE)[3];cat("method:",method,"\n")
alpha=as.numeric(commandArgs(TRUE)[4]);cat("alpha:",alpha,"\n")
#role="child";exposure="fresh_frozen";method="fdr";alpha=0.05

### Load info
info=as.data.frame(feather::read_feather(paste("result_data/info_",role,"_",exposure,".feather",sep="")))
cat("Sample size:", nrow(info),"\n\n\n")
cat("Distribution of exposure:\n");print(table(info[,exposure],useNA="always"))

### Focus on the ART-associated CpGs
ewas=as.data.frame(feather::read_feather(paste("result_data/ewas_",role,"_",exposure,".feather",sep="")));
cat("EWAS summary statistics\nThe number of CpGs:",nrow(ewas),"\n")
if(method=="bon"){sig=subset(ewas,ewas$p < alpha/nrow(ewas))}else if(method=="fdr"){sig=subset(ewas,ewas$q < alpha)}
cat("The number of CpGs after ",toupper(method)," ",alpha,":",nrow(sig),"\n\n\n")

### Start mediation analysis 
cat("Start mediation analysis\n")
target_cpgs=sig$Name
cat("The number of target CpGs:",length(target_cpgs),"\n")

### Define adjusting variables
adj_vars=c("MORS_ALDER","mom_edu2","mom_edu3","mom_edu4",
           "parity_b","bmi_before",
           "daily_before_during_preg1","daily_before_during_preg2","daily_before_during_preg3",
           "kjonn_epic_b","folate1","folate2")
cat("Print adjusting variables:",adj_vars,"\n")
### Define c_AMP_Plate as a factor variable
info$c_AMP_Plate=as.factor(info$c_AMP_Plate)

### Fit regressions with intereactions: BW ~ ART + CpG + ART*CpG + covariates
work_horse=function(i){
    model3=nlme::lme(as.formula(paste("VEKT~",target_cpgs[i],"+",exposure,"+",target_cpgs[i],"*",exposure,"+",
                                      paste(adj_vars,collapse = "+"),
                         sep="")),random=~1|c_AMP_Plate,data=info,na.action="na.exclude")
    tab=summary(model3)$tTable
    output=rbind(as.numeric(tab[which(rownames(tab)==exposure),]),
          as.numeric(tab[which(rownames(tab)==target_cpgs[i]),]),
          as.numeric(tab[which(rownames(tab)==paste(target_cpgs[i],":",exposure,sep="")),]))
    rownames(output)=c(exposure,target_cpgs[i],"interaction")
    colnames(output)=colnames(tab)
    return(output)
}

result=parallel::mclapply(1:length(target_cpgs),work_horse,mc.cores = 3);
names(result)=target_cpgs

result_dat=do.call("rbind",lapply(result,FUN=function(i){
    c(as.numeric(i[2,]),as.numeric(i[3,]))
}))
result_dat=data.frame(Name=rownames(result_dat),result_dat)
colnames(result_dat)=c("Name",
                       paste(c("b","se","df","t","p"),c("_bw_cpg"),sep=""),
                       paste(c("b","se","df","t","p"),c("_bw_inter"),sep=""))
result_dat$q_bw_cpg=p.adjust(result_dat$p_bw_cpg,"BH")
result_dat$q_bw_inter=p.adjust(result_dat$p_bw_inter,"BH")

### Select CpGs without interactions
if(method=="bon"){
    result_dat_no_inter=subset(result_dat,result_dat$p_bw_inter >= alpha/nrow(result_dat))
}else if(method=="fdr"){
    result_dat_no_inter=subset(result_dat,result_dat$q_bw_inter >= alpha)
}
cat("The number of cases with ",toupper(method)," significant interactions:",
    nrow(result_dat)-nrow(result_dat_no_inter),"\n")
### Save
saveRDS(result_dat,paste("result_data/list_inter_",role,"_",exposure,"_",method,alpha,".RData",sep=""))

### Define CpGs without interaction.
target_cpgs=result_dat_no_inter$Name

### Fit regressions without any interactions: BW ~ CpG + covariates
work_horse=function(i){
    model3=nlme::lme(as.formula(paste("VEKT~",target_cpgs[i],"+",
                                      paste(adj_vars,collapse = "+"),
                         sep="")),random=~1|c_AMP_Plate,data=info,na.action="na.exclude")
    tab=summary(model3)$tTable
    return(tab[rownames(tab)==target_cpgs[i],])
}

result_dat_no_inter=parallel::mclapply(1:length(target_cpgs),work_horse,mc.cores=3)
names(result_dat_no_inter)=target_cpgs;t1=do.call("rbind",result_dat_no_inter)
result_dat_no_inter=data.frame(Name=rownames(t1),t1);
colnames(result_dat_no_inter)=c("Name","b","se","df","t","p")

cat("Print BW-associated CpGs: \n")
if(method=="bon"){
    ### Bonferroni
    cat("Print significant CpGs after the Bonferroni correction at",alpha,":\n")
    print(subset(result_dat_no_inter,result_dat_no_inter$p < alpha/nrow(result_dat_no_inter)))    
}else if(method=="fdr"){
    ### Benjamini-Hochberg procedure
    result_dat_no_inter$q=p.adjust(result_dat_no_inter$p,"BH")
    cat("Print significant CpGs after FDR < ",alpha,":\n")
    print(subset(result_dat_no_inter,result_dat_no_inter$q < alpha))    
}
### Save
saveRDS(result_dat_no_inter,paste("result_data/list_",role,"_",exposure,"_",method,alpha,".RData",sep=""))
