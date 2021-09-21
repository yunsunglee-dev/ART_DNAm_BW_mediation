rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)
setwd("/mnt/work/yunsung.lee/06.ART_DNAm_BW")

role=commandArgs(TRUE)[1];cat("role:",role,"\n")
exposure=commandArgs(TRUE)[2];cat("exposure:",exposure,"\n")

info=as.data.frame(feather::read_feather(paste("data/info_",role,"_",exposure,".feather",sep="")))
expr=as.data.frame(feather::read_feather(paste("data/expr_",role,"_",exposure,".feather",sep="")))

if(role=="child"){
    adj_vars=c("MORS_ALDER",
               "daily_before_during_preg1",
               "daily_before_during_preg2",
               "daily_before_during_preg3",
               "bmi_before","parity_b","kjonn_epic_b")

}else if(role=="mother"){
    adj_vars=c("MORS_ALDER",
               "daily_before_during_preg1",
               "daily_before_during_preg2",
               "daily_before_during_preg3",
               "bmi_before","parity_b")
    
}else if(role=="father"){
    adj_vars=c("FARS_ALDER",
           "bmi_father","pat_spouse_smoke_preg")
}

cat("adjusting variables:",adj_vars,"\n")


work_horse=function(i){
  dat=na.omit(data.frame(expr[,i],
                         info[,c(exposure,"c_AMP_Plate",adj_vars)]))
  colnames(dat)=c("y","x","plate",adj_vars)
  mytry=try(solve(t(as.matrix(dat$x))%*%as.matrix(dat$x)),silent = T)
  if(class(mytry)[1]=="try-error"){
    output=c(i,rep(NA,3+length(names(table(dat$x)))))
  }else{
    model=Rfast::rint.reg(y=dat$y,x=as.matrix(dat[,c("x",adj_vars)]),id=as.factor(dat$plate))
    output=c(i,cbind(model$be,model$seb)[2,],
             nrow(dat)-(length(setdiff(colnames(dat),c("y","plate"))) + length(unique(dat$plate))),
             as.numeric(table(dat$x)))
  }
  return(output)        
}


work_horse(1)
ewas=parallel::mclapply(1:ncol(expr),work_horse,mc.cores=12)
print(length(ewas))
#ewas=data.frame(matrix(unlist(lapply(1:length(ewas), FUN=function(i){
#  if(length(ewas[[i]])!=4+length(names(table(info[,exposure])))){
#    output=c(ewas[[i]][1],rep(NA,3+length(names(table(info[,exposure])))))}else{output=ewas[[i]]}
#  return(output)
#})),ncol=4+length(names(table(info[,exposure]))),nrow=length(ewas),byrow = T))
ewas=data.frame(do.call("rbind",ewas));
print(dim(ewas));
ewas[,"Name"] = colnames(expr)[ewas[,1]]
ewas$t=ewas$X2/ewas$X3;ewas$p=pt(abs(ewas$t),df=ewas$X4,lower.tail = F)*2;ewas=ewas[order(ewas$p),]
ewas$X1=c();colnames(ewas)[1:5]=c("b","se","df","n0","n1")
ewas=ewas[,c("Name","b","se","df","n0","n1","t","p")]
ewas$q=p.adjust(ewas$p,"BH")
cat("EWAS result:\n");print(head(ewas,10))


feather::write_feather(ewas,path=paste("result_data/ewas_",role,"_",exposure,".feather",sep=""))