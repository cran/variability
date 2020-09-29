#' @title Estimation of Genetic Parameters
#' @param data traits to be analyzed
#' @param genotypevector vector containing genotypes
#' @param replicationvector vector containig replications
#' @importFrom stats anova lm qt
#' @return ANOVA, genotypic and phenotypic coefficient of variance, heritability, genetic advance and genetic advance as percentage of mean.
#' @export
#' @examples
#' data(vardata)
#' gen.var(vardata[3:11],vardata$Genotypes,vardata$Replication)
gen.var<-function (data, genotypevector, replicationvector)
{
  convert<-function (data1) {
    data1 <- as.data.frame(sapply(data1, as.numeric))
    data1 <- as.list(data1)
    return(data1)
  }
  analysis <- function(data1, genotypevector, replicationvector){
    data2<-as.numeric(data1)
    genotype <- as.factor(genotypevector)
    replication <- as.factor(replicationvector)
    r <- nlevels(replication)
    model <- lm(data2 ~ replication + genotype)
    anova.model <- anova(model)
    Maxi<-round(max(data2),4)
    Mini<-round(min(data2),4)
    GM<-round(mean(data2),4)
    EMS<-anova.model[3,3]
    Env.var<-round(EMS,4)
    SEm<-round(sqrt(EMS/r),4)
    CD5<-round(sqrt(EMS/r)*sqrt(2)*abs(qt(0.025,anova.model[3,1])),4)
    if(anova.model[2,5]>0.05){
      CD5<-paste(round(sqrt(EMS/r)*sqrt(2)*abs(qt(0.025,anova.model[3,1])),4),"NS")
    }
    CD1<-round(sqrt(EMS/r)*sqrt(2)*abs(qt(0.005,anova.model[3,1])),4)
    if(anova.model[2,5]>0.01){
      CD1<-paste(round(sqrt(EMS/r)*sqrt(2)*abs(qt(0.005,anova.model[3,1])),4),"NS")
    }
    GV<-round((anova.model[2,3]-EMS)/r,4)
    GVV<-round((anova.model[2,3]-EMS)/r,4)
    if(GV<0){
      GV<-paste("Note: GV is negative",round((anova.model[2,3]-EMS)/r,4))
    }
    PV<-round(GVV+EMS,4)
    if(GVV<0){
      GCV<-paste("Note: GV is negative,GCV calculated by using absolute GV",round((sqrt(abs(GVV))/mean(data2))*100,4))
    }else{
      GCV<-round((sqrt(GVV)/mean(data2))*100,4)
    }
    PCV<-round((sqrt(PV)/mean(data2))*100,4)
    ECV<-round((sqrt(EMS)/mean(data2))*100,4)
    hs<-round((GVV/PV),4)
    ga<-round((GVV/PV)*2.06*sqrt(PV),4)
    gam<-round((ga/mean(data2))*100,4)
    matri<-matrix(data = c(Maxi,Mini,GM,SEm,CD5,CD1,Env.var,GV,PV,ECV,GCV,PCV,hs,ga,gam),dimnames = list(c("Maximum","Minimum","Grand Mean","Standard Error of Mean (SEm)","Critical Difference (CD) 5%","Critical Difference (CD) 1%","Environmental Variance","Genotypic Variance","Phenotypic Variance","Environmental Coefficient of Variance","Genotypic Coefficient of Variance","Phenotypic Coefficient of Variance","Heritability (Broad Sense)","Genetic Advance","Genetic Advance as percentage of mean")),nrow = 15)
    table1<-as.table(matri,useNa=F)
    my.list<-list(anova.model,table1)
    return(my.list)
  }
  fiftn <- convert(data)
  colnumber <- ncol(data)
  output <- list()
  for (j in 1:colnumber) {
    output[[j]] <- analysis(fiftn[[j]], genotypevector, replicationvector)
  }
  names(output) <- names(data)
  return(output)
}
