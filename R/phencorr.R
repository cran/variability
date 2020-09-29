#' @title Phenotypic Correlation Analysis
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes
#' @param replication vector containing replications
#' @return Phenotypic correlation
#' @export
#' @importFrom stats anova lm pt
#' @examples
#' data(vardata)
#' pheno.corr(vardata[3:11],vardata$Genotypes,vardata$Replication)
pheno.corr<-function(data,genotypes,replication){
  convert<-function (data1) {
    data1 <- as.data.frame(sapply(data1, as.numeric))
    data1 <- as.list(data1)
    return(data1)
  }
  genotypes<-as.factor(genotypes)
  replication<-as.factor(replication)
  colnumber <- ncol(data)
  headings<-names(data)
  data2<-convert(data)
  phenotypic.cor1<-function(genotypes,replication,trait1,trait2){
    genotypes<-as.factor(genotypes)
    replication<-as.factor(replication)
    sumch1<-tapply(trait1,genotypes,sum)
    sumch2<-tapply(trait2,genotypes,sum)
    sumr1<-tapply(trait1,replication,sum)
    sumr2<-tapply(trait2,replication,sum)
    repli<-nlevels(replication)
    genotype<-nlevels(genotypes)
    tdf<-repli*genotype
    GT1<-sum(trait1)
    GT2<-sum(trait2)
    CF<-(GT1*GT2)/(repli*genotype)
    TSP<-round(sum(trait1*trait2)-CF,6)
    GSP<-round((sum(sumch1*sumch2)/repli)-CF,6)
    RSP<-round((sum(sumr1*sumr2)/genotype)-CF,6)
    ESP<-TSP-GSP-RSP
    DFR<-repli-1
    DFG<-genotype-1
    DFE<-DFR*DFG
    RMP<-round(RSP/DFR,6)
    GMP<-round(GSP/DFG,6)
    EMP<-round(ESP/DFE,6)
    EnvCov<-EMP
    GenCov<-round((GMP-EMP)/repli,6)
    PhenoCov<-EnvCov+GenCov
    model <- lm(trait1 ~ replication + genotypes)
    anova.model <- anova(model)
    EMS<-anova.model[3,3]
    GV<-round((anova.model[2,3]-EMS)/repli,6)
    PV<-GV+EMS
    model1 <- lm(trait2 ~ replication + genotypes)
    anova.model1 <- anova(model1)
    EMS1<-anova.model1[3,3]
    GV1<-round((anova.model1[2,3]-EMS1)/repli,6)
    PV1<-EMS1+GV1
    r<-round(PhenoCov/sqrt(PV*PV1),4)
    SEm<-sqrt((abs(1-(r^2)))/(tdf-2))
    tvalue<-r/SEm
    pvalue<-2*pt(abs(tvalue),tdf-2,lower.tail = F)
    ifelse(pvalue<=0.01,r<-paste(r,"**"),ifelse(pvalue>0.05,r<-paste(r,"NS"),r<-paste(r,"*")))
    return(r)
  }
  phenotypic.corr <- c()
  index=0
  for (i in 1:(colnumber)){
    for (j in 1:colnumber){
      index=index+1
      phenotypic.corr[index]<-phenotypic.cor1(genotypes,replication,data2[[i]],data2[[j]])
    }}
  matribapu<-noquote(matrix(phenotypic.corr,nrow = colnumber,dimnames =list(headings,headings)))
  matribapu1<-as.table(matribapu,useNa=F)
  Note<-"The sig of phenotypic correlation was tested using t test (two-tail). The degree of freedom used is (genotypes*replication) - 2"
  output<-list(PhenotypicCorrelation=matribapu1,Note=Note)
  return(output)
}
