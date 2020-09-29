#' @title Genotypic Path Analysis
#' @param dependent.var trait to be used a dependent variable
#' @param independent.var traits to be used as an indpendent variables
#' @param genotypes vector containing genotpes
#' @param replication vector containing replications
#' @return Direct effects, indirect effects and residual
#' @importFrom stats anova lm
#' @examples
#' data(vardata)
#' # Grain yield is considered as a dependent variable
#' geno.path(vardata[11],vardata[3:10],vardata$Genotypes,vardata$Replication)
#' @export
geno.path<-function(dependent.var,independent.var,genotypes,replication){
  convert<-function (data1) {
    data1 <- as.data.frame(sapply(data1, as.numeric))
    data1 <- as.list(data1)
    return(data1)
  }
  genotypes<-as.factor(genotypes)
  replication<-as.factor(replication)
  colnumber <- ncol(independent.var)
  totalnumber<-colnumber+1
  headings<-names(independent.var)
  totaldata<-data.frame(independent.var,dependent.var)
  data2<-convert(totaldata)
  gen.cor1<-function(genotypes,replication,trait1,trait2){
    sumch1<-tapply(trait1,genotypes,sum)
    sumch2<-tapply(trait2,genotypes,sum)
    sumr1<-tapply(trait1,replication,sum)
    sumr2<-tapply(trait2,replication,sum)
    repli<-nlevels(replication)
    genotype<-nlevels(genotypes)
    GT1<-sum(trait1)
    GT2<-sum(trait2)
    CF<-(GT1*GT2)/(repli*genotype)
    TSP<-round(sum(trait1*trait2)-CF,5)
    GSP<-round((sum(sumch1*sumch2)/repli)-CF,5)
    RSP<-round((sum(sumr1*sumr2)/genotype)-CF,5)
    ESP<-TSP-GSP-RSP
    DFR<-repli-1
    DFG<-genotype-1
    DFE<-DFR*DFG
    RMP<-round(RSP/DFR,5)
    GMP<-round(GSP/DFG,5)
    EMP<-round(ESP/DFE,5)
    EnvCov<-EMP
    GenCov<-round((GMP-EMP)/repli,5)
    model <- lm(trait1 ~ replication + genotypes)
    anova.model <- anova(model)
    EMS<-anova.model[3,3]
    GV<-abs(round((anova.model[2,3]-EMS)/repli,5))
    model1 <- lm(trait2 ~ replication + genotypes)
    anova.model1 <- anova(model1)
    EMS1<-anova.model1[3,3]
    GV1<-abs(round((anova.model1[2,3]-EMS1)/repli,5))
    r<-round(GenCov/sqrt(GV*GV1),5)
    return(r)
  }
  genetic.corr <- c()
  index=0
  for (i in 1:(totalnumber)){
    for (j in 1:totalnumber){
      index=index+1
      genetic.corr[index]<-gen.cor1(genotypes,replication,data2[[i]],data2[[j]])
    }}
  matribapu<-matrix(genetic.corr,nrow = totalnumber)
  corr.ind<-matribapu[1:colnumber,1:colnumber]
  corr.dep<-matribapu[1:colnumber,totalnumber]
  Direct <- solve(corr.ind, corr.dep)
  Coefficient <- corr.ind
  for (i in 1:colnumber) {
    for (j in 1:colnumber) {
      Coefficient[i, j] <- Direct[j] * corr.ind[i, j]
    }
  }
  Coefficient<-round(Coefficient,5)
  rownames(Coefficient)<-headings
  colnames(Coefficient)<-headings
  residual <- round(1 - t(Direct) %*% corr.dep,4)
  finaloutput<-list(Effects=Coefficient,Residual=residual)
  return(finaloutput)
}
