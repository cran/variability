#' @title Analysis of Covariance
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes
#' @param replication vector containing replications
#' @return ANCOVA, genotypic and phenotypic covariance
#' @importFrom stats pf
#' @export
#' @examples
#' data(vardata)
#' ancova(vardata[3:11],vardata$Genotypes,vardata$Replication)
ancova<-function(data,genotypes,replication){
  convert<-function (data1){
    data1 <- as.data.frame(sapply(data1, as.numeric))
    data1 <- as.list(data1)
    return(data1)}
  datam<-convert(data)
  num<-ncol(data)
  analysis<-function(genotypes,replication,trait1,trait2){
    genotypes<-as.factor(genotypes)
    replication<-as.factor(replication)
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
    GSP<-round((sum(sumch1*sumch2)/repli)-CF,4)
    RSP<-round((sum(sumr1*sumr2)/genotype)-CF,4)
    ESP<-TSP-GSP-RSP
    DFR<-repli-1
    DFG<-genotype-1
    DFE<-DFR*DFG
    RMP<-round(RSP/DFR,4)
    GMP<-round(GSP/DFG,4)
    EMP<-round(ESP/DFE,4)
    RFval<-round(RMP/EMP,4)
    GFval<-round(GMP/EMP,4)
    rpvalue<-round(pf(RFval,DFR,DFE,lower.tail = FALSE),4)
    gpvalue<-round(pf(GFval,DFG,DFE,lower.tail = FALSE),4)
    EnvCov<-EMP
    GenCov<-round((GMP-EMP)/repli,4)
    PhenCov<-round(EMP+((GMP-EMP)/repli),4)
    ANCOVA<-matrix(data=c(DFR,RSP,RMP,RFval,rpvalue,DFG,GSP,GMP,GFval,gpvalue,DFE,ESP,EMP,NA,NA,EMP,NA,NA,NA,NA,GenCov,NA,NA,NA,NA,PhenCov,NA,NA,NA,NA),dimnames = list(c("Replication", "Genotypes","Error","Environmental Covariance","Genotypic Covariance","Phenotypic Covariance"),c("DF", "SP","MP","F Cal","p-value")),nrow = 6,byrow = T)
    table1<-as.table(ANCOVA,useNa=F)
    return(table1)
  }
  list1 <- list()
  index=0
  for (i in 1:(num-1)){
    for (j in (i+1):num){
      index=index+1
      list1[[index]]<-analysis(genotypes,replication,datam[[i]],datam[[j]])
    }}
  naming<-names(datam)
  combi<-c()
  index=0
  for (i in 1:(num-1)){
    for (j in (i+1):num){
      index=index+1
      combi[index]<-paste(naming[i],naming[j])
    }}
  names(list1)<-combi
  return(list1)
}
