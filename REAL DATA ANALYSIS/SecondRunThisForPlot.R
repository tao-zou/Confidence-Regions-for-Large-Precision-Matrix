rm(list=ls())
library(lars)
library(MASS)
library(genlasso)
library(scalreg)
library(foreach)
library(doParallel)
library(igraph)

library(xtable)

library(GGally)
library(network)
library(sna)
library(ggplot2)

#registerDoParallel(cores=50)
#setwd('~/Desktop/casestudy1/Final')

#year




year=4

yearly=rbind(
  c(1,252),
  c(253,503),
  c(504,754),
  c(755,1007),
  c(1008,1259)
)

#--- Read data ---


zindex=read.csv("closeindex.csv")
A1 = read.csv("newdataclose.csv", head = TRUE)
A1 = A1[, -1]
Sector = A1[1, ]
Category = A1[2, ]
A1 = A1[-c(1, 2), ]
A1 = A1[yearly[year,1]:yearly[year,2],]


M = dim(A1)




categ = c(); sec = c()
for (j in 1 : M[2]){
  categ[j] = Category[1, j]
  sec[j] = Sector[1, j]
}
categ[131] = 59
sec[sec == 11] = 10

categ0 = categ[-which(categ == 99)]
sec0 = sec[-which(categ == 99)]
for (i in 1 : 10){
  cat(c(length(unique(categ0[sec0 == i])), unique(categ0[sec0 == i])), "\n")
}
IndSector1 = c(1 : 15)
IndSector2 = c(16 : 22)
IndSector3 = c(23 : (26-1))
IndSector4 = c(27 : 31)-1
IndSector5 = c(32 : 34)-1
IndSector6 = c(35 : 41)-1
IndSector7 = c(42 : 47)-1
IndSector8 = c(48 : 52)-1
IndSector9 = 53-1
IndSector10 = c(54, 55)-1

#IndSector = c(rep(1, 15), 2, 2, 3, rep(2, 5), rep(3, 3), rep(4, 5), rep(5, 3), rep(6, 8), rep(7, 5), rep(8, 5), 9, 10, 10)
IndSector = list(IndSector1, IndSector2, IndSector3, IndSector4, IndSector5, IndSector6, IndSector7, IndSector8, IndSector9, IndSector10)

secc=cbind(c(1 : 55), sort(unique(categ0)))







A1=A1[,-zindex[,2]]

M = dim(A1)
  
  

M

sp500company=read.csv("sp500company.csv", head = TRUE)


head(sp500company)

sort(unique(sp500company$Category...Note))

#sp500company$Category=99

#for (i in 1:length(sp500company[,1])){
#  if (sp500company$Category...Note[i] != 99){
#  sp500company$Category[i]=secc[which(secc[,2]==sp500company$Category...Note[i]),1]
#}
#}


sp500company=sp500company[-zindex[,2],]

G=length(sort(unique(sp500company$Category...Note)))-1


company=foreach (i = 1:G,.combine='rbind')%do%{
  stock=''
  for (j in 1:length(sp500company[,1])){
    if(sp500company$Category...Note[j]==(sort(unique(sp500company$Category...Note)))[i]){
      stock=paste0(stock,', ',as.character(sp500company$Listed.ticketer[j]))
      Industry=as.character(sp500company$GICS.Sub.Industry[j])
      IndustryNo=sp500company$Category...Note[j]
      Sector=as.character(sp500company$GICS.Sector[j])
      SectorNo=sp500company$Sectcor[j]
    }
  }
  cbind(data.frame(rbind(c(stock,Industry,Sector))),data.frame(rbind(c(IndustryNo,SectorNo))))
}

write.csv(company,file='company1.csv')


companyf=read.csv('companyf.csv')

companyf=companyf[,c(1,3,5,2,4)]

colnames(companyf)=c('Stock Symbols','Sectors', 'Sector No.', 'Sub Industries', 'Industry No.')

companyf=data.frame(companyf,row.names=NULL)

x=xtable(companyf)
print(x)

#sp500company=sp500company[,-4]


#-----Test for blocks, category-----

#OC = order(Category)
#SC = sort(Category)
categ=categ[-zindex[,2]]
G = length(unique(categ)) - 1


#-----Plot Graph of category-----

#PVBlock = PVBlock[, -G]
#PVBlock = PVBlock[-G, ]

filname=paste0("LogRCategoryPV_",year,".txt")

PVBlock=read.table(filname, header = FALSE)





PVBlocktemp = PVBlock
shuffle = c(1 : 17, 19 : 23, 18, 24 : 37,39:41,38,42:54)
PVBlocktemp1 = PVBlocktemp[, shuffle]
PVBlock = PVBlocktemp1[shuffle, ]

company=company[shuffle,]
write.csv(company,file='company.csv')

M = dim(PVBlock)
PV = c()
m = 1
for (i in 1 : (M[1] - 1)){
  for (j in (i + 1) : M[1]){
    PV[m] = PVBlock[i, j]
    m = m + 1
  }
}

#file=paste0('IndustryHist_',year,'.pdf')
#pdf(file,width=12,height=8)
hist(PV, xlab = "Pvalue", main = "Pvalue Histogram for Crossing Industry Blocks")
#dev.off()


alpha = 0.1
LenPV = length(PV)


judge=1-1*(sort(PV) < (alpha * c(1 : LenPV) / LenPV))


cri=1
for (i in length(judge):1){
  cri=cri*judge[i]
  if (cri==0)
    break
}

sort(PV)[i+1]


Connect = 1 * (PVBlock < sort(PV)[i+1])



M = dim(Connect)
edge = c()
for (i in 1 : (M[1] - 1)){
  for (j in (i + 1) : M[1]){
    if (Connect[i, j] == 1) edge = c(edge, c(i, j))
  }
}


head(Connect)

diag(Connect)=0
ConnectIndustry=matrix(0,nrow=10,ncol=10)
for (i in 1:length(IndSector)){
  for (j in 1:length(IndSector)){
    if (i == j){
    ConnectIndustry[i,j]=sum(Connect[IndSector[[i]],IndSector[[j]]])/2
    }else{
    ConnectIndustry[i,j]=sum(Connect[IndSector[[i]],IndSector[[j]]])
      }
  }
}

ConnectIndustry

result2=foreach(i = 1:length(IndSector),.combine='rbind')%do%{
  c(ConnectIndustry[i,i],sum(ConnectIndustry[i,-i]))
}

#result1

x=xtable(cbind(result2))
digits(x)<-0
print(x)





g = graph(edge, n = M[1], directed = FALSE)

file=paste0('Industry_',year,'.pdf')
pdf(file,width=11.5,height=10)
plot.igraph(g, layout = layout.circle, vertex.size = 10, vertex.label.cex = 2.0, vertex.label = c(1 : M[1]), 
            mark.groups = IndSector, mark.shape = 1, 
            #vertex.label = c("Consumer Discretionary", "Consumer Staples", "Energy", "Financials", "Health Care", "Industrials", "Information Technology", "Materials", "Telecom", "Utilities"), 
            vertex.color = NA, vertex.shape = "circle", vertex.frame.color = "black", vertex.label.color = "black", edge.width = 1.5, edge.color = "dimgrey", 
            main = "")
title("Partial Correlation Network of S&P 500 Sub Industries",font.main=1,cex.main=2.5)
text(1, 1, "Consumer Discretionary",cex=2)
text(-0.95, 1, "Consumer Staples",cex=2)
text(-1.2, 0.5, "Energy",cex=2)
text(-1.4, 0, "Financials",cex=2)
text(-1.3, -0.5, "Health Care",cex=2)
text(-0.8, -1, "Industrials",cex=2)
text(0.85, -1.1, "Information Technology",cex=2)
text(1.1, -0.7, "Materials",cex=2)
text(1.2, -0.35, "Telecom",cex=2)
text(1.25, -0.2, "Utilities",cex=2)
dev.off()














#--plotSector


sec[sec == 11] = 10
#OC = order(Sector)
#SC = sort(Sector)
sec=sec[-zindex[,2]]
G = length(unique(sec))


filname=paste0("LogRSectorPV_",year,".txt")

PVBlock=read.table(filname, header = FALSE)



M = dim(PVBlock)
PV = c()
m = 1
for (i in 1 : (M[1] - 1)){
  for (j in (i + 1) : M[1]){
    PV[m] = PVBlock[i, j]
    m = m + 1
  }
}

#file=paste0('SectorHist_',year,'.pdf')
#pdf(file,width=12,height=8)
hist(PV, xlab = "Pvalue", main = "Pvalue Histogram for Crossing Sector Blocks")
#dev.off()


alpha = 0.1
LenPV = length(PV)

judge=1-1*(sort(PV) < (alpha * c(1 : LenPV) / LenPV))


cri=1
for (i in length(judge):1){
  cri=cri*judge[i]
  if (cri==0)
    break
}

sort(PV)[i+1]


Connect = 1 * (PVBlock < sort(PV)[i+1])

Connect=1*(ConnectIndustry>0)

#--- Plot Network Graphs ---

M = dim(Connect)
edge = c()
for (i in 1 : (M[1] - 1)){
  for (j in (i + 1) : M[1]){
    if (Connect[i, j] == 1) edge = c(edge, c(i, j))
  }
}

g = graph(edge, n = 10, directed = FALSE)






##
mat=matrix(edge,ncol=2,byrow=TRUE)
mat=data.frame(mat)
colnames(mat)=c('Source','Target')


net = network(mat, directed = FALSE)

set.seed(7)

z=ggnet2(net, size='degree',max_size=30,
       label=c("Consumer Discretionary", "Consumer Staples", "Energy", 
               "Financials", "Health Care", "Industrials", "Information Technology", "Materials", "Telecom", "Utilities"),
       label.size=7,legend.position = "bottom",legend.size = 20)

coo=z$data

#write.csv(coo,file='position20171001.csv',row.names=FALSE)

coo=read.csv('position20171001.csv',header=T)

coo[9,6]=(coo[9,6]+coo[1,6])*2/3

coo[3,6]=coo[2,6]
coo[10,6]=coo[2,6]+0.05
coo[10,7]=(coo[3,7]+coo[8,6])/2

file=paste0('SectorA_',year,'.pdf')
pdf(file,width=14,height=10)
ggnet2(net, mode=as.matrix(coo[,c(6,7)]),size='degree',max_size=30,
       label=c("Consumer Discretionary", "Consumer Staples", "Energy", 
               "Financials", "Health Care", "Industrials", "Information Technology", "Materials", "Telecom", "Utilities"),color='lightgrey',alpha=0.75,edge.size=0.5,
       label.size=7,legend.position = "right",legend.size = 20)+ ggtitle('Partial Correlation Network of S&P 500 Sectors')+theme(plot.title=element_text(hjust=0.5, size=31))+theme(legend.background = element_rect(fill="white",
                                  size=0.5, linetype="solid", 
                                  colour ="black"))
#ggnet2(net,size='degree',max_size=30,
#       label=c("Consumer Discretionary", "Consumer Staples", "Energy", 
#               "Financials", "Health Care", "Industrials", "Information Technology", "Materials", "Telecom", "Utilities"),color='lightgrey',alpha=0.75,edge.size=0.5,
#       label.size=7,legend.position = "right",legend.size = 20)+ ggtitle('Partial Correlation Network of S&P 500 Sectors')+theme(plot.title=element_text(hjust=0.5, size=31))+theme(legend.background = element_rect(fill="white",
#                                  size=0.5, linetype="solid", 
#                                  colour ="black"))
dev.off()



#set.seed(5)
#x = gplot.layout.fruchtermanreingold(net, NULL)

#write.table(x,file='position.txt',col.names=FALSE,row.names=FALSE)


Connect

1*(ConnectIndustry>0)

Connect-1*(ConnectIndustry>0)

