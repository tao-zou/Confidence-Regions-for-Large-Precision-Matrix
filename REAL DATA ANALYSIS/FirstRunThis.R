rm(list=ls())
library(lars)
library(MASS)
library(genlasso)
library(scalreg)
library(foreach)
library(doParallel)
library(igraph)



###New
library(moments)
library(tseries)
library(fGarch)


registerDoParallel(cores=6)
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
IndSector6 = c(35 : 42)-1
IndSector7 = c(43 : 47)-1
IndSector8 = c(48 : 52)-1
IndSector9 = 53-1
IndSector10 = c(54, 55)-1

#IndSector = c(rep(1, 15), 2, 2, 3, rep(2, 5), rep(3, 3), rep(4, 5), rep(5, 3), rep(6, 8), rep(7, 5), rep(8, 5), 9, 10, 10)
IndSector = list(IndSector1, IndSector2, IndSector3, IndSector4, IndSector5, IndSector6, IndSector7, IndSector8, IndSector9, IndSector10)

cbind(c(1 : 55), sort(unique(categ0)))







A1=A1[,-zindex[,2]]

M = dim(A1)
  
  
####
R = matrix(0, (M[1] - 1), M[2])
P = matrix(0, (M[1] - 1), M[2])
####

##
for (j in 1 : M[2]){
  for (i in 1 : (M[1] - 1)){
    R[i, j] = log(A1[i + 1, j]) - log(A1[i, j])
    P[i, j] = 100 * (A1[i + 1, j] / A1[i, j] - 1)
  }
}

#RS = matrix(0, (M[1] - 1), M[2])
RS=foreach(j = 1 : M[2],.combine='cbind') %dopar% {
  RSj = ( R[, j] - mean(R[, j]) ) / sd(R[, j])
  RSj
}

###New
sum(apply(RS,2,kurtosis)>3)
RS=foreach(j = 1 : M[2],.combine='cbind') %dopar% {
  RSj = ( R[, j] - mean(R[, j]) ) 
  RSj
}
RSsd=foreach(j = 1 : M[2],.combine='cbind') %dopar% {
m1=garchFit(~ garch(1,1), data = RS[,j], trace = FALSE)
vectorm1=volatility(m1, type = c("sigma"))    
vectorm1
}
RS=RS/RSsd


head(RS)
X = RS
#setwd('~/Desktop/casestudy1/Final')
##
filname=paste0("RS_",year,".txt")

write.table(RS, filname, col.names = FALSE, row.names = FALSE)

n = dim(X)[1]; p = dim(X)[2]; smax = n / 2
Eresidual = matrix(0, n, p)
CoefMatrix = matrix(0, p, p - 1)

result=foreach (i = 1 : p) %dopar% {
  out = scalreg(X = X[, -i], y = X[, i], lam0 = sqrt(2 * log(p * log(p) / sqrt(n)) / n))
  Eresiduali = out$residuals
  CoefMatrixi = out$coefficients
  if ( sum(abs(CoefMatrixi) > 10^(-6)) > smax ){
    out = genlasso(X[, i], X = X[, -i], D = diag(1, p - 1))
    Coef = coef(out,  lambda = 2 * sqrt(var(X[, i]) * n * log(p)))
    Predict = predict(out, lambda = 2 * sqrt(var(X[, i]) * n * log(p)), Xnew = X[, -i])
    CoefMatrixi = t(Coef$beta)
    Eresiduali = X[, i] - Predict$fit
  }
  list(CoefMatrixi,Eresiduali)
}


for (i in 1:p){
  CoefMatrix[i, ]=result[[i]][[1]]
  Eresidual[, i]=result[[i]][[2]]
}

CovRes = t(Eresidual) %*% Eresidual / n

QS = function(u){
  if (u == 0) ker = 1
  else ker = 25 * ( sin(6 * pi * u / 5) / (6 * pi * u / 5) - cos(6 * pi * u / 5) ) / (12 * pi^2 * u^2)
  return(ker)
}

#--- Overall maximum ---

TestAll = matrix(0, n, p * (p - 1) / 2)
BTAll = matrix(0, n, p * (p - 1) / 2)
m = 1
for (i in 1 : (p - 1)){
  for (j in (i + 1) : p){
    TestAll[, m] = - ( Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1] ) / (diag(CovRes)[i] * diag(CovRes)[j])
    BTAll[, m] = - ( Eresidual[, i] * Eresidual[, j] ) / (diag(CovRes)[i] * diag(CovRes)[j]) + mean(TestAll[, m])
    m = m + 1
  }
}

BTAllMean = colMeans(BTAll)
BTAll = t(t(BTAll) - BTAllMean)

NumAll = 0
DenAll = 0
for (i in 1 : (p * (p - 1) / 2))  {
  AR1 = ar(BTAll[, i], aic = FALSE, order.max = 1)
  rhoEst = AR1$ar
  sigma2Est = AR1$var.pred
  NumAll = NumAll + 4 * (rhoEst * sigma2Est)^2 / (1 - rhoEst)^8
  DenAll = DenAll + sigma2Est^2 / (1 - rhoEst)^4
}

a2All = NumAll / DenAll
bandwidthAll = 1.3221 * (a2All * n)^(0.2)
BTcovAll = matrix(0, n, n)
for (i in 1 : n){
  for (j in 1 : n){
    BTcovAll[i, j] = QS(abs(i - j) / bandwidthAll)
  }
}

WdiagAllEmp = 0
for (i in 1 : n){
  WdiagAllEmp = WdiagAllEmp + BTAll[i, ]^2
}
WdiagAllEmp = WdiagAllEmp / n

M0 = 3000
#BTAllsim = rep(0,M0)

BTAllsim=foreach (i = 1 : M0,.combine='c') %dopar% {
  set.seed(i)
  temp = mvrnorm(1, rep(0, n), BTcovAll)
  BTAllsimi = (n)^(-0.5) * max(WdiagAllEmp^(-1/2) * abs(colSums(temp * BTAll)))
  BTAllsimi
}

QAll = sort(BTAllsim)[0.9 * M0]

Estimate = colMeans(TestAll)
uplimit = Estimate + WdiagAllEmp^(1/2) * QAll / sqrt(n)
lowlimit = Estimate - WdiagAllEmp^(1/2) * QAll / sqrt(n)

sum(uplimit * lowlimit > 0)



#-----Test for blocks, category-----

#OC = order(Category)
#SC = sort(Category)
categ=categ[-zindex[,2]]
G = length(unique(categ)) - 1
PVBlock = matrix(0, G, G)
StatBlock = matrix(0, G, G)

seedindex=0
for (Ci in 1 : (G - 1)){
  Bi = sort(unique(categ))[Ci]
  IndexI = which(categ == Bi)
  for (Cj in (Ci + 1) : G){
    seedindex=seedindex+1
    Bj = sort(unique(categ))[Cj]
    IndexJ = which(categ == Bj)
    BL = length(IndexI) * length(IndexJ)
    TestAll = matrix(0, n, BL)
    BTAll = matrix(0, n, BL)
    m = 1
    for (Ii in 1 : length(IndexI)){
      i0 = IndexI[Ii]
      for (Ij in 1 : length(IndexJ)){
        j0 = IndexJ[Ij]
        if (j0 < i0){
          i = j0
          j = i0
        }
        else {
          i = i0
          j = j0
        }
        TestAll[, m] = - ( Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1] ) / (diag(CovRes)[i] * diag(CovRes)[j])
        BTAll[, m] = - ( Eresidual[, i] * Eresidual[, j] ) / (diag(CovRes)[i] * diag(CovRes)[j]) + mean(TestAll[, m])
        m = m + 1
        #cat("i, j = ", c(Ii, Ij, i0, j0, i, j), "\n")
      }
    }
    BTAllMean = colMeans(BTAll)
    BTAll = t(t(BTAll) - BTAllMean)
    
    NumAll = 0
    DenAll = 0
    for(i in 1 : BL){
      AR1 = ar(BTAll[, i], aic = FALSE, order.max = 1)
      rhoEst = AR1$ar
      sigma2Est = AR1$var.pred
      NumAll = NumAll + 4 * (rhoEst * sigma2Est)^2 / (1 - rhoEst)^8
      DenAll = DenAll + sigma2Est^2 / (1 - rhoEst)^4
    }
    a2All = NumAll / DenAll
    bandwidthAll = 1.3221 * (a2All * n)^(0.2)
    BTcovAll = matrix(0, n, n)
    for (i in 1 : n){
      for (j in 1 : n){
        BTcovAll[i, j] = QS(abs(i - j) / bandwidthAll)
      }
    }
    
    WdiagAllEmp = 0
    for (i in 1 : n){
      WdiagAllEmp = WdiagAllEmp + BTAll[i, ]^2
    }
    WdiagAllEmp = WdiagAllEmp / n
    TestStat = sqrt(n) * max(WdiagAllEmp^(-1/2) * abs(colMeans(TestAll)))
    
    M0 = 1000
    #BTAllsim = rep(0,M0)
    BTAllsim=foreach (i = 1 : M0, .combine='c') %dopar% {
      set.seed(i+10000*seedindex)
      temp = mvrnorm(1, rep(0, n), BTcovAll)
      BTAllsimi = (n)^(-0.5) * max(WdiagAllEmp^(-1/2) * abs(colSums(temp * BTAll)))
      BTAllsimi
    }
    QAll = sort(BTAllsim)[0.95 * M0]
    PVBlock[Ci, Cj] = mean(BTAllsim > TestStat)
    PVBlock[Cj, Ci] = PVBlock[Ci, Cj]
    StatBlock[Ci, Cj] = TestStat
    StatBlock[Cj, Ci] = StatBlock[Ci, Cj]
    
    cat("Blocks = ", c(Bi, Bj, 1 * (TestStat > QAll)), "\n")
  }
}

#-----Plot Graph of category-----

#PVBlock = PVBlock[, -G]
#PVBlock = PVBlock[-G, ]

filname=paste0("LogRCategoryPV_",year,".txt")

write.table(PVBlock, filname, col.names = FALSE, row.names = FALSE)

filname=paste0("LogRCategoryStat_",year,".txt")

write.table(StatBlock, filname, col.names = FALSE, row.names = FALSE)




PVBlocktemp = PVBlock
shuffle = c(1 : 17, 19 : 23, 18, 24 : 54)
PVBlocktemp1 = PVBlocktemp[, shuffle]
PVBlock = PVBlocktemp1[shuffle, ]


M = dim(PVBlock)
PV = c()
m = 1
for (i in 1 : (M[1] - 1)){
  for (j in (i + 1) : M[1]){
    PV[m] = PVBlock[i, j]
    m = m + 1
  }
}

file=paste0('IndustryHist_',year,'.pdf')
pdf(file,width=12,height=8)
hist(PV, xlab = "Pvalue", main = "Pvalue Histogram for Crossing Industry Blocks")
dev.off()


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

g = graph(edge, n = M[1], directed = FALSE)

file=paste0('Industry_',year,'.pdf')
pdf(file,width=12,height=8)
plot.igraph(g, layout = layout.circle, vertex.size = 10, vertex.label.cex = 0.6, vertex.label = c(1 : M[1]), 
            mark.groups = IndSector, mark.shape = 1, 
            #vertex.label = c("Consumer Discretionary", "Consumer Staples", "Energy", "Financials", "Health Care", "Industrials", "Information Technology", "Materials", "Telecom", "Utilities"), 
            vertex.color = NA, vertex.shape = "circle", vertex.frame.color = "blue", vertex.label.color = "black", edge.width = 1.5, edge.color = "red", 
            main = "Conditional Dependence Graph of S&P500 Industries")

text(1, 1, "Consumer Discretionary")
text(-0.95, 1, "Consumer Staples")
text(-1.2, 0.5, "Energy")
text(-1.4, 0, "Financials")
text(-1.3, -0.5, "Health Care")
text(-0.8, -1, "Industrials")
text(0.85, -1.1, "Information Technology")
text(1.1, -0.7, "Materials")
text(1.2, -0.35, "Telecom")
text(1.25, -0.2, "Utilities")
dev.off()














#--plotSector


sec[sec == 11] = 10
#OC = order(Sector)
#SC = sort(Sector)
sec=sec[-zindex[,2]]
G = length(unique(sec))
PVBlock = matrix(0, G, G)
StatBlock = matrix(0, G, G)

for (Ci in 1 : (G - 1)){
  Bi = sort(unique(sec))[Ci]
  IndexI = which(sec == Bi)
  for (Cj in (Ci + 1) : G){
    seedindex=seedindex+1
    Bj = sort(unique(sec))[Cj]
    IndexJ = which(sec == Bj)
    BL = length(IndexI) * length(IndexJ)
    TestAll = matrix(0, n, BL)
    BTAll = matrix(0, n, BL)
    m = 1
    for (Ii in 1 : length(IndexI)){
      i0 = IndexI[Ii]
      for (Ij in 1 : length(IndexJ)){
        j0 = IndexJ[Ij]
        if (j0 < i0){
          i = j0
          j = i0
        }
        else {
          i = i0
          j = j0
        }
        TestAll[, m] = - ( Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1] ) / (diag(CovRes)[i] * diag(CovRes)[j])
        BTAll[, m] = - ( Eresidual[, i] * Eresidual[, j] ) / (diag(CovRes)[i] * diag(CovRes)[j]) + mean(TestAll[, m])
        m = m + 1
        #cat("i, j = ", c(Ii, Ij, i0, j0, i, j), "\n")
      }
    }
    BTAllMean = colMeans(BTAll)
    BTAll = t(t(BTAll) - BTAllMean)
    
    NumAll = 0
    DenAll = 0
    for(i in 1 : BL){
      AR1 = ar(BTAll[, i], aic = FALSE, order.max = 1)
      rhoEst = AR1$ar
      sigma2Est = AR1$var.pred
      NumAll = NumAll + 4 * (rhoEst * sigma2Est)^2 / (1 - rhoEst)^8
      DenAll = DenAll + sigma2Est^2 / (1 - rhoEst)^4
    }
    a2All = NumAll / DenAll
    bandwidthAll = 1.3221 * (a2All * n)^(0.2)
    BTcovAll = matrix(0, n, n)
    for (i in 1 : n){
      for (j in 1 : n){
        BTcovAll[i, j] = QS(abs(i - j) / bandwidthAll)
      }
    }
    
    WdiagAllEmp = 0
    for (i in 1 : n){
      WdiagAllEmp = WdiagAllEmp + BTAll[i, ]^2
    }
    WdiagAllEmp = WdiagAllEmp / n
    TestStat = sqrt(n) * max(WdiagAllEmp^(-1/2) * abs(colMeans(TestAll)))
    
    M0 = 10000
    #BTAllsim = rep(0,M0)
    BTAllsim=foreach (i = 1 : M0 , .combine='c') %dopar% {
      set.seed(i+10000*seedindex)
      temp = mvrnorm(1, rep(0, n), BTcovAll)
      BTAllsimi = (n)^(-0.5) * max(WdiagAllEmp^(-1/2) * abs(colSums(temp * BTAll)))
      BTAllsimi
      }
    QAll = sort(BTAllsim)[0.95 * M0]
    PVBlock[Ci, Cj] = mean(BTAllsim > TestStat)
    PVBlock[Cj, Ci] = PVBlock[Ci, Cj]
    StatBlock[Ci, Cj] = TestStat
    StatBlock[Cj, Ci] = StatBlock[Ci, Cj]
    
    cat("Sectors = ", c(Bi, Bj, 1 * (TestStat > QAll)), "\n")
  }
}

filname=paste0("LogRSectorPV_",year,".txt")

write.table(PVBlock, filname, col.names = FALSE, row.names = FALSE)

filname=paste0("LogRSectorStat_",year,".txt")

write.table(StatBlock, filname, col.names = FALSE, row.names = FALSE)

M = dim(PVBlock)
PV = c()
m = 1
for (i in 1 : (M[1] - 1)){
  for (j in (i + 1) : M[1]){
    PV[m] = PVBlock[i, j]
    m = m + 1
  }
}

file=paste0('SectorHist_',year,'.pdf')
pdf(file,width=12,height=8)
hist(PV, xlab = "Pvalue", main = "Pvalue Histogram for Crossing Sector Blocks")
dev.off()


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

#--- Plot Network Graphs ---

M = dim(Connect)
edge = c()
for (i in 1 : (M[1] - 1)){
  for (j in (i + 1) : M[1]){
    if (Connect[i, j] == 1) edge = c(edge, c(i, j))
  }
}

g = graph(edge, n = 10, directed = FALSE)

file=paste0('Sector_',year,'.pdf')
pdf(file,width=12,height=8)
plot.igraph(g, layout = layout_with_fr, vertex.label = c("Consumer Discretionary", "Consumer Staples", "Energy", 
                                                         "Financials", "Health Care", "Industrials", "Information Technology", "Materials", "Telecom", "Utilities"), vertex.color = NA, 
            vertex.shape = "circle", vertex.frame.color = "blue", vertex.label.color = "black", edge.width = 1.5, edge.color = "red", 
            main = "Conditional Dependence Graph of S&P500 Sectors")
dev.off()

