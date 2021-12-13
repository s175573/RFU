# load('trimerEncodingWorkplace0901.Rdata')
# training data:
#     lung cancer: Creelan et al, blood samples, n=20
#                  Cascone et al., blood samples, n=15
#     control: HC blood samples from Carter Blood (Li lab inhouse), n=34
#              Adaptive covid INCOV cohort, n=157

# testing data:
#     lung cancer: MD anderson 2020 cohort, n=121
#                  Federic et al., blood??, n=25
#     control: Adaptive covid adaptive cohort, n=160
#              Adaptive Control 2020 cohort, n=88
## Lung cancer analysis
lung.d2=cbind(RFU.lung.Creelan, RFU.control.inhouse)  ## control.inhouse: Carter blood 

tmp.col.ld2=c(rep(1,20), rep(2,34))

lung.d3=cbind(RFU.lung.cascone[,grep('blood', colnames(RFU.lung.cascone))], RFU.covid.INCOV)  ## control.HC: Adaptive control 2020
tmp.col.ld3=c(rep(1,15),rep(2, 157))

lung.test = cbind(RFU.lung.MDA, RFU.covid.adaptive)
tmp.col.test=c(rep(1,121), rep(2, 160))

lung.test2 = cbind(RFU.lung.Federic.Nf, RFU.control.HC)
tmp.col.test2=c(rep(1,25),rep(2, 88))

#lung.test3 = cbind(RFU.lung.JHU.b, RFU.batch2) ## JHU 2020 cohort; control: Emerson 2017 batch2
#tmp.col.test3=c(rep(1,27), rep(2,120))

# test for batch effect
tmp.dd=lung.test2
tmp.lab = tmp.col.test2
tmp.pca=prcomp(t(tmp.dd[rr.hk,]))$x
par(mfrow=c(2,2), mar=c(4,4,2,0))
for(ii in 1:4)boxplot(tmp.pca[,ii] ~ tmp.lab, main=paste('PC',ii,sep=''), ylab='PC score', xlab='')  ## make sure no significance for the top 4 PCs


  ## training accuracy
tmp.cor=cor(lung.d2, lung.d3)
knn.pred(sqrt(1-tmp.cor), tmp.col.ld2, n=3) -> cMat
tmp.rr=cMat[,1]/cMat[,2]
tmp.yy = 2-tmp.col.ld3
plot(roc(tmp.yy ~ tmp.rr)) ## AUC = 0.905

  ## independent test 1
tmp.cor=cor(cbind(lung.d2, lung.d3), lung.test)
tmp.cols=c(tmp.col.ld2, tmp.col.ld3)
tmp.cols=2-tmp.cols  ## make control = 0, cancer = 1

knn.pred(sqrt(1-tmp.cor), tmp.cols, n=3) -> cMat
tmp.rr=cMat[,1]/cMat[,2]
tmp.yy=2-tmp.col.test
plot(roc(tmp.yy ~ tmp.rr))  ## AUC=0.87

  ## independent test 2
tmp.cor=cor(cbind(lung.d2, lung.d3), lung.test2)
tmp.cols=c(tmp.col.ld2, tmp.col.ld3)
tmp.cols=2-tmp.cols  ## make control = 0, cancer = 1

knn.pred(sqrt(1-tmp.cor), tmp.cols, n=3) -> cMat
tmp.rr=cMat[,1]/cMat[,2]
tmp.yy=2-tmp.col.test2
plot(roc(tmp.yy ~ tmp.rr))  ## AUC=0.97
