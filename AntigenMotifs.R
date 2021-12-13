## Antigen motif encoding data analysis
library(Rtsne)
DIR='tmp/'
ffs=dir(DIR, full.names=T)
YYs=c()

for(ff in ffs){
  print(ff)
  dd=read.table(ff,header=F,sep='\t',stringsAsFactors = F)
  dd=unique(dd)
  TMP=Rtsne(dd[,5:20], check_duplicates=FALSE)
  YYs=c(YYs, list(TMP$Y))
}

names(YYs)=ffs


## tri-mer matrix
trimers=c()
AAstring='ACDEFGHIKLMNPQRSTVWY'
AAstringList=unlist(strsplit(AAstring,''))
for(A1 in AAstringList){
  for(A2 in AAstringList){
    for(A3 in AAstringList){
      trimers=c(trimers, paste(A1, A2, A3, sep=''))
    }
  }
}

trimerMat=matrix(0, 8000, 8000)
rownames(trimerMat)=colnames(trimerMat)=trimers

dd=read.table('TrainData_AntigenMotif_AntigenMotif2threePosition_Beta.py_1000_3_928990.txt',header=F,sep='\t', stringsAsFactors = F)
dd0=dd[which(dd[,3]==0),]
for(jj in 1:nrow(dd0)){
  s1=dd0[jj,1]
  s2=dd0[jj,2]
  s1=unlist(strsplit(s1, ':'))
  s2=unlist(strsplit(s2, ':'))
  s1=s1[[1]]
  s2=s2[[1]]
  n1=nchar(s1)
  s1.tri=sapply(1:(n1-2), function(x)substr(s1, x, x+2))
  n2=nchar(s2)
  s2.tri=sapply(1:(n2-2), function(x)substr(s2, x, x+2))
  trimerMat[s1.tri, s2.tri]=trimerMat[s1.tri, s2.tri]+1
}

## Using GIANA clustered big TCR file to build tri-mer replacement matrix
## Split trimers based on the positions of the first letter in the CDR3 sequence

# use the biggest clustered data, with 14 million sequences
dd=read.table('HC600_COVID600_AllCancer_MS_ARD--RotationEncodingBL62.txt',header=F,sep='\t',stringsAsFactors=F)

## select small-sized clusters, to exclude small-world effect
tt=table(dd[,2])
nn=names(which(tt<=6))
dd0=dd[which(dd[,2] %in% nn),]

## remove clusters with identical CDR3s
tmp0=split(dd0, dd0[,2])
sapply(tmp0, function(x)length(unique(x[,1])))->LL
vvL=which(LL>1)
tmp0=tmp0[vvL]

## consider a total of 5 positions, from 5 to 9, skip the last ones if CDR3 is too short, leaving at least 4 letters at the C end
pos.trimers=c()  ## All 40,000 positional trimers
AAstring='ACDEFGHIKLMNPQRSTVWY'
AAstringList=unlist(strsplit(AAstring,''))
for(pos in 5:10){
  for(A1 in AAstringList){
    for(A2 in AAstringList){
      for(A3 in AAstringList){
        pos.trimers=c(pos.trimers, paste(pos,'_',A1, A2, A3, sep=''))
      }
    }
  }
}


getTrimerPos <- function(cdr3, st=5, ed=4, MAX=6){
  nL=nchar(cdr3)
  posL=c(st:min(nL-ed, st+MAX-1))
  ss=sapply(posL, function(x)substr(cdr3, x, x+2))
  ss=unlist(ss)
  ns=length(ss)+st-1
  return(paste(st:ns, ss, sep='_'))
}

getTrimerPairs <- function(x){
  n=nrow(x)
  Tpairs=c()
  for(ii in 1:n){
    s1=getTrimerPos(x[ii,1])
    for(jj in ii:n){
      if(jj==ii)next
      s2=getTrimerPos(x[jj,1])
      vv=which(s1 != s2)
      Tpairs=rbind(Tpairs, cbind(s1[vv], s2[vv]))
    }
  }
  return(Tpairs)
}

TpairMat=c()
count=0
for(x in tmp0){
  if(count %% 10000 == 0)print(count)
  Tp=getTrimerPairs(x)
  TpairMat=rbind(TpairMat, Tp)
  count=count+1
}


vv=grep('10_',TpairMat[,1])
TT=TpairMat0[vv,]
TT=apply(TT, 1, sort)
TT=t(TT)

tmp.tt=table(TT[,1], TT[,2])
ss=unique(c(TT[,1], TT[,2]))
ns=length(ss)
tmp.mat=matrix(0, ns, ns)
colnames(tmp.mat)=rownames(tmp.mat)=ss
for(kk in 1:nrow(tmp.tt)){
  tmp.mat[rownames(tmp.tt)[kk], colnames(tmp.tt)]=tmp.tt[kk,]
}


tmp.matp=t(tmp.mat)
tmp.mat0=tmp.mat+tmp.matp
cor.mat=cor(tmp.mat0, method='s')
fit = cmdscale(sqrt(1-cor.mat), k=100)


if(1){  ## make graph plot
tmp.mat1=tmp.mat
tmp.mat1[which(tmp.mat1<=50)]=0
tmp.mat1[which(tmp.mat1>50)]=1

#vv=grep('^G|^R', rownames(tmp.mat1))
#tmp.mat1.sub=tmp.mat1[vv,vv]

gg=graph_from_adjacency_matrix(tmp.mat1, mode='undirected')

Isolated=which(degree(gg)<1)  ## 2 for cancer, 5 for all

g3=delete.vertices(gg, Isolated)
#V(g3)$name -> ss
#tmp.ss=apply(tmp.mat1,1,sum)[ss]
#vv=which(tmp.ss>=5 & tmp.ss<=5)[1:60]
#ss[-vv]=NA
#V(g3)$name <- ss
par(mar=c(0,0,0,0))
plot(g3, vertex.size=0.1, vertex.label.cex=1, vertex.label.font=2)

}

tmp.matp=t(tmp.mat)
tmp.mat0=tmp.mat+tmp.matp

## perform imputation  ## doesn't work
library(softImpute)
tmp.mat0[which(tmp.mat0==0)]=NA

fit4=softImpute(tmp.mat0,rank=50,lambda=10)
ximp=complete(xna,fit4)
impute(fit4,i=c(1,3,7),j=c(2,5,10))
impute(fit4,i=c(1,3,7),j=c(2,5,10),unscale=FALSE)#ignore scaling and centering

## perform MDS
cor.mat=cor(tmp.mat0, method='s')
fit = cmdscale(cor.mat, k=100)


## calculate distance
getDist <- function(x1, x2){
  delta = x1 - x2
  return(sum(delta^2))
}

getTrimerDist <- function(s1, s2){
  x1=fit[s1,]
  x2=fit[s2,]
  return(getDist(x1, x2))
}

EncodeRepertoire <- function(ff){
  dd1=read.table(ff, header=F,sep='\t',stringsAsFactors = F)
  vv=grep('^C',dd1[,1])
  dd1=dd1[vv,]
  ssL=sapply(dd1[,1], getTrimerPos, st=5,ed=4)
  #sapply(sapply(dd1[,1], getTrimerPos, st=5,ed=4), function(x)tt0[intersect(x, names(tt0))]) -> tmp
  #sapply(tmp, function(x)x^0/sum(x^0)) ->tmp1
  # for(x in tmp1){
  #   if(length(which(is.na(x)))>0){
  #     print(x)
  #     x.mat=cbind(x.mat, NA)
  #     next
  #   }
  #   ss=gsub('.+_','',names(x))
  #   xs=fit[ss,]
  #   x %*% xs -> xx
  #   x.mat=rbind(x.mat, xx)
  # }
  foo <- function(ss){
    ss=gsub('.+_','',ss)
    ss=intersect(ss, rownames(fit))
    if(length(ss)==1)return(fit[ss,])
    xs=t(fit[ss,])
    xx = rowMeans(xs)
    xx
  }
  x.mat=sapply(ssL, foo)
  #x.mat=do.call(rbind, x.mat)
  return(t(x.mat))
}

dda=read.table('TCRantigen_Processed_selectedFeb25_2021.txt',header=F,sep='\t')
tt=table(c(TpairMat0[,1], TpairMat0[,2]))
tt0=1/(1+tt)
tmp.ss=sapply(dda[,1], getTrimerPos, st=5,ed=4)
sapply(sapply(tmp.ss, function(x)gsub('.+_','',x)), function(x)tt0[intersect(x, names(tt0))]) -> tmp
sapply(tmp, function(x)x^0/sum(x^0)) ->tmp1

getCoord<-function(ss, tt=tt0, Beta=fit){
  xs=getTrimerPos(ss, ed=4)
  xs=intersect(xs, names(tt))
  xx=tt0[xs]
  xx=xx/sum(xx)
  xs0=gsub('.+_','',xs)
  return( xx %*% Beta[xs0,])
}

x.mat=c()
for(x in tmp1){
  if(length(which(is.na(x)))>0){
    print(x)
    x.mat=cbind(x.mat, NA)
    next
  }
  ss=gsub('.+_','',names(x))
  xs=fit[ss,]
  x %*% xs -> xx
  x.mat=rbind(x.mat, xx)
}

vv=order(dda[,5], decreasing=F)
dd1=dda[vv,]
x.mat1=x.mat[vv,]

Nmat=nrow(x.mat1)
tmpDist=c()
for(ii in 1:Nmat){
  if(ii %% 100 ==0)print(ii)
  x1=x.mat1[ii,]
  n1=nchar(dd1[ii,1])
  for(jj in ii:Nmat){
    if(jj==ii)next
    x2=x.mat1[jj,]
    n2=nchar(dd1[jj,1])
    if(abs(n1-n2)>=3)next
    tmp.d=getDist(x1, x2)
    tmpDist=rbind(tmpDist, c(tmp.d, dd1[ii,5], dd1[jj, 5]))
  }
}

vv=which(tmpDist1000[,2]==tmpDist1000[,3])
xx1=as.numeric(tmpDist1000[vv,1])
xx2=as.numeric(tmpDist1000[-vv,1])
boxplot(xx1,xx2)

tmp.x=c(xx1, xx2)
tmp.y=c(rep(1, length(xx1)), rep(0,length(xx2)))

plot(roc(tmp.y ~ tmp.x))

TS=Rtsne(x.mat, check_duplicates=F)


## Visualize differential densities of two repertoires
ff1='/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/Shipp/Samples/for_GIANA/TestReal-Donor_ZC.tsv'
ff2='../../TCR_peptide/data/Adaptive/Shipp/Lymphoma/PBMC/for_GIANA/TestReal-ID12_P121_PreCycle2_PBMC.tsv'
#ff1='/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/COVID19/Samples/for_GIANA/TestReal-KHBR20-00205_TCRB.tsv'
#ff2='/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/Control2020/Samples/for_GIANA/TestReal-Subject_41.tsv'
#ff2='../../TCR_peptide/data/Adaptive/Shipp/Lymphoma/PBMC/for_GIANA/TestReal-ID26_P161_PreCycle2_PBMC.tsv'
#ff1='../../TCR_peptide/data/Adaptive/EarlyBRCA/Sample/for_GIANA/TestReal-BR01T.tsv'
#ff2='../../TCR_peptide/data/Adaptive/EarlyBRCA/Sample/for_GIANA/TestReal-BR01N.tsv'
d1=EncodeRepertoire(ff1)
d2=EncodeRepertoire(ff2)
TMP0=Rtsne(rbind(d1, d2), check_duplicates=F)
n1=nrow(d1)
n2=nrow(d2)
kde2d(TMP0$Y[1:n1,1], TMP0$Y[1:n1,2],n=400)->dens1
kde2d(TMP0$Y[(n1+1):(n1+n2),1], TMP0$Y[(n1+1):(n1+n2),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens2
image(dens1$x, dens1$y, dens1$z - dens2$z, col=redblue,zlim=c(-0.0004, 0.0004))
range(dens1$z - dens2$z)

  ## Figure 1 contour plots
kde2d(TMP0$Y[1:n1,1], TMP0$Y[1:n1,2],n=200)->dens1
image(dens1, col=redblue, axes=F)
contour(dens1, lwd=3, labcex=1, col=1, add=T)
kde2d(TMP0$Y[(n1+1):(n1+n2),1], TMP0$Y[(n1+1):(n1+n2),2],n=200, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens2
image(dens2, col=redblue, axes=F)
contour(dens2, lwd=3, labcex=1, col=1, add=T)
image(dens1$x, dens1$y, dens1$z - dens2$z, col=redblue,zlim=c(-0.0002, 0.0002))
contour(dens1$x, dens1$y, dens1$z - dens2$z, col=1, lwd=3, labcex=1, add=T)
vv=which(TMP0$Y[,1]>=-5 & TMP0$Y[,1]<=0 & TMP0$Y[,2]>=36 & TMP0$Y[,2]<=40)
rownames(d1)[vv]

vv=which(TMP0$Y[,1]>=0 & TMP0$Y[,1]<=8 & TMP0$Y[,2]>=-17 & TMP0$Y[,2]<= -12)
rownames(rbind(d1,d2))[vv]


ffs=dir('../../TCR_peptide/data/Adaptive/OV2021/for_GIANA/',full.names=T)
ff.control=ffs[grep('Carter',ffs)][1:10]
dd.control=c()

for(ff in ff.control){
  print(ff)
  d0=EncodeRepertoire(ff)
  dd.control=rbind(dd.control, d0)
}

ff.ov=ffs[grep('GO',ffs)][52:55]

dd.ov=c()

for(ff in ff.ov){
  print(ff)
  d0=EncodeRepertoire(ff)
  dd.ov=c(dd.ov, list(d0))
}

dd.c1=EncodeRepertoire(ffs[grep('OV1072',ffs)])
dd.c2=EncodeRepertoire(ffs[grep('OV1043',ffs)])

TMP0=Rtsne(rbind(dd.control, dd.c1, dd.c2, dd.ov[[1]],dd.ov[[2]], dd.ov[[3]], dd.ov[[4]]), check_duplicates=F, num_threads=0)
#TMP0=umap(rbind(dd.control, dd.c1, dd.c2, dd.ov[[1]],dd.ov[[2]], dd.ov[[3]], dd.ov[[4]]))

n1=nrow(dd.control)
n2=nrow(dd.c1)
n3=nrow(dd.c2)
n4=nrow(dd.ov[[1]])
n5=nrow(dd.ov[[2]])
n6=nrow(dd.ov[[3]])
n7=nrow(dd.ov[[4]])


kde2d(TMP0$Y[1:n1,1], TMP0$Y[1:n1,2],n=400)->dens1
kde2d(TMP0$Y[(n1+1):(n1+n2),1], TMP0$Y[(n1+1):(n1+n2),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens2
kde2d(TMP0$Y[(n1+n2+1):(n1+n2+n3),1], TMP0$Y[(n1+n2+1):(n1+n2+n3),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens3
kde2d(TMP0$Y[(n1+n2+n3+1):(n1+n2+n3+n4),1], TMP0$Y[(n1+n2+n3+1):(n1+n2+n3+n4),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens4
kde2d(TMP0$Y[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),1], TMP0$Y[(n1+n2+n3+n4+1):(n1+n2+n3+n4+n5),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens5
kde2d(TMP0$Y[(n1+n2+n3+n4+n5+1):(n1+n2+n3+n4+n5+n6),1], TMP0$Y[(n1+n2+n3+n4+n5+1):(n1+n2+n3+n4+n5+n6),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens6
kde2d(TMP0$Y[(n1+n2+n3+n4+n5+n6+1):(n1+n2+n3+n4+n5+n6+n7),1], TMP0$Y[(n1+n2+n3+n4+n5+n6+1):(n1+n2+n3+n4+n5+n6+n7),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens7


par(mar=c(0,0,0,0),mfrow=c(2,3))
image(dens2$x, dens2$y, dens2$z - dens1$z, col=redblue, zlim=c(-0.0004,0.0004))
image(dens3$x, dens3$y, dens3$z - dens1$z, col=redblue, zlim=c(-0.0004,0.0004))
image(dens4$x, dens4$y, dens4$z - dens1$z, col=redblue, zlim=c(-0.0004,0.0004))
image(dens5$x, dens5$y, dens5$z - dens1$z, col=redblue, zlim=c(-0.0004,0.0004))
image(dens6$x, dens6$y, dens6$z - dens1$z, col=redblue, zlim=c(-0.0004,0.0004))
image(dens7$x, dens7$y, dens7$z - dens1$z, col=redblue, zlim=c(-0.0004,0.0004))


range(dens1$z - dens2$z)

par(mar=c(0,0,0,0),mfrow=c(2,2))
plot(TMP0$Y[1:n1,], col=rgb(50,50,50,alpha=20,max=255), cex=0.7, axes=F)


  ## Using Shipp data to perform a systematic differential repertoire analysis
ffs.HC=dir('../../TCR_peptide/data/Adaptive/Shipp/Samples/for_GIANA/', full.names=T)
ffs.PBMC=dir('../../TCR_peptide/data/Adaptive/Shipp/Lymphoma/PBMC/for_GIANA/', full.names=T)
ffs.PBMC=ffs.PBMC[grep('.tsv',ffs.PBMC)]

par(mar=c(0,0,2,0), mfrow=c(4,4))
for(f1 in ffs.HC[1:4]){
  f1.root=unlist(strsplit(f1, '/'))
  f1.root=f1.root[length(f1.root)]
  print(f1.root)
  for(f2 in ffs.PBMC[1:4]){
    f2.root=unlist(strsplit(f2, '/'))
    f2.root=f2.root[length(f2.root)]   
    TMP0=DiffDensityPlot(f1, f2)
    title(paste(f1.root, f2.root, sep='_'))
  }
}

  ## pre-cycle comparison
ffs.PBMC=dir('../../TCR_peptide/data/Adaptive/Shipp/Lymphoma/PBMC/for_GIANA/', full.names=T)
ffs.PBMC=ffs.PBMC[grep('.tsv',ffs.PBMC)]
ID.tt=table(str_extract(ffs.PBMC, 'ID[0-9]{2}'))
IDs=names(which(ID.tt==3))
#par(mar=c(0,0,2,0), mfrow=c(3,5))
dDO.list=c()
for(ii in IDs){
  print(ii)
  vv=grep(ii, ffs.PBMC)
  f1=ffs.PBMC[vv[1]]
  f2=ffs.PBMC[vv[2]]
  f3=ffs.PBMC[vv[3]]
#  png(file=paste('DiffPlot/',ii,'.png'), width=1000, height=1000, res=300)
#  par(mar=c(2,2,2,0))
#  dDO=DiffDensityPlot(f1, f2, ZLIM=0.0001)
#  title(ii)
#  dev.off()
#  dDO.list=c(dDO.list, list(dDO))
  d1=EncodeRepertoire(f1)
  d2=EncodeRepertoire(f2)
  d3=EncodeRepertoire(f3)
  d1=na.omit(d1)
  d2=na.omit(d2)
  d3=na.omit(d3)
  dd=rbind(d1, d2, d3)
  TMP0=Rtsne(dd, check_duplicate=F)
  n1=nrow(d1)
  n2=nrow(d2)
  n3=nrow(d3)
  kde2d(TMP0$Y[1:n1,1], TMP0$Y[1:n1,2],n=400)->dens1
  kde2d(TMP0$Y[(n1+1):(n1+n2),1], TMP0$Y[(n1+1):(n1+n2),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens2
  kde2d(TMP0$Y[(n1+n2+1):(n1+n2+n3),1], TMP0$Y[(n1+n2+1):(n1+n2+n3),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens3
  dDO=c(tsne=list(TMP0), dens1=list(dens1), dense2=list(dens2), dens3=list(dens3),n1=n1, n2=n2)
  dDO.list=c(dDO.list, list(dDO))
}

par(mfrow=c(5,8), mar=c(0,0,2,0))
names(dDO.list)=IDs
for(ii in 1:40){
  dDO=dDO.list[[ii]]
  if(dDO$n1<=5000)next
  image(dDO$dens1, col=redblue, zlim=c(0,0.0005), main=IDs[ii], cex.main=1)
  #image(dDO$dens1$x, dDO$dens1$y, dDO$dense2$z - dDO$dens1$z, col=redblue, main=IDs[ii], zlim=c(-0.0001,0.0001))
}

par(mfrow=c(3,9),mar=c(0,0,2,0))
for(id in IDs[which(shipp.info[IDs,4]==2)]){
  ii=dDO.list[[id]]
  image(ii$dense2$x, ii$dense2$y, ii$dens3$z-ii$dense2$z, col=redblue, main=id, zlim=c(-0.0001,0.0001))
}

boxplot(sapply(dDO.list, function(x)max(x$dense2-x$dens3)) ~ shipp.info[IDs,3])
  ## Put 4 CR, 4 PD and 4 HC together to visualize on tsne
tmp.CRs=c('ID12','ID13','ID14','ID15')
tmp.PDs=c('ID28','ID29','ID30','ID31')
tmp.HCs=c('BC','KW','SB','RK')

dd.CRPD=c()
NNs=c()
for(id in c(tmp.CRs, tmp.PDs)){
  ff=ffs.PBMC[grep(id,ffs.PBMC)]
  ff=ff[grep('baseline',ff)]
  dd=EncodeRepertoire(ff)
  NNs=c(NNs, nrow(dd))
  dd.CRPD=rbind(dd.CRPD, dd)
}

for(id in tmp.HCs){
  ff=ffs.HC[grep(id, ffs.HC)]
  dd=EncodeRepertoire(ff)
  NNs=c(NNs, nrow(dd))
  dd.CRPD=rbind(dd.CRPD, dd)
}

tsne.CRPD=Rtsne(dd, check_duplicates=F)
cNNs=cumsum(NNs)
kde2d(tmp$Y[72111:112086,1], tmp$Y[72111:112086,2], n=400, lims =c(range(tmp$Y[,1]), range(tmp$Y[,2]))) -> dens0
par(mfrow=c(2,4), mar=c(0,0,2,0))
for(ii in 1:8){
  if(ii==1)vv=1:cNNs[1] else vv=cNNs[ii-1]:cNNs[ii]
  tmp.dens=kde2d(tmp$Y[vv,1], tmp$Y[vv,2], n=400, lims =c(range(tmp$Y[,1]), range(tmp$Y[,2])))
  image(dens0$x, dens0$y, tmp.dens$z-dens0$z, col=redblue, zlim=c(-0.0004,0.0004), main=names(NNs)[ii])
}


tmp.CRs=c('ID32','ID33','ID60','ID62')
tmp.PDs=c('ID65','ID66','ID36','ID28')
dd.CRPD=c()
NNs=c()
ffs.PBMC.CD4=dir('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/Shipp/Lymphoma/for_GIANA/', full.names=T)
ffs.PBMC.CD4=ffs.PBMC.CD4[grep('baseline_CD4', ffs.PBMC.CD4)]
for(id in c(tmp.CRs, tmp.PDs)){
  ff=ffs.PBMC.CD4[grep(id,ffs.PBMC.CD4)]
  dd=EncodeRepertoire(ff)
  NNs=c(NNs, nrow(dd))
  dd.CRPD=rbind(dd.CRPD, dd)
}

tsne.CRPD.CD4=Rtsne(dd.CRPD, check_duplicates=F)

names(NNs)=c(tmp.CRs, tmp.PDs)
cNNs=cumsum(NNs)
par(mar=c(0,0,2,0), mfrow=c(2,4))
for(ii in 1:8){
  if(ii==1)vv=1:cNNs[1] else vv=cNNs[ii-1]:cNNs[ii]
  plot(tsne.CRPD.CD4$Y[vv,], pch=19, cex=0.5, col=rgb(50,50,50,alpha=50,max=255), main=names(NNs)[ii])
}

## Calculate KL divergence between patients and HC
KL.mat=matrix(NA, ncol=length(ffs.HC), nrow=length(ffs.PBMC))
rownames(KL.mat)=ffs.PBMC
colnames(KL.mat)=ffs.HC
for(ff in ffs.HC[2:10]){
  print(ff)
  ddc=EncodeRepertoire(ff)
  for(id in IDs){
    print(id)
    vv=grep(id, ffs.PBMC)
    flag=0
    for(ii in vv){
      if(!is.na(KL.mat[ffs.PBMC[ii], ff]))next
      ddp=EncodeRepertoire(ffs.PBMC[ii])
      if(nrow(ddp)<=5000){
        flag=1
        break
      }
      dd=rbind(ddc, ddp)
      tmp=Rtsne(dd, check_duplicates=FALSE)
      n1=nrow(ddc)
      n2=nrow(ddp)
      xy=tmp$Y
      dens1=kde2d(xy[1:n1,1], xy[1:n1,2], n=400, lims=c(range(xy[,1]), range(xy[,2])))
      dens2=kde2d(xy[(n1+1):(n1+n2),1], xy[(n1+1):(n1+n2),2], n=400, lims=c(range(xy[,1]), range(xy[,2])))
      tmp.KL=KLD(dens1$z, dens2$z)
      KL.mat[ffs.PBMC[ii], ff]=tmp.KL$mean.sum.KLD
    }
    if(flag==1)next
  }
}

plot(0,0,xlim=c(1,3),ylim=c(0.01,0.08), cex=0)
for(id in IDs){
  vv.id=grep(id, rownames(KL.mat))
  tmp.col=shipp.info[id,4]
  if(shipp.info[id,3]=='PD')tmp.col=3
  if(tmp.col==2)next
  points(KL.mat[vv.id, 2], pch=19, type='b', col=tmp.col)
}

plot(0,0,xlim=c(1,3),ylim=c(0.01,0.08), cex=0)
tmp.cycleDiff=c()
for(id in IDs){
  vv.id=grep(id, rownames(KL.mat))
  tmp.col=shipp.info[id,4]
  if(shipp.info[id,3]=='PD')tmp.col=3
  tmp.cycleDiff=c(tmp.cycleDiff, KL.mat[vv.id[3],1] - KL.mat[vv.id[2],1])
}
names(tmp.cycleDiff)=IDs
boxplot(tmp.cycleDiff ~ shipp.info[IDs, 3])

## Position new points onto existing tSNE: doesn't work well
PositionNew <- function(dn, dr, ref.coord, k=10){
  ## dn: new data
  ## dr: ref data
  tmp.cor=cor(t(dn), t(dr))
  k.nn <- function(x){
    xs=sort(x, decreasing=T)
    v=which(x >= xs[k])
    new.d=ref.coord[v,]
    return(apply(new.d, 2, median))
  }
  new.coord= apply(tmp.cor, 1, k.nn)
  return(t(new.coord))
}


## Build mega control to define "genes of TCR repertoire"
ffs=dir('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/sampleExport/batch1/for_GIANA/',full.names=T)
#ff.control=ffs[grep('Carter',ffs)][1:10]
ff.control0=sample(ffs, 50)
dd.control0=c()

for(ff in ff.control0){
  print(ff)
  d0=EncodeRepertoire(ff)
  if(nrow(d0)>10000)break
  dd.control0=rbind(dd.control0, d0)
}

dd.control0=na.omit(dd.control0)

TMP0c=Rtsne(dd.control0, check_duplicates=F)

AssignRFUs <- function(ff, CL=km5000, THR=0.6){
  ## Assign repertoire functional units defined by k-means cluster centroids from pooled control samples
  dd=EncodeRepertoire(ff)
  ddc=CL$centers
  dd.cor = cor(t(dd), t(ddc), method='s')
  tmp=apply(dd.cor, 1, max)
  vv=which(tmp>=THR)
  tmpv=apply(dd.cor, 1, function(x)which.max(x)[1])
  tt=table(tmpv)
  tt0=table(CL$cl)
  RFU = rep(0, length(tt0))
  names(RFU)=1:length(RFU)
  RFU[names(tt)] = tt
  #RFU = RFU / as.numeric(tt0)  ## normalize by original cluster sizes
  RFU = RFU / nrow(dd) * 10000 ## normalize by total TCR in the data
  Nmiss= nrow(dd) - length(vv)
  list(RFU= RFU, N=Nmiss, COR=tmp)
}

## Calling RFUs from different datasets
dd.RFU =c()
N.miss=c()
ffs=dir('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/',full.names=T)
for(ff in ffs[1:length(ffs)]){
  print(ff)
  RFU = AssignRFUs(ff)
  dd.RFU=cbind(dd.RFU, RFU$RFU)
  N.miss=c(N.miss, RFU$N)
}

colnames(dd.RFU)=gsub('.+TestReal-','',ffs)
colnames(dd.RFU)=gsub('.tsv','',colnames(dd.RFU))
vvc=which(OVlabels[colnames(dd.RFU)]!='highGrade')
vvo=which(OVlabels[colnames(dd.RFU)]=='highGrade')
pList=c()
FClist=c()
for(ii in 1:nrow(dd.RFU)){
  x1=dd.RFU[ii,vvc]
  x2=dd.RFU[ii,vvo]
  tmp.test=t.test(x1,x2 )
  pList=c(pList, tmp.test$p.value)
  FClist=c(FClist, mean(x1)/mean(x2))
}

ovRFU = dd.RFU

dd.RFU =c()
tmp.N=c()
ffs=dir('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/Melanoma_Valpione_2020/Samples/PBMC/for_GIANA/', full.names=T)
for(ff in ffs[1:length(ffs)]){
  print(ff)
  if(length(grep('[tc]sv|txt',ff))==0)next
  RFU = AssignRFUs(ff)
  dd.RFU=cbind(dd.RFU, RFU$RFU)
  tmp.N=c(tmp.N, RFU$N)
}

colnames(dd.RFU)=gsub('.+TestReal-','',ffs)
colnames(dd.RFU)=gsub('.[tc]sv|txt','',colnames(dd.RFU))

## Joining gene usage
J.RFU =c()
ffs=dir('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/sampleExport/', full.names=T)
ffs=ffs[grep('.tsv',ffs)]
for(ff in ffs[1:length(ffs)]){
  print(ff)
  if(length(grep('.tsv',ff))==0)next
  tmp.dd=read.table(ff, header=T, sep='\t')
  tmp.xx=table(tmp.dd[,'jMaxResolved'])
  tmp.xx=as.numeric(tmp.xx)/sum(tmp.xx)
  J.RFU=cbind(J.RFU, tmp.xx)
}

colnames(J.RFU)=gsub('.+TestReal-','',ffs)
colnames(J.RFU)=gsub('.[tc]sv|txt','',colnames(dd.RFU))


pList=c()
FClist=c()
vvc=which(OVlabels[colnames(ovRFU)]!='highGrade')
vvo=which(OVlabels[colnames(ovRFU)]=='highGrade')
for(ii in 1:nrow(ovRFU)){
  x1=ovRFU[ii,vvc]
  x2=ovRFU[ii,vvo]
  tmp.test=t.test(x1,x2 )
  pList=c(pList, tmp.test$p.value)
  FClist=c(FClist, mean(x1)/mean(x2))
}


pList=c()
FClist=c()
vvc=which(unlist(hcmvData[colnames(dd.RFU),5])==0)
vvo=which(unlist(hcmvData[colnames(dd.RFU),5])==1)
for(ii in 1:nrow(dd.RFU)){
  x1=dd.RFU[ii,vvc]
  x2=dd.RFU[ii,vvo]
  tmp.test=wilcox.test(x1,x2 )
  pList=c(pList, tmp.test$p.value)
  FClist=c(FClist, mean(x1)/mean(x2))
}

## Unique IDs for COVID data
RFU.covid=cbind(RFU.covid.Adaptive, RFU.covid.KHBR, RFU.covid.INCOV)
vvn=which(table(km5000$cl)>=20)
tt.ID=table(covid.info[colnames(RFU.covid),'subject_id'])
sel.IDs=c()
for(ii in names(tt.ID)){
  vv.id=which(covid.info[,'subject_id']==ii)
  tmp.id=covid.info[vv.id,'sample_name']
  tmp.id=intersect(tmp.id, colnames(RFU.covid))
  sel.IDs=c(sel.IDs,tmp.id[1])
}

RFU.covid.u=RFU.covid[vvn,sel.IDs]


pList=c()
FClist=c()
vvc=grep('Caucasian',covid.info[rownames(pca),8])
vvo=grep('Black',covid.info[rownames(pca),8])
for(ii in 1:nrow(RFU.covid.u)){
  x1=RFU.covid.u[ii,vvc]
  x2=RFU.covid.u[ii,vvo]
  tmp.test=wilcox.test(x1,x2 )
  pList=c(pList, tmp.test$p.value)
  FClist=c(FClist, mean(x1)/mean(x2))
}


## Annotate RFUs with TCRs of known specificity
dda=read.table('TCRantigenData_unique.txt', header=F,sep='\t',stringsAsFactors = F)
epitopes=gsub('.+_','',dda[,5])
tta=table(epitopes)
vva=names(which(tta>=100))

for(aa in vva){
  vv=grep(aa, dda[,5])
  fname=paste('KnownAntigens/',aa,'.txt',sep='')
  write.table(dda[vv,],file=fname, quote=F, sep='\t', row.names=F, col.names=F)
}

ffs=dir('KnownAntigens/',full.names=T)
aa.RFU=c()
for(ff in ffs){
  print(ff)
  TMP=AssignRFUs(ff)
  aa.RFU=cbind(aa.RFU, TMP$RFU)
}

colnames(aa.RFU)=gsub('.+//','',ffs)
colnames(aa.RFU)=gsub('.txt','',colnames(aa.RFU))

RFUannotationTable <- matrix(NA, ncol=3, nrow=5000)
colnames(RFUannotationTable)=c('ID','epitope','species')
RFUannotationTable[,1]=1:5000
for(ii in 1:ncol(aa.RFU)){
  tmp.xx=aa.RFU[,ii]
  tmp.aa=colnames(aa.RFU)[ii]
  tmp.vv=grep(tmp.aa,dda[,5])[1]
  tmp.species=gsub('Antigen:','',dda[tmp.vv,5])
  tmp.species=gsub('_.+','',tmp.species)
  vvx=which(tmp.xx>=1)
  for(vv in vvx){
    if(is.na(RFUannotationTable[vv,2])){
      RFUannotationTable[vv,2]=tmp.aa
      RFUannotationTable[vv,3]=tmp.species
    }else{
      RFUannotationTable[vv,2]=paste(RFUannotationTable[vv,2], tmp.aa, sep=';')
      RFUannotationTable[vv,3]=paste(RFUannotationTable[vv,3], tmp.species, sep=';')
    }
  }
}


## Shipp data
pList.shipp=c()
FC.shipp=c()

#vv.res=rownames(Shipp.info[which(Shipp.info[,3] == 'CR'),])
#vv.non=rownames(Shipp.info[which(Shipp.info[,3] == 'PR'),])
vv1=names(which(Shipp.labels=='HC'))
vv2=names(which(Shipp.labels=='Cancer-baseline'))

for(ii in 1:nrow(RFU.shipp)){
  x1=RFU.shipp[ii,vv1]
  x2=RFU.shipp[ii,vv2]
  tmp=t.test(x1, x2)
  pList.shipp=c(pList.shipp, tmp$p.value)
  FC.shipp=c(FC.shipp, mean(x1)/mean(x2))
}


vv1=grep('CD4',colnames(RFU.shipp.CD4.CD8))
vv2=grep('CD8', colnames(RFU.shipp.CD4.CD8))

pList.shipp=c()
FC.shipp=c()
for(ii in 1:nrow(dd.RFU)){
  x1=RFU.shipp.CD4.CD8[ii,vv1]
  x2=RFU.shipp.CD4.CD8[ii,vv2]
  tmp=wilcox.test(x1, x2, paired=T)
  pList.shipp=c(pList.shipp, tmp$p.value)
  FC.shipp=c(FC.shipp, median(x1)/median(x2))
}
plot(log(FC.shipp), -log(pList.shipp))


vv1=grep('CD4',colnames(dd.RFU))
vv2=grep('CD8', colnames(dd.RFU))

pList=c()
FClist=c()
for(ii in 1:nrow(RFU.covid)){
  x1=RFU.covid[ii,]
  x2=RFU.lung.MDA[ii,]
  tmp=t.test(x1, x2)
  pList=c(pList, tmp$p.value)
  FClist=c(FClist, mean(x1)/mean(x2))
}
plot(log(FClist), -log(pList))

pList.irep=c()
FClist.irep=c()
for(ii in 1:nrow(RFU.iRep.rcc)){
  x1=RFU.iRep.rcc[ii,]
  x2=RFU.iRep.Normal[ii,]
  tmp=t.test(x1, x2)
  pList.irep=c(pList.irep, tmp$p.value)
  FClist.irep=c(FClist.irep, mean(x1)/mean(x2))
}
plot(log(FClist.irep), -log(pList.irep))


pList.ov=c()
FClist.ov=c()
vv1=which(OVlabels[colnames(RFU.ov)]=='Normal')
vv2=which(OVlabels[colnames(RFU.ov)]=='highGrade')

for(ii in 1:nrow(RFU.ov)){
  x1=RFU.ov[ii,vv1]
  x2=RFU.ov[ii,vv2]
  tmp=t.test(x1, x2)
  pList.ov=c(pList.ov, tmp$p.value)
  FClist.ov=c(FClist.ov, mean(x1)/mean(x2))
}
plot(log(FClist.ov), -log(pList.ov))


pList.ov=c()
FClist.ov=c()
vv1=which(OVlabels[colnames(RFU.ov)]=='Normal')
vv2=which(OVlabels[colnames(RFU.ov)]=='highGrade')

for(ii in 1:nrow(RFU.ov)){
  x1=RFU.ov[ii,vv1]
  x2=RFU.ov[ii,vv2]
  tmp=t.test(x1, x2)
  pList.ov=c(pList.ov, tmp$p.value)
  FClist.ov=c(FClist.ov, mean(x1)/mean(x2))
}
plot(log(FClist.ov), -log(pList.ov))

pList.cerv=c()
FClist.cerv=c()
vv1=which(cerv.lab[colnames(RFU.cerv)]==1)
vv2=which(cerv.lab[colnames(RFU.cerv)]==2)

for(ii in 1:nrow(RFU.cerv)){
  x1=RFU.cerv[ii,vv1]
  x2=RFU.cerv[ii,vv2]
  tmp=t.test(x1, x2)
  pList.cerv=c(pList.cerv, tmp$p.value)
  FClist.cerv=c(FClist.cerv, mean(x1)/mean(x2))
}
plot(log(FClist.cerv), -log(pList.cerv))


pList.lung=c()
FClist.lung=c()

for(ii in 1:nrow(RFU.lung.MDA)){
  x1=RFU.covid[ii,]
  x2=RFU.lung.MDA[ii,]
  tmp=t.test(x1, x2)
  pList.lung=c(pList.lung, tmp$p.value)
  FClist.lung=c(FClist.lung, mean(x1)/mean(x2))
}
names(pList.lung)=1:5000
plot(log(FClist.lung), -log(pList.lung))

#tmp.dd=cbind(RFU.covid[,601:ncol(RFU.covid)], RFU.lung.MDA)

## predict ovarian cancer
tmp.mm.cerv=names(which(p.adjust(pList.cerv,'BH')<=0.05 & FClist.cerv<1))
tmp.mm.shipp=names(which(p.adjust(pList.shipp,'BH')<=0.1 & FClist.shipp<1))
tmp.mm.ov=names(which(p.adjust(pList.ov,'BH')<=0.1 & FClist.ov>1))
tmp.mm.lung=names(which(p.adjust(pList.lung,'BH')<=0.1 & FClist.lung<1))


tmp.mm.shipp.u=names(which(p.adjust(pList.shipp,'BH')<=0.1 & FClist.shipp>1))
tmp.mm.ov.u=names(which(p.adjust(pList.ov,'BH')<=0.1 & FClist.ov<1))
tmp.mm.lung.u=names(which(p.adjust(pList.lung,'BH')<=0.1 & FClist.lung>1))


tmp.mm.ov=names(sort(pList.ov)[1:30])
tmp.mm.shipp=names(sort(pList.shipp)[1:30])
tmp.mm.lung=names(sort(pList.lung)[1:200])
tmp.mm.cerv=names(sort(pList.cerv)[1:30])
tmp.mm.irep=names(sort(pList.irep))[1:30]

## train and predict lung cancer
vv.train.n=sample(colnames(RFU.covid), 700)
vv.train.t=colnames(RFU.lung.MDA)

vv.test.n=colnames(RFU.covid)[which(! colnames(RFU.covid) %in% vv.train.n)]
vv.test.t=colnames(RFU.lung.cascone)[grep('blood', colnames(RFU.lung.cascone))]

pList.train=c()
for(ii in 1:nrow(RFU.covid)){
  tmp=t.test(RFU.covid[ii, vv.train.n], RFU.lung.MDA[ii,])
  pList.train=c(pList.train, tmp$p.value)
}
names(pList.train)=1:5000

tmp.mm = names(sort(pList.train))[1:50]

#tmp.dd.test=cbind(RFU.covid[,vv.test.n], RFU.lung.cascone[,vv.test.t], RFU.control[,vv.lung])

#tmp.lab=c(rep(1, length(vv.test.n)), rep(2, length(vv.test.t)+22))
tmp.dd.test=cbind(RFU.covid[,vv.test.n], RFU.lung.cascone[,vv.test.t])

tmp.lab=c(rep(1, length(vv.test.n)), rep(2, length(vv.test.t)))
tmp.pca=prcomp(t(tmp.dd.test[tmp.mm,]))$x

boxplot(tmp.pca[,1] ~ tmp.lab)
roc(tmp.lab ~ tmp.pca[,1])

## Cross validation based prediction
  ## training HC and lung cancer samples
vv.train.n=sample(colnames(RFU.covid), 700)
vv.train.t=sample(colnames(RFU.lung.MDA),60)
  ## validation HC and lung cancer samples
vv.test.n0=colnames(RFU.covid)[which(! colnames(RFU.covid) %in% vv.train.n)]
vv.test.n=sample(vv.test.n0,500)
vv.test.t=colnames(RFU.lung.MDA)[which(! colnames(RFU.lung.MDA) %in% vv.train.t)]
  ## independent test HC and lung cancer samples
vv.test.n.ind = vv.test.n0[which( ! vv.test.n0 %in% vv.test.n)]
vv.test.t.ind=colnames(RFU.lung.cascone)[grep('blood', colnames(br.RFU.lung.cascone))]

  ## obtain predictive markers
pList.train=c()
for(ii in 1:nrow(RFU.covid)){
  tmp=t.test(RFU.covid[ii, vv.train.n], RFU.lung.MDA[ii,vv.train.t])
  pList.train=c(pList.train, tmp$p.value)
}
names(pList.train)=1:5000

tmp.mm = names(sort(pList.train))[1:50]

tmp.dd.train=cbind(RFU.covid[,vv.train.n], RFU.lung.MDA[,vv.train.t])
tmp.lab.train=c(rep(0, length(vv.train.n)), rep(1, length(vv.train.t)))

svm.fit=svm(y=tmp.lab.train, x= t(tmp.dd.train[tmp.mm,]))
predict(svm.fit, t(tmp.dd.train[tmp.mm,])) -> val.train  ## training accuracy
roc(tmp.lab.train ~ val.train)

predict(svm.fit, t(tmp.dd.val[tmp.mm,])) -> val.predict  ## validation accuracy
roc(tmp.lab.val ~ val.predict)


tmp.dd.test=cbind(RFU.covid[,vv.test.n.ind], br.RFU.lung.cascone[,vv.test.t.ind])

tmp.lab.test=c(rep(1, length(vv.test.n.ind)), rep(2, length(vv.test.t.ind)))
predict(svm.fit, t(tmp.dd.test[tmp.mm,])) -> val.test
roc(tmp.lab.test ~ val.test)  ## test accuracy

  ## validation performance
tmp.dd.val = cbind(RFU.covid[,vv.test.n], RFU.lung.MDA[,vv.test.t])
tmp.lab.val=c(rep(1, length(vv.test.n)), rep(2, length(vv.test.t)))
tmp.pca=prcomp(t(tmp.dd.val[tmp.mm,]))$x
boxplot(tmp.pca[,1] ~ tmp.lab.val)
roc(tmp.lab.val ~ tmp.pca[,1])


  ## independent test
tmp.dd.test=cbind(RFU.covid[,vv.test.n.ind], br.RFU.lung.cascone[,vv.test.t.ind])

tmp.lab.test=c(rep(1, length(vv.test.n.ind)), rep(2, length(vv.test.t.ind)))
tmp.pca=prcomp(t(tmp.dd.test[tmp.mm,]))$x

boxplot(tmp.pca[,1] ~ tmp.lab)
roc(tmp.lab ~ tmp.pca[,1])

tmp.dd.test=cbind(RFU.covid[,vv.test.n.ind], br.RFU.lung.jhu)

tmp.lab=c(rep(1, length(vv.test.n.ind)), rep(2, ncol(br.RFU.lung.jhu)))
tmp.pca=prcomp(t(tmp.dd.test[tmp.mm,]))$x

boxplot(tmp.pca[,1] ~ tmp.lab)
roc(tmp.lab ~ tmp.pca[,1])

## using cascone as training
vv.train.n=sample(colnames(RFU.covid), 700)
vv.train.t=colnames(RFU.lung.cascone[,grep('blood',colnames(RFU.lung.cascone))])

vv.test.n=colnames(RFU.covid)[which(! colnames(RFU.covid) %in% vv.train.n)]
vv.test.t=colnames(RFU.lung.MDA)

pList.train=c()
for(ii in 1:nrow(RFU.covid)){
  tmp=t.test(RFU.covid[ii, vv.train.n], RFU.lung.cascone[ii,])
  pList.train=c(pList.train, tmp$p.value)
}
names(pList.train)=1:5000

tmp.mm = names(sort(pList.train))[1:200]

#tmp.dd.test=cbind(RFU.covid[,vv.test.n], RFU.lung.cascone[,vv.test.t], RFU.control[,vv.lung])

#tmp.lab=c(rep(1, length(vv.test.n)), rep(2, length(vv.test.t)+22))
tmp.dd.test=cbind(RFU.covid[,vv.test.n], RFU.lung.MDA)

tmp.lab=c(rep(1, length(vv.test.n)), rep(2, length(vv.test.t)))
tmp.pca=prcomp(t(tmp.dd.test[tmp.mm,]))$x

boxplot(tmp.pca[,1] ~ tmp.lab)
roc(tmp.lab ~ tmp.pca[,1])

## prepare JHU lung cancer data
ffs=dir('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/inHouse/lung/Samples/', full.names=T)
ffs=ffs[grep('.txt',ffs)]
for(ff in ffs){
  print(ff)
  dd=read.table(ff, header=F,sep='\t',stringsAsFactors = F)
  dd=ProcessAdaptiveVgenes(dd)
  N0=nrow(dd)
  N1=min(10000,N0)
  dd1=cbind(dd[1:N1,c(1,2,5,5,6)])
  tmp=unlist(strsplit(ff,'/'))
  ff0=tmp[length(tmp)]
  write.table(dd1, file=paste('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/inHouse/lung/Samples/for_GIANA/',ff0, sep=''), quote=F, sep='\t', row.names=F, col.names=F)
}



## Train and test for ovarian cancer
# Use OV_TIL vs batch2 to train model
tmp.TIL=cbind(RFU.OV_TIL)

pList.TIL=c()
for(ii in 1:5000){
  x2=tmp.TIL[ii,]
  x1=RFU.batch2[ii,]
  tmp=t.test(x1, x2)
  pList.TIL=c(pList.TIL, tmp$p.value)
}
names(pList.TIL)=1:5000
tmp.mm.TIL=names(sort(pList.TIL))[1:50]

tmp.lab.TIL=c(rep(1, 120), rep(2, ncol(tmp.TIL)))
tmp.dd.TIL=cbind(RFU.batch2, tmp.TIL)
svm.TIL = svm(y=tmp.lab.TIL, t(tmp.dd.TIL[tmp.mm.TIL,]))

predict(svm.TIL, t(tmp.dd.TIL[tmp.mm.TIL,])) -> TIL.train

boxplot(TIL.train ~ tmp.lab.TIL)

predict(svm.TIL, t(RFU.ov[tmp.mm.TIL, ])) -> TIL.test
boxplot(TIL.test ~ OVlabels[colnames(RFU.ov)])


predict(svm.TIL, t(RFU.ov[tmp.mm.TIL, vv.hgsc])) -> TIL.test
roc(OVlabels[vv.hgsc] ~ TIL.test)


## Process stage-I HGSC samples
ffs=dir('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/OV2021/stage_I/',full.names=T)
ffs=ffs[grep('.csv',ffs)]
for(ff in ffs){
  dd=read.table(ff, header=T,sep=',', stringsAsFactors = F)
  vv=which(nchar(dd[,1])<=5)
  dd=dd[-vv,]
  dd[,2]=gsub('h','',dd[,2])
  ff0=unlist(strsplit(ff,'/'))
  ff0=ff0[length(ff0)]
  dd=dd[order(dd[,ncol(dd)],decreasing=T),c(1,2,24,24)]
  dd=cbind(dd, ID=ff0)
  N0=min(nrow(dd),10000)
  dd=dd[1:N0,]
  dd[,1]=paste('C',dd[,1],'F',sep='')
  write.table(dd, file=paste('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/Adaptive/OV2021/stage_I/for_GIANA/',ff0,sep=''), sep='\t', quote=F, row.names=F, col.names=F)
}


## Study if TRBV gene allele frequency is a predictor of age or race in the COVID-19 dataset

getTRBVmat <- function(DIR){
  ffs=dir(DIR, full.names=T)
  ffs=ffs[grep('tsv',ffs)]
  TRBV.list=c()
  for(ff in ffs){
    dd=read.table(ff, header=T, sep='\t',stringsAsFactors = F)
    tmp.tt=table(dd[,2])/nrow(dd)
    TRBV.list=c(TRBV.list, list(tmp.tt))
  }
  names(TRBV.list)=gsub('.+TestReal-','',ffs)
  
  unique(unlist(sapply(TRBV.list, names))) -> tmp.uu
  
  TRBV.mat=matrix(0, ncol=length(ffs), nrow=length(tmp.uu))
  colnames(TRBV.mat)=names(TRBV.list)
  rownames(TRBV.mat)=tmp.uu
  
  for(ID in names(TRBV.list)){
    Txx=unlist(TRBV.list[ID])
    names(Txx)=gsub('.+.tsv.','',names(Txx))
    TRBV.mat[names(Txx), ID]=Txx
  }
  
  colnames(TRBV.mat)=gsub('.tsv','',colnames(TRBV.mat))
  
  TRBV.mat
}


trbv.pca=prcomp(t(TRBV.mat))$x
boxplot(trbv.pca[,1] ~ covid.info[rownames(trbv.pca),8])


Prepare.iRepIGH <- function(indir, outdir){
  ffs=dir(indir, full.names=T)
  ffs=ffs[grep('.csv|.tsv',ffs)]
  for(ff in ffs){
    ff0=strsplit(ff0,'/')
    ff0=ff0[length(ff0)]
    print(ff0)
    dd=read.table(ff, header=T,sep=',', stringsAsFactors = F)
    vv.rm=grep('\\*',dd[,1])
    dd=dd[-vv.rm,]
    dd1=dd[,c('CDR3.pep.','V','C','copy')]
    dd1[,1]=paste('C',dd1[,1],'W', sep='')
    dd1[,2]=gsub('h','',dd1[,2])
    dd1[,3]=gsub('h','',dd1[,3])
    dd1=dd1[order(dd1[,4],decreasing=T),]
    Nm=min(10000, nrow(dd1))
    dd1=dd1[1:Nm,]
    dd1=cbind(dd1, ID=gsub('.csv','',ff0))
    write.table(dd1, file=paste(outdir, ff0, sep=''), quote=F, sep='\t', row.names=F, col.names=F)
  }
}


## Merged lung MDA, OV, Shipp and Control2020:
tmp.pred=tmp.pca[,1]+tmp.pca[,2]
tmp.cor=cor(t(tmp.dd), tmp.pred, method='s')
vv=which(abs(tmp.cor)>=0.3)
tmp.dd.vv=tmp.dd[vv,]
colnames(tmp.dd.vv)=tmp.labels
heatmap(tmp.dd.vv, col=redblue) -> tmp.heatmap

png(file='LungShippOV_merged.png',height=4000, width=8000, res=400)
par(mar=c(3,2,0,0),las=2)
image(t(tmp.dd.vv[tmp.heatmap$rowInd, tmp.heatmap$colInd]), col=redblue, axes=F)
axis(1, at=seq(0,1,length.out=468), tmp.labels[tmp.heatmap$colInd], cex.axis=0.2)
axis(2, at=seq(0,1,length.out=length(vv)), vv[tmp.heatmap$rowInd], cex.axis=0.2)
dev.off()

## Tina's Wish grant figure 4
vv.train.n=which(OVlabels[colnames(RFU.ov)] %in% c('Benign','Normal'))
vv.train.t=which(OVlabels[colnames(RFU.ov)] %in% c('highGrade'))
pList.train=c()
for(ii in 1:nrow(RFU.ov)){
  tmp=t.test(RFU.ov[ii, vv.train.n], RFU.ov[ii,vv.train.t])
  pList.train=c(pList.train, tmp$p.value)
}
names(pList.train)=1:5000

tmp.mm = names(sort(pList.train))[1:50]

tmp.dd.train=cbind(RFU.ov[,vv.train.n], RFU.ov[,vv.train.t])
tmp.lab.train=c(rep(0, length(vv.train.n)), rep(1, length(vv.train.t)))

svm.fit=svm(y=tmp.lab.train, x= t(tmp.dd.train[tmp.mm,]))
predict(svm.fit, t(tmp.dd.train[tmp.mm,])) -> val.train  ## training accuracy
roc(tmp.lab.train ~ val.train)

## encoding repertoire with variable gene
xx=EncodeRepertoire(f0, w=0)
dist(xx[,1:5000]) -> tmp

fit = cmdscale(sqrt(1-cor.mat), k=800)
x.mat.tmp=EncodeRepertoire(ff)
rownames(x.mat.tmp)=NULL
tmp=dist(x.mat.tmp)
tmp=as.matrix(tmp)
tmp.d=melt(tmp)
tmp.d=tmp.d[which(tmp.d[,3]>0),]
#vv.rm=which(abs(nchar(as.character(tmp.d[,1])) - nchar(as.character(tmp.d[,2])))>3)
#tmp.d=tmp.d[-vv.rm,]
tmp.yy=as.numeric(ddam[tmp.d[,1],5]==ddam[tmp.d[,2],5])
tmp.roc.v=roc(tmp.yy ~ tmp.d[,3])
tmp.roc.v

tmp.aucs=c()

## Test distance to center in kmeans with increasing number of centers
load('km1000.Rdata')
NpercList=c()
for(ii in 1:5000){
  cat('.')
tmp.dd=dd.control0[which(km.obj$cluster==ii),]
tmp.dd=rbind(tmp.dd, km.obj$centers[ii,])
Nd=nrow(tmp.dd)
tmp.dist=as.matrix(dist(tmp.dd))
#plot(density(tmp.dist[nrow(tmp.dd),]))
Nwithin=length(which(tmp.dist[Nd, 1:(Nd-1)]<=0.25))
Nperc = Nwithin/(Nd-1)
NpercList=c(NpercList, Nperc)
}
NdList=table(km.obj$cluster)
plot(NdList, NpercList)

NpL5000=NpercList

## Main Figure 2: use RFU 5000
require(ape)
ii=172
km.obj=km5000
tmp.dd=dd.control0[which(km.obj$cluster==ii),]
tmp.dd=rbind(tmp.dd, km.obj$centers[ii,])
Nd=nrow(tmp.dd)
rownames(tmp.dd)[Nd]='Centroid'
nj(as.matrix(dist(tmp.dd))) -> tmp.tree
root(tmp.tree, outgroup=Nd) -> tmp.tree.r
par(mar=c(0,0,0,0))
plot(tmp.tree.r, cex=1.5, font=2, tip.color=c(rep(1,Nd-1),2))
add.scale.bar(0.18,1, lwd=2, cex=1.5)

## Main Figure 4: RFU heatmap for Shipp data
vv=apply(RFU.shipp.f, 1, sd)
vvs=names(which(vv>=3))
tmp.RFU=RFU.shipp.f[vvs,]
tmp.shipp.col=rep(2, ncol(RFU.shipp.f))
tmp.shipp.col[grep('Donor', colnames(RFU.shipp.f))]=1
heatmap(tmp.RFU, col=redblue, ColSideColors=as.character(tmp.shipp.col))

## Main Figure 4: PCA plot
  ## lymphoma
tmp.pca.shipp=prcomp(t(RFU.shipp.f))$x
par(mar=c(5,5,1,1))
plot(tmp.pca.shipp[,1], tmp.pca.shipp[,2], col=tmp.shipp.col, pch=19, cex=1.2, xlab='PC1', ylab='PC2', cex.lab=1.5, cex.axis=1.5)
legend('bottomleft', legend=c('Control','Lymphoma'), pch=19, col=1:2, cex=1.5, bty='n')

vioplot(tmp.pca.shipp[,1] ~ tmp.shipp.col, names=c('Lymphoma','HC'), cex.names=1.5, cex.axis=1.5, col='gold', ylim=c(min(tmp.pca.shipp[,1]), max(tmp.pca.shipp[,1])*1.2))
vioplot(tmp.pca.shipp[,2] ~ tmp.shipp.col, names=c('Lymphoma','HC'), cex.names=1.5, cex.axis=1.5, col='gold', ylim=c(min(tmp.pca.shipp[,2]), max(tmp.pca.shipp[,2])*1.2))

tmp.xx=tmp.pca.shipp[,1] + tmp.pca.shipp[,2]
roc.shipp=roc(tmp.shipp.col ~ tmp.xx)

par(mar=c(4,4,1,1))
plot(roc.shipp, lwd=3, cex.lab=1.5, cex.axis=1.5)
legend('bottomright', legend='AUC = 93.2%', cex=1.5, bty='n')

  ## lung cancer
tmp.pca.lung=prcomp(t(cbind(RFU.covid.adaptive, RFU.lung.MDA)))$x
tmp.lung.col=c(rep(1, ncol(RFU.covid.adaptive)), rep(2, ncol(RFU.lung.MDA)))

par(mar=c(5,5,1,1))
plot(tmp.pca.lung[,1], tmp.pca.lung[,2], col=tmp.lung.col, pch=19, cex=1.2, xlab='PC1', ylab='PC2', cex.lab=1.5, cex.axis=1.5)
legend('bottomleft', legend=c('COVID-19','Lung Cancer'), pch=19,col=1:2, cex=1.5, bty='n')

vioplot(tmp.pca.lung[,1] ~ tmp.lung.col, names=c('Lung','HC'), cex.names=1.5, cex.axis=1.5, col='gold', ylim=c(min(tmp.pca.lung[,1]), max(tmp.pca.lung[,1])*1.2))
vioplot(tmp.pca.lung[,2] ~ tmp.lung.col, names=c('Lung','HC'), cex.names=1.5, cex.axis=1.5, col='gold', ylim=c(min(tmp.pca.lung[,2]), max(tmp.pca.lung[,2])*1.2))


tmp.xx=tmp.pca.lung[,2]
roc.lung=roc(tmp.lung.col ~ tmp.xx)

par(mar=c(4,4,1,1))
plot(roc.lung, lwd=3, cex.lab=1.5, cex.axis=1.5)
legend('bottomright', legend='AUC = 96.3%', cex=1.5, bty='n')


  ## ovarian cancer
vv.ov=names(which(OVlabels[colnames(RFU.ov)]!='lowGrade'))
tmp.pca.ov=prcomp(t(RFU.ov[,vv.ov[-37]]))$x  ## remove GO223_TCRB 
tmp.ov.col=rep(1, nrow(tmp.pca.ov))
tmp.ov.col[which(OVlabels[rownames(tmp.pca.ov)]=='highGrade')]=2
tmp.ov.col[which(OVlabels[rownames(tmp.pca.ov)]=='Benign')]='gray'

par(mar=c(5,5,1,1))
plot(tmp.pca.ov[,1], tmp.pca.ov[,2], col=tmp.ov.col,pch=19, cex=1.2, xlab='PC1', ylab='PC2', cex.lab=1.5, cex.axis=1.5)
legend('bottomleft',legend=c('Control','Benign','Ovarian Cancer'), pch=19, col=c('black','gray','red'), cex=1.5, bty='n')

tmp.ov.col.f=factor(tmp.ov.col, levels=c(1,'gray',2))

par(mar=c(5,5,1,1))
vioplot(tmp.pca.ov[,1] ~ tmp.ov.col.f, names=c('HC','Benign','HGSC'), cex.names=1.5, cex.axis=1.5, col='gold', ylim=c(min(tmp.pca.ov[,1]), max(tmp.pca.ov[,1])*1.2))

tmp.xx=tmp.pca.ov[,1] 
tmp.yy=tmp.ov.col
tmp.yy[which(tmp.yy=='gray')]=NA
roc.ov1=roc(tmp.yy ~ tmp.xx)

tmp.xx=tmp.pca.ov[,1] 
tmp.yy=tmp.ov.col
tmp.yy[which(tmp.yy=='gray')]=1
roc.ov2=roc(tmp.yy ~ tmp.xx)

par(mar=c(4,4,1,1))
plot(roc.ov1, lwd=3, cex.lab=1.5, cex.axis=1.5)
legend('bottomright', legend='AUC = 79.6%', cex=1.5, bty='n')

plot(roc.ov2, lwd=3, cex.lab=1.5, cex.axis=1.5)
legend('bottomright', legend='AUC = 77.9%', cex=1.5, bty='n')

## Main Figure 4: Volcano plot: ggrepel
tmp.p.list=c()
tmp.fc.list=c()
vv.shipp.cancer=grep('^ID', colnames(RFU.shipp.f))
vv.shipp.HC=grep('Donor', colnames(RFU.shipp.f))
for(ii in 1:5000){
  xx=RFU.shipp.f[ii, vv.shipp.cancer]
  yy=RFU.shipp.f[ii, vv.shipp.HC]
  tmp.pp=t.test(xx, yy)
  tmp.p.list=c(tmp.p.list, tmp.pp$p.value)
  tmp.fc.list=c(tmp.fc.list, mean(xx)/mean(yy))
}
names(tmp.p.list)=names(tmp.fc.list)=rownames(RFU.shipp.f)

tmp.dd.rare=data.frame(x=log(tmp.fc.list),y=-log10(p.adjust(tmp.p.list,'BH')),NAME=rownames(RFU.shipp.f))
tmp.dd.rare=na.omit(tmp.dd.rare)
tmp.dd.rare[which(!is.finite(tmp.dd.rare[,1])),1]=-4
g=ggplot(tmp.dd.rare,aes(x,y,label=NAME))
tmp.dd.rare.s=subset(tmp.dd.rare, y>=2)
tmp.dd.rare.s=cbind(tmp.dd.rare.s,COL=rep(2,nrow(tmp.dd.rare.s)))
tmp.dd.rare.s[which(tmp.dd.rare.s[,1]<0),4]=3
#tmp.to.label.u=c('ITGAE','GZMB','ZNF683','CCL3','CCL4','PRF1','LAG3','TIGIT','IFNG','PDCD1','CTLA4')
#tmp.to.label.d=c('KLF2','SELL','KLRG1')
#tmp.dd.rare.s[which(tmp.dd.rare.s[,'NAME']%in% tmp.to.label.u),4]=2
#tmp.dd.rare.s[which(tmp.dd.rare.s[,'NAME']%in% tmp.to.label.d),4]=3
tmp.rare.color=c(rgb(150,50,50,alpha=150,max=255),'darkred','darkgreen')
g+geom_point(size=1,alpha=0.1)+xlim(-1.5,2)+
  geom_label_repel(data=tmp.dd.rare.s,aes(x,y,label=NAME, max.overlaps=20),size=3,label.size=0.6,color=tmp.rare.color[tmp.dd.rare.s[,'COL']],segment.colour=rgb(50,50,50,alpha=70,max=255))+
  theme_bw()+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))+xlab('Mean Fold Change (log)')+ylab('Adjusted P value (-log10)')


## Figure S1
tmp.rocs=c(list(tmp.roc.100), list(tmp.roc.200), list(tmp.roc.300),list(tmp.roc.400), list(tmp.roc.500),
           list(tmp.roc.600), list(tmp.roc.700), list(tmp.roc.800),list(tmp.roc.900), list(tmp.roc.1000))

sens.b=0.2413676
tmp.sens=c()
for(rr in tmp.rocs){
  vv=min(which(rr$specificities>=0.9))
  tmp.sens=c(tmp.sens, rr$sensitivities[vv])
}
names(tmp.sens)=seq(100,1000,by=100)
par(mar=c(5,5,2,0), las=2)
barplot(tmp.sens, cex.lab=1.5, cex.axis=1.5, cex.names=1.5, col='cornflowerblue', ylim=c(0,0.25))
abline(h=sens.b, lwd=2, col='orange')

## Main Figure 3: TCGA data analysis
## create a control dataset from healthy donors
tmp.RFU.control = apply(RFU.batch2, 1, sum)
tmp.tt.tcga = table(RFU.tcgaNew[,7])
Nall=sum(tmp.RFU.control)
Ntcga=sum(tmp.tt.tcga)
tmp.pe.tcga=c()
tmp.ORe.tcga=c()

for(ii in names(tmp.tt.tcga)){
  pp0=tmp.RFU.control[ii]/Nall
  pp1=tmp.tt.tcga[ii]/Ntcga
  tmp.p=binom.test(tmp.tt.tcga[ii], Ntcga, pp0)
  tmp.pe.tcga=c(tmp.pe.tcga, tmp.p$p.value)
  tmp.ORe.tcga=c(tmp.ORe.tcga, pp1/pp0)
}

tmp.pe.tcga.p=p.adjust(tmp.pe.tcga,'BH')

#rr.selected.tcga= names(which(tmp.pe.tcga.p<=0.01 & tmp.ORe.tcga >= 2 & tmp.tt.tcga>=20))
rr.selected.tcga= names(which(tmp.pe.tcga.p<=0.05 & tmp.ORe.tcga >= 1))  ## n=252
rr.depleted.tcga= names(which(tmp.pe.tcga.p<=0.05 & tmp.ORe.tcga <= 0.5))  ## n=70 
## Identify cancer-type specific RFU groups
## use the more lenient criterion to select RFU
tmp.tt.cc=table(RFU.tcgaNew[,4], RFU.tcgaNew[,7])
#rr.selected.tcga= names(which(tmp.pe.tcga.p<=0.05 & tmp.ORe.tcga >= 2 & tmp.tt.tcga>=10))
tmp.tt.cancer=table(RFU.tcgaNew[,4])
ss.cancer=names(tmp.tt.cancer[which(tmp.tt.cancer>=1000)])
apply(t(tmp.tt.cc[ss.cancer, rr.selected.tcga]), 1, function(x)x/tmp.tt.cancer[ss.cancer]) -> tmp.tt.cc.f
tmp.tt.cc.denoise=tmp.tt.cc
tmp.tt.cc.denoise[which(tmp.tt.cc.denoise<=2)]=0
apply(t(tmp.tt.cc.denoise[ss.cancer, rr.selected.tcga]), 1, function(x)x/tmp.tt.cancer[ss.cancer]) -> tmp.tt.cc.denoise.f
heatmap(tmp.tt.cc.denoise.f[ss.cancer, rr.selected.tcga], col=redblue) -> tmp.heatmap
par(mar=c(3,4,0,0),las=2)
image(sqrt(t(tmp.tt.cc.denoise.f[ss.cancer, rr.selected.tcga][tmp.heatmap$rowInd, tmp.heatmap$colInd])), col=redblue, axes=F)
axis(2, at=seq(0,1,length.out=length(ss.cancer)), ss.cancer[tmp.heatmap$rowInd], cex.axis=1, font=2)
axis(1, at=seq(0,1,length.out=length(rr.selected.tcga)), rr.selected.tcga[tmp.heatmap$colInd], cex.axis=0.8)

## define cancer-type-specific RFUs
cancer.RFU.list=c()
for(cc in ss.cancer){
tmp.p.cc.list=c()
for(ii in 1:ncol(tmp.tt.cc.denoise.f)){
  vvc=grep(cc, ss.cancer)
  xx=tmp.tt.cc.denoise.f[-vvc,ii]-tmp.tt.cc.denoise.f[vvc,ii]
  if(sd(xx)<=0.000000001){
    tmp.p.cc.list=c(tmp.p.cc.list, 1)
    next
  }
  tmp.p=wilcox.test(xx, alternative='less')
  tmp.p.cc.list=c(tmp.p.cc.list, tmp.p$p.value)
}
names(tmp.p.cc.list)=colnames(tmp.tt.cc.denoise.f)
tmp.p.cc.list.a=p.adjust(tmp.p.cc.list,'BH')
tmp.ss=names(which(tmp.p.cc.list.a<=0.05))
cancer.RFU.list=c(cancer.RFU.list, list(tmp.ss))
}
names(cancer.RFU.list)=ss.cancer



ssAll=colnames(All.expr)
pDEGmat=c()
pFCmat=c()
for(ii in 1:length(rr.selected.tcga)){
  rr=rr.selected.tcga[ii]
  print(rr)
  tmp.p.DEG=c()
  tmp.FC=c()
  ss=RFU.tcga[which(RFU.tcga[,7]==rr),3]
  ss=intersect(ss, ssAll)
  ssn=ssAll[which(!ssAll %in% ss)]
  for(ii in 1:nrow(All.expr)){
    x1=na.omit(All.expr[ii,ss])
    x2=All.expr[ii,ssn]
    if(length(x1)<=10){
      tmp.p.DEG=c(tmp.p.DEG, NA)
      tmp.FC=c(tmp.FC, NA)
      next
    }
    tmp.p=wilcox.test(x1, x2)
    tmp.p.DEG=c(tmp.p.DEG, tmp.p$p.value)
    if(mean(x2, na.rm=T)==0)tmp.FC=c(tmp.FC,0) else tmp.FC=c(tmp.FC, mean(x1, na.rm=T)/mean(x2, na.rm=T))
  }
  tmp.p.DEG.a=p.adjust(tmp.p.DEG, 'BH')
  pDEGmat=cbind(pDEGmat, tmp.p.DEG.a)
  pFCmat=cbind(pFCmat, tmp.FC)
}
colnames(pDEGmat)=colnames(pFCmat)=rr.selected.tcga
rownames(pDEGmat)=rownames(pFCmat)=rownames(All.expr)

apply(pDEGmat, 1, function(x)length(which(x<=0.05))) -> tmp
vv.rm=which(tmp<=10)
pDEGmat.keep=pDEGmat[-vv.rm,]
pFCmat.keep=pFCmat[-vv.rm,]

## mutation correlation : NOT useful!
if(1){
mutation.pos.tt=table(TCGAmaf.f[,4])
selected.mutations=names(which(mutation.pos.tt>=20))
ss.mut=unique(TCGAmaf.f[,1])
Ns=length(unique(ss.mut))
Nm=length(selected.mutations)
MutMat=matrix(0, ncol=Nm, nrow=Ns)
rownames(MutMat)=ss.mut
colnames(MutMat)=selected.mutations
for(ii in selected.mutations){
  vv=grep(ii, TCGAmaf.f[,4])
  tmp.ss=unique(TCGAmaf.f[vv,1])
  MutMat[tmp.ss,ii]=1
}

ss.tcga=substr(unique(RFU.tcga[,3]), 1, 16)
ss.mut=substr(ss.mut, 1, 16)
rownames(MutMat)=ss.mut
ss.tcga=intersect(ss.tcga, ss.mut)

MutationTCRmatP=c()
MutationTCRmatOR=c()
for(rr in rr.selected.tcga){
  ss.rr=substr(unique(RFU.tcga[which(RFU.tcga[,7]==rr),3]), 1, 16)
  tmp.vec=rep(0,length(ss.tcga))
  names(tmp.vec)=ss.tcga
  ss.rr=intersect(ss.rr, ss.tcga)
  tmp.vec[ss.rr]=1
  tmp.OR=c()
  tmp.pp=c()
  if(length(unique(tmp.vec))==1){
    tmp.OR=c(tmp.OR, rep(NA, Nm))
    tmp.pp=c(tmp.pp, rep(NA, Nm))
    next
  }
  print(rr)
  for(ii in selected.mutations){
    tmp.yy=MutMat[ss.tcga,ii]
    if(length(unique(tmp.yy))==1){
      tmp.OR=c(tmp.OR, NA)
      tmp.pp=c(tmp.pp,NA)
      next
    }
    tmp.test=fisher.test(tmp.vec, tmp.yy)
    tmp.OR=c(tmp.OR, tmp.test$estimate)
    tmp.pp=c(tmp.pp, tmp.test$p.value)
  }
  MutationTCRmatP=cbind(MutationTCRmatP, tmp.pp)
  MutationTCRmatOR=cbind(MutationTCRmatOR, tmp.OR)
}
}

## Silhouette plot from randomly selected centers
Nclus=20
rand.clus=sample(1:5000, Nclus)
tmp.cluster=km5000$cluster
vv.clus=which(tmp.cluster %in% rand.clus)

tmp.cluster.s=tmp.cluster[vv.clus]
tmp.Dist=daisy(dd.control0[vv.clus,])

plot(silhouette(tmp.cluster.s, tmp.Dist), col=c('cyan','purple'))

## Using Shugay 2017 dataset to define age-associated RFUs as housekeeping RFUs: not sure if useful
shugay.age=shugay.info[colnames(RFU.shugay),4]
sh.age.p=c()
sh.age.cor=c()
for(ii in 1:5000){
  tmp=cor.test(shugay.age, RFU.shugay[ii,], method='s')
  sh.age.p=c(sh.age.p, tmp$p.value)
  sh.age.cor=c(sh.age.cor, tmp$estimate)
}
names(sh.age.p)=names(sh.age.cor)=rownames(RFU.shugay)
sh.age.pa=p.adjust(sh.age.p,'BH')

rr.housekeeping = names(which(sh.age.pa<=0.05 & sh.age.cor<=-0.2))

## Very interesting analysis: subgroups of cancer-associated RFUs
heatmap(cor(RFU.covid[rr.selected.tcga,], method='s'), keep.dendro=T)->tmp.heatmap  ## correlation structure preserved across multiple datasets
cut(tmp.heatmap$Rowv, 2.35)$lower -> tmp.subtrees

rr.selected.tcga.gr1=rr.selected.tcga[unlist(tmp.subtrees[[1]])]
rr.selected.tcga.gr2=rr.selected.tcga[unlist(tmp.subtrees[[2]])]
rr.selected.tcga.gr3=rr.selected.tcga[unlist(tmp.subtrees[[3]])]

rr.selected.tcga.groups=rep(1,length(rr.selected.tcga))
names(rr.selected.tcga.groups)=rr.selected.tcga
rr.selected.tcga.groups[rr.selected.tcga.gr1]=1  ## helper CD4
rr.selected.tcga.groups[rr.selected.tcga.gr2]=2  ## indolent or autoreactive
rr.selected.tcga.groups[rr.selected.tcga.gr3]=3  ## cytotoxic, exhausted


apply(tmp.RFU, 2, rank)-> tmp.rank

apply(tmp.rank[rr.selected.tcga.gr1,], 2, mean)-> tmp.rs1
apply(tmp.rank[rr.selected.tcga.gr2,], 2, mean)-> tmp.rs2
apply(tmp.rank[rr.selected.tcga.gr3,], 2, mean)-> tmp.rs3

apply(tmp.rank[rr.housekeeping,],2, mean) -> tmp.hks

## Using 10x single cell data to explore the cell phenotypes of selected RFUs
require(Seurat)
# Donor 1
#InputDir='/Users/bo/Desktop/UTSW/Research/CatchBond/data/singleCell/new20K/'
#donor1=Read10X(InputDir)
#d1=CreateSeuratObject(donor1, project='RFU')
  ## somehow Seurat doesn't work for this dataset. 
InputDir='/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/singleCell/DanaPeer/GEX/BC10/'
mtx=readMM(paste(InputDir,'GSM3148577_BC10_TUMOR1_matrix.mtx', sep=''))  
fm=read.table(paste(InputDir,'GSM3148577_BC10_TUMOR1_genes.tsv', sep=''),header=F,sep='\t')
mtx=as.matrix(mtx)
rownames(mtx)=fm[,2]
gm=read.table(paste(InputDir,'GSM3148577_BC10_TUMOR1_barcodes.tsv', sep=''),header=F,sep='\t')
colnames(mtx)=gm[,1]
sc.dd1=mtx
sc.tcr1=read.table('/Users/bo/Desktop/UTSW/Research/TCR_peptide/data/singleCell/DanaPeer/TCR/GSM3148582_BC10_TUMOR1_filtered_contig_annotations.csv',header=T,sep=',',stringsAsFactors=F)
vv=which(sc.tcr1[,'chain']=='TRB' & nchar(sc.tcr1[,'cdr3'])>=8 & sc.tcr1[,'umis']>=2)
sc.tcr1b=sc.tcr1[vv,]
tmp.tt=table(sc.tcr1b[,1])
vvn=names(which(tmp.tt>1))
sc.tcr1b=sc.tcr1b[which(!sc.tcr1b[,1] %in% vvn),]
sc.tcr1b.s=sc.tcr1b[,c('cdr3','v_gene','barcode')]
write.table(sc.tcr1b.s, file='donor1_tcr_forRFU.txt', quote=F, sep='\t', row.names=F, col.names=F)
TMP=AssignRFUs('donor1_tcr_forRFU.txt')

RFU.sc.d1=TMP$TCR

## link TCR data to RFU and GEX
sc.tcr1b.s=cbind(sc.tcr1b.s, RFU=RFU.sc.d1)

ss.barcode=intersect(sc.tcr1b.s[,3], colnames(sc.dd1))
sc.dd1=sc.dd1[,ss.barcode]

tmp.dd=sc.dd1[,sample(colnames(sc.dd1),10000)]  ## using all the cells will overwhelm the memory

apply(tmp.dd, 1, sd) -> tmp
gg.selected=names(which(tmp>=0.5))

## do a DEG analysis on cancer-associated RFUs vs others
bb.cancer=sc.tcr1b.s[which(sc.tcr1b.s[,4] %in% rr.selected.tcga.gr1),3]
vv1=which(colnames(sc.dd1) %in% bb.cancer)
tmp.pp=c()
tmp.fc=c()
for(ii in 1:nrow(sc.dd1)){
  xx1=sc.dd1[ii,vv1]
  xx2=sc.dd1[ii,-vv1]
  tmp=wilcox.test(xx1,xx2)
  tmp.pp=c(tmp.pp, tmp$p.value)
  fc=mean(xx1)/mean(xx2)
  tmp.fc=c(tmp.fc, fc)
}

## Main figure 5A: 
xx=tmp.RFU.control/Nall
yy=tmp.tt.tcga/Ntcga
xx=xx[names(yy)]
nnc=names(xx)
xx=as.numeric(xx)
yy=as.numeric(yy)
names(xx)=names(yy)=nnc
par(mar=c(5,5,1.5,1.5))
plot(xx,yy, cex=0.8, col=rgb(50,50,50,alpha=50,max=255), xlab='Expected RFU Frequency', ylab='TCGA RFU Frequency', cex.lab=1.5, cex.axis=1.5)
points(xx[rr.selected.tcga], yy[rr.selected.tcga], pch=1, col='purple',lwd=2)
points(xx[rr.depleted.tcga], yy[rr.depleted.tcga], pch=1, col='cornflowerblue',lwd=2)
abline(0,1,lwd=3, lty=2, col='gold')
legend('bottomright',legend=c('Cancer Enriched (n=252)','Cancer Depleted (n=70)'), pch=1, col=c('purple','cornflowerblue'), bty='n', cex=1.3)

  ## violin plot showing selected T cell biomarkers
tmp.ss.enriched=unique(RFU.tcga[which(RFU.tcga[,7] %in% rr.selected.tcga),3])
tmp.ss.other=unique(RFU.tcga[which(!RFU.tcga[,7] %in% rr.selected.tcga),3])
tmp.ss.enriched=intersect(tmp.ss.enriched, colnames(All.expr))
tmp.ss.other=intersect(tmp.ss.other, colnames(All.expr))

## to match purity distribution, design an accept/reject sampling procedure: importance
load('~/Desktop/DFCI/microEnv/scripts/AGPall.Rdata')
Purity.enriched=na.omit(AGP.all[substr(tmp.ss.enriched, 1, 12),3])
Purity.other=na.omit(AGP.all[substr(tmp.ss.other, 1, 12),3])
NPe=length(Purity.enriched)
NPo=length(Purity.other)
tmp.ss.other.s=c()
tmp.purity.other=c()
for(ss in tmp.ss.other){
  tmp.purity=AGP.all[substr(ss, 1, 12),3]
  if(is.na(tmp.purity))next
  tmp.qq=length(which(abs(Purity.enriched - tmp.purity)<=0.1))/NPe
  tmp.qqo=length(which(abs(Purity.other - tmp.purity)<0.1))/NPo
  Y=runif(1)
  if(Y<= tmp.qq/tmp.qqo){
    tmp.ss.other.s=c(tmp.ss.other.s, ss)
    tmp.purity.other=c(tmp.purity.other, tmp.purity)
  }
}

tmp.P.enriched=c()
tmp.FC.enriched=c()

for(ii in 1:nrow(All.expr)){
  xx=All.expr[ii, tmp.ss.enriched]
  yy=All.expr[ii, tmp.ss.other.s]
  xx=na.omit(xx)
  if(length(xx)==0){
    tmp.P.enriched=c(tmp.P.enriched, NA)
    tmp.FC.enriched=c(tmp.FC.enriched, NA)
    next    
  }
  if(sd(xx, na.rm=T)==0 | sd(yy, na.rm=T)==0){
    tmp.P.enriched=c(tmp.P.enriched, NA)
    tmp.FC.enriched=c(tmp.FC.enriched, NA)
    next
  }
  tmp=wilcox.test(xx,yy)
  tmp.P.enriched=c(tmp.P.enriched, tmp$p.value)
  tmp.FC.enriched=c(tmp.FC.enriched, mean(xx, na.rm=T)/mean(yy, na.rm=T))
}

tmp.P.enriched.a=p.adjust(tmp.P.enriched.a,'BH')
names(tmp.P.enriched.a)=names(tmp.FC.enriched)=rownames(All.expr)
plot(log2(tmp.FC.enriched),-log10(tmp.P.enriched.a), pch=19, col=rgb(50,50,50,alpha=50,max=255), cex=0.5)

par(mar=c(1,3,2,0), mfrow=c(1,4))
for(gg in c('TCF7','SIT1','CD96','SIRPG')){
  vioplot(All.expr[gg,tmp.ss.enriched], All.expr[gg,tmp.ss.other.s],names=NA, xaxt='n', cex.axis=1.5, main=gg, ylim=c(0,max(All.expr[gg,], na.rm=T)*1.05),col=c('orange','darkgreen'), cex.lab=1.5, cex.main=1.5)
}

## Main figure 5
tmp.mat=-log10(pDEGmat.keep.na)
tmp.tt=table(RFU.tcga[,7])
tmp.tt=sort(tmp.tt[rr.selected.tcga],decreasing=T)
tmp.selected.ggs=c('PDCD1','LAG3','GZMB','IFNG','CD3E','CRTAM','ZNF683','TOX','TCF7','CD96','CXCR3','TIGIT','CD27','IL2RA','ICOS','SIRPG','SLA2','CD2','CD28','HAVCR2','CTLA4','FOXP3')
tmp.selected.ags=rownames(tmp.mat)[grep('^OR|^KR|TTT',rownames(tmp.mat))]
tmp.mat.s=tmp.mat[tmp.selected.ggs, names(tmp.tt[1:30])]
#tmp.mat.sa=tmp.mat[tmp.selected.ags, names(tmp.tt[1:30])]
heatmap(tmp.mat.s) -> tmp.hg
#heatmap(tmp.mat.sa) -> tmp.ag

par(mar=c(3.1,3,0,0),las=2, mgp=c(3,0.1,0))
image(tmp.mat.s[tmp.hg$rowInd, tmp.hg$colInd], axes=F)
axis(1, at=seq(0,1,length.out=length(tmp.selected.ggs)), tmp.selected.ggs[tmp.hg$rowInd], cex.axis=0.9, tck=0, col = NA, col.ticks = NA)
axis(2, at=seq(0,1,length.out=length(colnames(tmp.mat.s))), colnames(tmp.mat.s)[tmp.hg$colInd], cex.axis=0.9,tck=0, font=2, col = NA, col.ticks = NA)

par(mar=c(0,0,3,4),las=2, mgp=c(3,0.1,0))
image(t(tmp.mat.sa[tmp.ag$rowInd, tmp.hg$colInd]), axes=F, col=redblue)
axis(4, at=seq(0,1,length.out=length(tmp.selected.ags)), tmp.selected.ags[tmp.ag$rowInd], cex.axis=0.9, col = NA, col.ticks = NA)
axis(3, at=seq(0,1,length.out=length(colnames(tmp.mat.s))), colnames(tmp.mat.s)[tmp.hg$colInd], cex.axis=1.1, font=2, col = NA, col.ticks = NA)



fan.RFU <- function(x, centroid=NA, is.plot=TRUE, MaxR=NA, is.grid=TRUE){
  ## x is a n-by-m encoding matrix, with each row the coordinates of a sequence in an RFU
  require(plotrix)
  if(is.na(centroid)){
    centroid=colMeans(x)
  }
  nx=nrow(x)
  x1=rbind(x, centroid)
  x.dist=as.matrix(dist(x1))
  vv=order(x.dist[nx+1,1:nx])
  pos.angles=seq(0,360,length.out=nx+1)[1:nx]
  pos.radius=c(x.dist[nx+1,vv[1]])
  for(ii in 2:nx){
    new.pos=x.dist[nx+1, vv[ii]] + x.dist[vv[ii-1], vv[ii]]
    #new.pos=x.dist[nx+1, vv[ii]]
    pos.radius=c(pos.radius, new.pos)
  }
  if(is.na(MaxR))MaxR=max(pos.radius)
  if(is.plot){
    if(is.grid)
    polar.plot(pos.radius, pos.angles, radial.lim=c(0,MaxR), rp.type='polygon',
               poly.col=rgb(50,20,120,alpha=150,max=255),lwd=3, mar=c(1,1,2,1), 
               grid.col='lightgreen', rad.col='lightblue',line.col=rgb(50,50,50,alpha=150,max=255), 
               show.grid.labels=F,show.grid=is.grid,show.radial.grid=is.grid)
    else     polar.plot(pos.radius, pos.angles, radial.lim=c(0,MaxR), rp.type='polygon',
                        poly.col=rgb(50,20,120,alpha=150,max=255),lwd=3, mar=c(0,0,0,0), 
                        grid.col='lightgreen', rad.col='lightblue',line.col=rgb(50,50,50,alpha=150,max=255), 
                        show.grid.labels=F,show.grid=is.grid,show.radial.grid=is.grid)
  }
  SU=cbind(Radius=pos.radius, Angles=pos.angles)
  return(SU)
}


plot.fan.RFU <- function(RFUnumber, Data=dd.control0, Km=km5000, is.plot=TRUE, is.grid=TRUE){
  vv=which(Km$cluster==RFUnumber)
  x=Data[vv,]
  SU=fan.RFU(x, is.plot=is.plot, is.grid=is.grid)
  return(SU)
}


AWCs=c()
for(ii in 1:5000){
  SU=plot.fan.RFU(ii, is.plot=F)
  AWCs=c(AWCs, mean(SU[,1]))
}


## Figure 5: Putting together another visualization for RFU: sequence Logo plot by ClustalOmega
require(msa)
library(ggseqlogo)

RFUnum=1199
SS=names(which(km5000$cluster==RFUnum))
tmp.aln=msa(SS, type='protein', method='ClustalOmega')
tmp.SS=msaConvert(tmp.aln)

p1=ggseqlogo(tmp.SS$seq)
p1 + theme(axis.text.x = element_text(size=20))

## Figure 3: Cancer type enrichment for selected RFUs  ## not significant
if(1){
COLs= c("dodgerblue2","gold", # red
        "green4",
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "cyan",
        "skyblue2","#FB9A99", # lt pink
        "palegreen2",
        "gray70", "khaki2",
        "maroon","orchid1","deeppink1","blue1","steelblue4",
        "darkturquoise","green1","yellow4",
        "darkorange4","brown")

names(COLs)[1:14]=ss.cancer
RFUnum=4739

tt.RFU0=table(RFU.tcga[,4])
pp.RFU0=tt.RFU0[ss.cancer]/sum(tt.RFU0[ss.cancer])

tt.RFU=table(RFU.tcga[which(RFU.tcga[,7]==RFUnum),4])
Nr=sum(tt.RFU)

pp.RFU=rep(1,length(ss.cancer))
names(pp.RFU)=ss.cancer
OR.RFU=pp.RFU

for(cc in ss.cancer){
  if(!cc %in% names(tt.RFU)){
    cc.count=0
  }else cc.count=tt.RFU[cc]
  tmp=binom.test(cc.count, Nr, pp.RFU0[cc])
  pp.RFU[cc]=tmp$p.value
  OR.RFU[cc]=cc.count/Nr/pp.RFU0[cc]
}

pp.RFU.a=p.adjust(pp.RFU,'BH')

par(mar=c(1,2.5,1,0))
barplot(OR.RFU, ylim=c(0,max(OR.RFU)*1.05), col=COLs[1:14],names=NA, cex.axis=1.5, border=NA)
abline(h=1,lwd=2, col='darkred', lty=2)

vvs=which(pp.RFU.a<=0.1)  ## set FDR level 0.1

start.pos=seq(0,1.4*14, by=1.4)
xxs=start.pos[vvs]+0.6
}

## Figure 5: single cell
#RFU.zemin = AssignRFUs(ff='/Users/bo/Desktop/UTSW/Research/CatchBond/data/singleCell/Zemin/TCR_data_selected_forGIANA.txt')
load('scAnalysisWorkplace0924.Rdata')
par(mar=c(0,0,0,0))
plot(sc.tsne$Y, col=rgb(50,50,50,alpha=50, max=255), pch=19, cex=0.6)
ss=RFU.zemin.dd[which(RFU.zemin.dd[,6]==1199),3]
points(sc.tsne$Y[ss,], col='firebrick', pch=19, cex=0.9)
ss=RFU.zemin.dd[which(RFU.zemin.dd[,6]==4739),3]
points(sc.tsne$Y[ss,], col='orange', pch=19, cex=0.9)
ss=RFU.zemin.dd[which(RFU.zemin.dd[,6]==4033),3]
points(sc.tsne$Y[ss,], col='purple', pch=19, cex=0.9)
ss=RFU.zemin.dd[which(RFU.zemin.dd[,6]==1909),3]
points(sc.tsne$Y[ss,], col='cornflowerblue', pch=19, cex=0.9)



## DEG analysis of single cells within a RFU group
Allss=colnames(scData.s)
GGtop=names(sort(tmp.sd, decreasing=T))[1:3000]

tmpPmat=c()
tmpFCmat=c()

for(ii in c(1199, 4739, 4033, 1909)){

ss=RFU.zemin.dd.s[which(RFU.zemin.dd.s[,6]==ii),3]
vss=Allss[which(!Allss %in% ss)]
tmp.P.list=c()
tmp.FC.list=c()
for(gg in GGtop){
  xx=scData.s[gg,ss]
  yy=scData.s[gg,vss]
  tmp=wilcox.test(xx,yy)
  tmp.P.list=c(tmp.P.list, tmp$p.value)
  if(mean(yy, na.rm=T)==0)tmp.fc=NA else tmp.fc = mean(xx, na.rm=T)/mean(yy, na.rm=T)
  tmp.FC.list=c(tmp.FC.list, tmp.fc)
}
names(tmp.P.list)=names(tmp.FC.list)=GGtop

tmp.P.list=p.adjust(tmp.P.list)

tmpPmat=cbind(tmpPmat, tmp.P.list)
tmpFCmat=cbind(tmpFCmat, tmp.FC.list)
}

colnames(tmpPmat) = colnames(tmpFCmat) = c(1199, 4739, 4033, 1909) 


## volcano plot for 1199 and 4739
ii='4739'

gg.tcga=names(which(pDEGmat[,ii] <= 0.01 & pFCmat[,ii]>=1))

tmp.p.list=tmpPmat[,ii]
tmp.fc.list=tmpFCmat[,ii]

tmp.dd.rare=data.frame(x=log(tmp.fc.list),y=-log10(tmp.p.list),NAME=GGtop)
tmp.dd.rare=na.omit(tmp.dd.rare)
tmp.dd.rare[which(!is.finite(tmp.dd.rare[,1])),1]=-4
g=ggplot(tmp.dd.rare,aes(x,y,label=NAME))
tmp.dd.rare.s=subset(tmp.dd.rare, y>=2)
tmp.dd.rare.s=cbind(tmp.dd.rare.s,COL=rep('gray25',nrow(tmp.dd.rare.s)))
tmp.ss=intersect(gg.tcga, rownames(tmp.dd.rare.s))
tmp.dd.rare.s[which(tmp.dd.rare.s$NAME %in% tmp.ss),4]='red'
#tmp.to.label.u=c('ITGAE','GZMB','ZNF683','CCL3','CCL4','PRF1','LAG3','TIGIT','IFNG','PDCD1','CTLA4')
#tmp.to.label.d=c('KLF2','SELL','KLRG1')
#tmp.dd.rare.s[which(tmp.dd.rare.s[,'NAME']%in% tmp.to.label.u),4]=2
#tmp.dd.rare.s[which(tmp.dd.rare.s[,'NAME']%in% tmp.to.label.d),4]=3
tmp.rare.color=c(rgb(150,50,50,alpha=150,max=255),'darkred','darkgreen')
g+geom_point(size=1,alpha=0.1)+xlim(0.5,2)+ylim(0,20)+
  geom_label_repel(data=tmp.dd.rare.s,aes(x,y,label=NAME, max.overlaps=100),size=4,force=10, label.size=1,color=tmp.dd.rare.s[,'COL'],segment.colour=rgb(50,50,50,alpha=70,max=255))+
  theme_bw()+theme(axis.text=element_text(size=14),axis.title=element_text(size=14))+xlab('Mean Fold Change (log)')+ylab('Adjusted P value (-log10)')

## Lung cancer Dejima
#RFU.lung.dejima = dd.RFU
#RFU.lung.dejima.f=dd.RFU[,which(tmp.N>=2000)]
EarlyLung.labels=rep('ADC', ncol(RFU.lung.Dejima.f))
tmp.ss=colnames(RFU.lung.Dejima.f)
names(EarlyLung.labels)=tmp.ss
vv.N=grep('N',tmp.ss)
EarlyLung.labels[vv.N]='Normal'
vv.H=grep('H',tmp.ss)
EarlyLung.labels[vv.H]='AAH'
vv.S=grep('S',tmp.ss)
EarlyLung.labels[vv.S]='AIS'
vv.M=grep('M',tmp.ss)
EarlyLung.labels[vv.M]='MIA'


vv=which(EarlyLung.labels %in% c('Normal','ADC'))
lung.d1=RFU.lung.Dejima.f[,vv]
tmp.col.ld1=as.numeric(as.factor(EarlyLung.labels[vv]))

lung.d2=cbind(RFU.lung.Creelan, RFU.control.inhouse)  ## control.inhouse: Carter blood 

tmp.col.ld2=c(rep(1,20), rep(2,34))

lung.d3=cbind(RFU.lung.cascone[,grep('blood', colnames(RFU.lung.cascone))], RFU.control.HC)  ## control.HC: Adaptive control 2020
tmp.col.ld3=c(rep(1,15),rep(2, 88))

lung.test = cbind(RFU.lung.MDA, RFU.covid.adaptive)
tmp.col.test=c(rep(1,121), rep(2, 160))

p1.list = c()
for(ii in 1:nrow(lung.d1)){
  tmp=wilcox.test(lung.d1[ii,] ~ EarlyLung.labels[vv])
  p1.list=c(p1.list, tmp$p.value)
}
fc1.list=aggregate(t(lung.d1), list(EarlyLung.labels[vv]), mean)

p2.list = c()

for(ii in 1:nrow(lung.d2)){
  tmp=wilcox.test(lung.d2[ii,] ~ tmp.col.ld2)
  p2.list=c(p2.list, tmp$p.value)
}
fc2.list=aggregate(t(lung.d2), list(tmp.col.ld2), mean)

p3.list = c()

for(ii in 1:nrow(lung.d3)){
  tmp=wilcox.test(lung.d3[ii,] ~ tmp.col.ld3)
  p3.list=c(p3.list, tmp$p.value)
}
fc3.list=aggregate(t(lung.d3), list(tmp.col.ld3), mean)


names(p1.list)=names(p2.list)=names(p3.list)=rownames(lung.d1)

tmp.fc3=fc3.list[1,tmp.ss]/fc3.list[2,tmp.ss]
tmp.fc2=fc2.list[1,tmp.ss]/fc2.list[2,tmp.ss]

knn.pred <- function(x, y, n=5){
  ## x: distance matrix between reference and query samples, row: ref, col: query
  ## y: binary labels of each reference class
  y.classes=unique(y)
  cMat=c()
  for(cc in y.classes){
    vvc=which(y==cc)
    yc=apply(x[vvc,], 2, function(x)median(sort(x)[1:n]))
    cMat=cbind(cMat, yc)
  }
  colnames(cMat)=y.classes
  cMat
}

tmp.cor=cor(cbind(lung.d2, lung.d3), lung.test)
tmp.cols=c(tmp.col.ld2, tmp.col.ld3)
tmp.cols=2-tmp.cols  ## make control = 0, cancer = 1

knn.pred(sqrt(1-tmp.cor), tmp.cols, n=3) -> cMat
tmp.rr=cMat[,1]/cMat[,2]
tmp.yy=c(rep(1,121), rep(0,88))
plot(roc(tmp.yy ~ tmp.rr))  ## AUC=0.926, SP=0.95, SN=0.758; SP=0.99, SN=0.458; SP=1, SN=0.45

  ## lung additional validation
RFU.lung.Federic.Nf=RFU.lung.Federic.N[,which(tmp.N[vv]>=2000)]



## For ovarian cancer, KNN based prediction
ref.d1=cbind(RFU.inhouse.ov.rcc)
ref.d2=cbind(RFU.inhouse.pda,RFU.covid.adaptive)
tmp.ref=cbind(ref.d1, ref.d2)
tmp.yy = c(rep(1,ncol(ref.d1)), rep(0,ncol(ref.d2)))

dis.d1=RFU.ov.hgsc
#dis.d2=cbind(RFU.control.inhouse)
dis.d2=RFU.ov.benign
tmp.dis = cbind(dis.d1, dis.d2)

tmp.cor=cor(tmp.ref, tmp.dis)  ## Jaya cohort
tmp.yy.dis = c(rep(1, ncol(dis.d1)), rep(0, ncol(dis.d2)))

knn.pred(sqrt(1-tmp.cor), tmp.yy, n=3) -> tmp.cMat
tmp.rr=tmp.cMat[,1]/tmp.cMat[,2]
boxplot(tmp.rr ~ tmp.yy.dis)  ## 
tmp.roc = roc(tmp.yy.dis ~ tmp.rr)  ## AUC=0.819, SP: 96%, SN: 53%

## Benign sample selection criteria for SHARE biobank
dd=read.csv('Benign Ovarian Cases-Bo.txt',header=T,sep='\t')
vv=grep('tumor', dd[,'Histologic.Type...First.letter.Capitalize.'])
ddv=dd[-vv,]
newdd=ddv[which(as.numeric(ddv[,9])>=50 & as.numeric(ddv[,9])<=80 & ddv[,'Race']=='White' & ddv[,'Reviewed..Yes..Y...Pending.']!='Y')[1:42],]
write.table(newdd, file='SHAREbiobank_benign_ovarian_selected_BL.txt', quote=F, sep='\t', row.names=F)



## Multi-cancer separation
RFU.melanoma.riaz.pre.f=RFU.melanoma.riaz.pre[,which(tmp.N[colnames(RFU.melanoma.riaz.pre)]>=1000)]
  ## pairwise
  # lung vs melanoma
tmp.dd=cbind(RFU.lung.MDA, RFU.melanoma.riaz.pre.f)
tmp.lab=c(rep(1,121),rep(2,26))
tmp.pca=prcomp(t(tmp.dd[rr.hk,]))$x  ## check batch effect
boxplot(tmp.pca[,1] ~ tmp.lab)
  # ov vs melanoma
tmp.dd=cbind(RFU.ov.hgsc, RFU.melanoma.riaz.pre.f)
tmp.pca=prcomp(t(tmp.dd[rr.hk,]))$x
tmp.lab=c(rep(1,45), rep(2,26))
boxplot(tmp.pca[,1] ~ tmp.lab)
  # ov vs lung
tmp.dd=cbind(RFU.ov.hgsc, RFU.lung.MDA)
tmp.pca=prcomp(t(tmp.dd[rr.hk,]))$x
tmp.lab=c(rep(1,45), rep(2,121))
boxplot(tmp.pca[,1] ~ tmp.lab)
  # ov vs 



## CD4 vs CD8 RFUs
tmp.p.CD48.shipp=c()
tmp.fc.CD48.shipp=c()
for(ii in 1:5000){
  xx1=mean(RFU.shipp.CD48.f[ii,which(tmp.lab.cd4==0)],na.rm=T)
  xx2=mean(RFU.shipp.CD48.f[ii,which(tmp.lab.cd4==1)],na.rm=T)
  tmp=wilcox.test(RFU.shipp.CD48.f[ii,] ~ tmp.lab.cd4)
  tmp.p.CD48.shipp=c(tmp.p.CD48.shipp,tmp$p.value)
  tmp.fc.CD48.shipp=c(tmp.fc.CD48.shipp, xx2/xx1)
}

tmp.p.CD48.ms=c()
tmp.fc.CD48.ms=c()
for(ii in 1:5000){
  xx1=mean(RFU.MS.CD48[ii,which(tmp.ms.cd4==0)],na.rm=T)
  xx2=mean(RFU.MS.CD48[ii,which(tmp.ms.cd4==1)],na.rm=T)
  tmp=wilcox.test(RFU.MS.CD48[ii,] ~ tmp.ms.cd4)
  tmp.p.CD48.ms=c(tmp.p.CD48.ms,tmp$p.value)
  tmp.fc.CD48.ms=c(tmp.fc.CD48.ms, xx2/xx1)
}

rr.CD4=rownames(RFU.MS.CD48)[which(tmp.fc.CD48.shipp> 2 & tmp.fc.CD48.ms > 2)]
rr.CD8=rownames(RFU.MS.CD48)[which(tmp.fc.CD48.shipp< 0.5 & tmp.fc.CD48.ms < 0.5)]

tmp.pca=prcomp(t(RFU.shipp.CD48.f))$x
par(mar=c(5,5,1,1))
plot(tmp.pca[,1], tmp.pca[,2], col=tmp.lab.cd4+1, pch=19, xlab='PC1', ylab='PC2', cex.axis=1.5, cex.lab=1.5)
legend('bottomleft', legend=c('CD4', 'CD8'), pch=19, cex=1.5, bty='n', col=c(1,2))

tmp.xx = tmp.pca[,1] + tmp.pca[,2]
tmp.roc=roc(tmp.lab.cd4 ~ tmp.xx)
par(mar=c(4,4,1,1))
plot(tmp.roc, lwd=3, cex.lab=1.5, cex.axis=1.5)
legend('bottomright', legend='AUC = 94.1%', cex=1.5, bty='n')

tmp.pca=prcomp(t(RFU.MS.CD48))$x
par(mar=c(5,5,1,1))
plot(tmp.pca[,1], tmp.pca[,2], col=tmp.ms.cd4+1, pch=19, xlab='PC1', ylab='PC2', cex.axis=1.5, cex.lab=1.5)
legend('bottomleft', legend=c('CD4', 'CD8'), pch=19, cex=1.5, bty='n', col=c(1,2))

tmp.xx = tmp.pca[,1] + tmp.pca[,2]
tmp.roc=roc(tmp.ms.cd4 ~ tmp.xx)
par(mar=c(4,4,1,1))
plot(tmp.roc, lwd=3, cex.lab=1.5, cex.axis=1.5)
legend('bottomright', legend='AUC = 93.8%', cex=1.5, bty='n')


## Twin data (Cystic Fibrosis)
## Genetic factors contribute to repertoire similarity, which decreases with age, due to more diversified environmental exposure
tmp.cor=cor(RFU.cf.twin)
tmp=gsub('[_AB]{1,2}','',colnames(RFU.cf.twin))

tmp.ii=seq(1,32,by=2)
intratwin.cor=diag(tmp.cor[tmp.ii, tmp.ii+1])
twin.age=twin.info[tmp[tmp.ii],2]
cor.test(twin.age, intratwin.cor, method='s')  # p value = 0.044, cor = -0.51









