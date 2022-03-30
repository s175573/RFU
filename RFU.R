## Repertoire Functional Units
load('RFU_housekeepking.Rdata')  ## housekeeping RFUs for batch detection
load('RFUcancer_list0901.Rdata') ## cancer related RFUs selected by TCGA data
f='CDRtable.txt'
x=read.table(f, header=F, sep='\t', stringsAsFactors = F)
vv.rm=c(grep('[*.]',x[,2]), grep('[*.]',x[,3]))
x=x[-vv.rm,]
CDR1=x[,2]
names(CDR1)=x[,1]
CDR2=x[,3]
names(CDR2)=x[,1]

ProcessAdaptiveVgenes <- function(dd){
  Vgene_o=c(1,2,9,13:19,26:28,30)
  Vgene_o=as.character(Vgene_o)
  dd=dd[which(nchar(as.character(dd[,2]))>0),]
  gsub('TCRBV[0]{0,1}','TRBV',dd[,2])->tmpV
  ## Multiple calls
  vv.m=grep('/',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],'/')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  vv.m=grep(',',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],',')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  tmpV=gsub('-0','-',tmpV)
  v_digit=grep('\\*',tmpV)
  if(length(v_digit)>0)tmpV[- v_digit] = paste(tmpV[- v_digit],'*01',sep='') else tmpV=paste(tmpV,'*01',sep='')
  ## 1. Orphan V genes do not have "-1", need to remove
  Vnumbers=gsub('TRBV','',tmpV)
  Vnumbers=gsub('\\*.+','',Vnumbers)
  Vnumbers=gsub('-.+','',Vnumbers)
  vv.o=which(Vnumbers %in% Vgene_o)
  tmpV1=tmpV
  tmpV1[vv.o]=gsub('-1','',tmpV1[vv.o])
  ## 2. Non-orphan V genes but without "-1", need to add
  vv.no=which(! Vnumbers %in% Vgene_o)
  vv.non=grep('-',tmpV1[vv.no])
  if(length(vv.non)>0)tmpV1[vv.no][-vv.non]=gsub('\\*01','-1*01',tmpV1[vv.no][-vv.non])
  dd[,2]=tmpV1
  return(dd)
}

PrepareAdaptiveFile <- function(indir,outdir,AlternativeFormat=FALSE,thr=10000, naive=FALSE, pattern='.tsv'){
  ffs=dir(indir,full.names=T)
  ffs=ffs[grep(pattern, ffs)]
  for(ff in ffs){
    if(length(grep('\\.tsv|\\.txt',ff))==0 ) next
    ff0=unlist(strsplit(ff,'\\/'))
    ff0=ff0[length(ff0)]
    if(!file.exists(outdir))dir.create(outdir)
    newff=paste(outdir,'TestReal-',ff0,sep='')
    if(exists(newff))next
    print(ff)
    ddnew=read.table(ff,header=T,sep='\t',stringsAsFactors=F)
    if(AlternativeFormat){
      ii.aa=grep('amino_acid', colnames(ddnew))[1]
      ii.vv=grep('v_gene', colnames(ddnew))[1]
      ii.ff=grep('max_productive_frequency',colnames(ddnew))
    }else{
      ii.aa=2
      ii.vv=6
      ii.ff=4
    }
    tmp.vv.rm=grep('\\*',ddnew[,ii.aa])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.vv.rm=grep('X',ddnew[,ii.aa])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.nn=nchar(ddnew[,ii.aa])
    tmp.vv=which(tmp.nn>=10 & tmp.nn<=24)
    #tmp.vv=which(tmp.nn>=1)
    ddnew=ddnew[tmp.vv,]
    ddnew=ddnew[which(ddnew[,ii.vv]!='unresolved'),]
    ddnew=ddnew[grep('^C.+[FV]$',ddnew[,ii.aa]),]
    #ss=sum(ddnew[,3])
    #ddnew[,4]=ddnew[,3]/ss
    if(naive)ddnew=ddnew[order(ddnew[,ii.ff]),] else ddnew=ddnew[order(ddnew[,ii.ff],decreasing=T),]
    #q75=quantile(ddnew[,3],0.75)
    #THR=min(ddnew[,3]+1)
    #tmp.vv=which(ddnew[,3]>max(q75,THR))
    #if(length(tmp.vv)<=thr)tmp.vv=1:thr
    if(nrow(ddnew)>thr)ddnew=ddnew[1:thr,c(ii.aa,ii.vv,ii.ff)] else ddnew=ddnew[,c(ii.aa,ii.vv,ii.ff)]
    #ddnew=ddnew[1:5000,]
    ddnew=ddnew[1:min(thr,nrow(ddnew)),]
    ddnew=ProcessAdaptiveVgenes(ddnew)
    ddnew=cbind(ddnew, RANK=rank(ddnew[,3])/nrow(ddnew))
    TMP=gsub('.+/TCR_peptide/data/Adaptive/','',indir)
    tmp=strsplit(TMP,'/')[[1]][[1]]
    ddnew=cbind(ddnew, Info=paste(tmp,ff0,sep=':'))
    write.table(ddnew,file=newff,quote=F,sep='\t',row.names=F)
  }
}

getTrimers <- function(cdr3, st=3, ed=3){
  nL=nchar(cdr3)
  posL=st:(nL-ed)
  ss=sapply(posL, function(x)substr(cdr3, x, x+2))
  ss=unlist(ss)
  return(ss)
}

getTrimerPairs <- function(x){
  n=nrow(x)
  Tpairs=c()
  for(ii in 1:n){
    s1=getTrimers(x[ii,1])
    for(jj in ii:n){
      if(jj==ii)next
      s2=getTrimers(x[jj,1])
      vv=which(s1 != s2)
      Tpairs=rbind(Tpairs, cbind(s1[vv], s2[vv]))
    }
  }
  return(Tpairs)
}

getTrimerMatrix <- function(f, thr.c=6){
  ## f is the input TCR clustering data, desirably with more than 10M sequences
  ## tri-mer matrix
  print('Caution! This function is a one-time training implementation which takes long time to finish.')
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
  
  dd=read.table(f,header=FALSE,sep='\t',stringsAsFactors=F)
  
  ## select small-sized clusters, to exclude small-world effect
  tt=table(dd[,2])
  nn=names(which(tt<= thr.c))
  dd0=dd[which(dd[,2] %in% nn),]
  
  ## remove clusters with identical CDR3s
  tmp0=split(dd0, dd0[,2])
  sapply(tmp0, function(x)length(unique(x[,1])))->LL
  vvL=which(LL>1)
  tmp0=tmp0[vvL]
  
  TpairMat=c()
  count=0
  for(x in tmp0){
    if(count %% 10000 == 0)print(count)
    Tp=getTrimerPairs(x)
    TpairMat=rbind(TpairMat, Tp)
    count=count+1
  }
  
  TT= table(TpairMat[,1], TpairMat[,2])
  for(ss in rownames(TT)){
    trimerMat[ss,colnames(TT)]=TT[ss,]
  }
  return(trimerMat)
}

getTrimerEncoding <- function(x, DIM=100){
  ## Input: trimer-replacement matrix
  x0 <- x + t(x)
  ss <- apply(x0, 1, sum)
  vv.rm <- which(ss==0)
  x0 <- x0[-vv.rm, -vv.rm]
  x0c <- cor(x0, method='s')
  fit <- cmdscale(sqrt(1-x0c), k=DIM)
  fit
}

#EncodeRepertoire <- function(ff, Vgene=TRUE, w=1){
  dd1=read.table(ff, header=F,sep='\t',stringsAsFactors = F)
  vv=grep('^C',dd1[,1])
  dd1=dd1[vv,]
  if(Vgene){
    Vs=which(dd1[,2] %in% names(CDR1))
    dd1=dd1[Vs,]
  }
  if(is.null(ncol(dd1)))ssL=sapply(dd1, getTrimers) else ssL=sapply(dd1[,1], getTrimers)
  foo <- function(ss){
    ss=gsub('.+_','',ss)
    ss=intersect(ss, rownames(fit))
    if(length(ss)==1)return(fit[ss,])
    xs=t(fit[ss,])
    xx = rowMeans(xs)    
  }
  x.mat=sapply(ssL, foo)
  if(Vgene){
    Vss=dd1[,2]
    cdr1=CDR1[Vss]
    cdr2=CDR2[Vss]
    ssV1=lapply(cdr1, getTrimers, st=1, ed=0)
    ssV2=lapply(cdr2, getTrimers, st=1, ed=0)
    x.mat1=sapply(ssV1, foo)
    x.mat2=sapply(ssV2, foo)
    x.matV=rbind(x.mat1, x.mat2)
    colnames(x.matV)=colnames(x.mat)
    x.mat0= rbind(x.mat, w*x.matV)
  }else x.mat0=x.mat
  return(t(x.mat0))
}

EncodeRepertoire<-function(ff){
  dd1=read.table(ff, header=F,sep='\t',quote='',stringsAsFactors = F)
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
  list(RFU= RFU, N=Nmiss, COR=tmp, TCR=tmpv)
}

#AssignRFUs.d <- function(ff, CL=km5000, THR=0.25, BIG=FALSE){
  ## Assign repertoire functional units defined by k-means cluster centroids from pooled control samples
  require(Rfast)
  require(parallel)
  dd=EncodeRepertoire(ff)
  ddc=CL$centers
  Nd=nrow(dd)
  Nc=nrow(ddc)
  
  if(BIG){
    cat(" Switching to big data mode: slow but memory friendly")
    foo <- function(x){
      DIFF=apply(ddc, 1, function(y)y-x)
      apply(DIFF, 2, function(x)sum(x^2)) -> tmp.dist
      vvx=which.min(tmp.dist)
      minx=min(tmp.dist)
      return(c(minx, vvx))
    }
    TMP = mclapply(1:Nd, function(x)foo(dd[x,]))
    TMPm=do.call(cbind, TMP)
    tmp= TMPm[1,]
    tmpv = as.integer(TMPm[2,])
  }else{ 
    tmp.dist = Dist(rbind(dd,ddc))
    tmp.dist=as.matrix(tmp.dist)
    tmp.mat=tmp.dist[1:Nd, (Nd+1):(Nd+Nc)]
    tmp=apply(tmp.mat, 1, min)
    tmpv=apply(tmp.mat, 1, function(x)which.min(x)[1])
  }
  vv=which(tmp>THR)
  tmpv[vv]=NA
  
  tt=table(tmpv)
  tt0=table(CL$cl)
  RFU = rep(0, length(tt0))
  names(RFU)=1:length(RFU)
  RFU[names(tt)] = tt
  #RFU = RFU / as.numeric(tt0)  ## normalize by original cluster sizes
  RFU = RFU / nrow(dd) * 10000 ## normalize by total TCR in the data
  Pmiss= length(vv)/nrow(dd)
  list(RFU= RFU, N=nrow(dd),P=Pmiss, COR=tmp, TCR=tmpv)
}

RFUbatch <- function(DIR){
  dd.RFU =c()
  tmp.N=c()
  ffs=dir(DIR, full.names=T)
  ffs=ffs[grep('[tc]sv|txt',ffs)]
  for(ff in ffs[1:length(ffs)]){
    print(ff)
    if(length(grep('[tc]sv|txt',ff))==0)next
    RFU = AssignRFUs(ff)
    dd.RFU=cbind(dd.RFU, RFU$RFU)
    tmp.N=c(tmp.N, RFU$N)
  }
  
  ffs=gsub('.+/','',ffs)
  
  colnames(dd.RFU)=gsub('.+TestReal-','',ffs)
  colnames(dd.RFU)=gsub('.[tc]sv|txt','',colnames(dd.RFU))
  
  TMP=c(RFU=list(dd.RFU), Nmiss=list(tmp.N))
  return(TMP)
}

batch.check <- function(d1, d2, mode=1, TYPE='median'){
  # d1, d2 are RFU matrices
  if(mode==1){
    if(TYPE=='median'){
      m1=apply(d1, 1, median)
      m2=apply(d2, 1, median)
    }
    if(TYPE=='mean'){
      m1=apply(d1, 1, mean)
      m2=apply(d2, 1, mean)      
    }
    r1=rank(m1)
    r2=rank(m2)
  }
  if(mode==2){
    m1=apply(d1, 2, rank)
    m2=apply(d2, 2, rank)
    if(TYPE=='median'){
      r1=apply(m1, 1, median)
      r2=apply(m2, 1, median)
    }
    if(TYPE=='mean'){
      r1=apply(m1, 1, mean)
      r2=apply(m2, 1, mean)      
    }
  }
  par(mar=c(5,4,0,0))
  plot(r1, r2, pch=19, xlab='Rank of data1', ylab='Rank of data2', cex.lab=1.5, cex.axis=1.5)
  cor(r1, r2, method='s')
}

batch.check2 <- function(d1, d2){
  # d1, d2 are RFU matrices
  n1 = ncol(d1)
  n2 = ncol(d2)
  labels = c(rep(1, n1), rep(2, n2))
  dd= cbind(d1, d2)
  dd = dd[rr.hk,]
  pca0 = prcomp(t(dd))
  pca = pca0$x
  par(mar=c(4,5,2,0), mfrow=c(2,2))
  tmp=wilcox.test(pca[,1] ~ labels)
  boxplot(pca[,1] ~ labels, main=paste('PC1: P=', round(tmp$p.value,3),sep=''), cex.axis=1.3, cex.lab=1.3, ylab='PC1', xlab='')
  tmp=wilcox.test(pca[,2] ~ labels)
  boxplot(pca[,2] ~ labels, main=paste('PC2: P=', round(tmp$p.value,3),sep=''), cex.axis=1.3, cex.lab=1.3, ylab='PC2', xlab='')
  tmp=wilcox.test(pca[,3] ~ labels)
  boxplot(pca[,3] ~ labels, main=paste('PC3: P=', round(tmp$p.value,3),sep=''), cex.axis=1.3, cex.lab=1.3, ylab='PC3', xlab='')
  tmp=wilcox.test(pca[,4] ~ labels)
  boxplot(pca[,4] ~ labels, main=paste('PC4: P=', round(tmp$p.value,3),sep=''), cex.axis=1.3, cex.lab=1.3, ylab='PC4', xlab='')
}


DiffDensityPlot <- function(f1, f2, TMP=NULL, n1=-1, ZLIM=0.0002){
  require(MASS)
  require(Rtsne)
  if(!is.null(TMP)){
    ## just make the plot with given Rtsne object
    n2=nrow(TMP$Y)-n1
    kde2d(TMP$Y[1:n1,1], TMP$Y[1:n1,2],n=400)->dens1
    kde2d(TMP$Y[(n1+1):(n1+n2),1], TMP$Y[(n1+1):(n1+n2),2],n=400, lims = c(range(TMP$Y[1:n1,1]), range(TMP$Y[1:n1,2])))->dens2
    image(dens2$x, dens2$y, dens2$z - dens1$z, col=redblue, zlim=c(-ZLIM,ZLIM))
    return()
  }
  d1=EncodeRepertoire(f1)
  d2=EncodeRepertoire(f2)
  d1=na.omit(d1)
  d2=na.omit(d2)
  TMP0=Rtsne(rbind(d1, d2), check_duplicate=FALSE)
  n1=nrow(d1)
  print(n1)
  n2=nrow(d2)
  kde2d(TMP0$Y[1:n1,1], TMP0$Y[1:n1,2],n=400)->dens1
  kde2d(TMP0$Y[(n1+1):(n1+n2),1], TMP0$Y[(n1+1):(n1+n2),2],n=400, lims = c(range(TMP0$Y[1:n1,1]), range(TMP0$Y[1:n1,2])))->dens2
  image(dens2$x, dens2$y, dens2$z - dens1$z, col=redblue, zlim=c(-ZLIM,ZLIM))
  diffDensObj=c(tsne=list(TMP0), dens1=list(dens1), dens2=list(dens2), n1=n1)
  return(diffDensObj)
}

DiffRFUAnalysis <- function(RFU, labels, TEST='w'){
  classes = unique(labels)
  Nc=length(classes)
  Nr=nrow(RFU)
  DRA.mat=c()
  for(ii in 1:Nc){
    vv1=which(labels==classes[ii])
    for(jj in ii:Nc){
      if(jj==ii)next
      vv2=which(labels==classes[jj])
      tmp.mat=c()
      for(kk in 1:Nr){
        x1=RFU[kk, vv1]
        x2=RFU[kk, vv2]
        if(TEST=='w'){
          tmp=wilcox.test(x1, x2)
        }
        if(TEST=='t'){
          tmp=t.test(x1, x2)
        }
        tmp.fc=mean(x1)/mean(x2)
        tmp.mat=rbind(tmp.mat, c(Pval=tmp$p.value, FC=tmp.fc))
      }
      tmp.mat=cbind(tmp.mat, FDR=p.adjust(tmp.mat[,1]))
      Ci=classes[ii]
      Cj=classes[jj]
      Cij=paste(Ci, Cj, sep='_')
      colnames(tmp.mat)=paste(Cij,':', colnames(tmp.mat), sep='')
      DRA.mat=cbind(DRA.mat, tmp.mat)
    }
  }
  rownames(DRA.mat)=rownames(RFU)
  DRA.mat
}

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

ssRSEA <- function(x, labels){
  ## single sample RFU-set-enrichment analysis
  ## x: input RFU vector for a given sample
  ## label: binary vector with 1 being in the RFU set of interest
  vv=order(x,decreasing=T)
  xs=x[vv]
  labels=labels[vv]
  yy=c()
  y.cur=0
  n0=length(which(labels==0))
  S0=sum(xs[which(labels==1)])
  # for(ii in 1:length(xs)){
  #   if(labels[ii]==0){
  #     y.cur=y.cur-1/n0
  #   }else{
  #     y.cur=y.cur+xs[ii]/S0
  #   }
  #   yy=c(yy,y.cur)
  # }
  yList=rep(-1/n0, length(xs))
  yList[which(labels==1)]=xs[which(labels==1)]/S0
  yy=cumsum(yList)
  ddy=cbind(curve=yy,scores=xs,labels)
  rownames(ddy)=names(labels)
  return(ddy)   ## mean(ddy[,1]) is a predictor
}



