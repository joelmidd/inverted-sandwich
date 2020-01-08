rm(list=ls())
try(source("C:\\Users\\joelm\\Dropbox\\R.Functions\\FUNCTIONS.R"))
try(source("C:\\Users\\joel.middleton\\Dropbox\\R.Functions\\FUNCTIONS.R"))
try(source("C:\\Users\\M\\Dropbox\\R.Functions\\FUNCTIONS.R"))

#THIS WORKS IN R STUDIO TO GET THE LOCAL DIRECTORY WHERE FILE IS SAVED:
flnm.path<-get.filename(set.dir=T, extension=F)
setwd(flnm.path$path)

sim.var.ests<-function(y0, y1, x.vals, treat, strat, perms=NULL, dmat=NULL, p=NULL, recurs=6, interact=T){
  require(MASS)
  
  n<-length(y0)
  if(length(j.dim(x.vals))==1){
    k<-1
  }else{
    k<-ncol(x.vals)
  }
  
  xx<-gimmie.xx(x.vals
                , intercept=T, contrast=c(-1,1)
                , interact=interact, block=F)
  
  xx.sep<-gimmie.xx(x.vals
                    , intercept=F, contrast=c(-1,1)
                    , interact=interact, block=T)
  
  xx.mats<-list(xx=xx, xx.sep=xx.sep)
  yy<-c(-y0, y1)
  
  if(is.null(perms)){
    require(ri)
    perms<-genperms(Z=treat, blockvar = strat)
  }
  pperms<-rbind(perms==0, perms==1)
  
  if(is.null(dmat)){
    res.dmat<-get.d.mat(treat=treat, strat=strat)
    dmat<-res.dmat$dmat
    p <- res.dmat$p
  }
  d01<-gimmie.quad(dmat, 1,2)
  d.t <- bdiag(gimmie.quad(dmat,1,1)-(d01+t(d01))/2
               , gimmie.quad(dmat,2,2)-(d01+t(d01))/2)
  pmat<-(p%*%t(p))*(dmat+1)
  
  pmat.tilde<-pmat
  pmat.tilde[pmat==0]<- -99999999999
  sv.pmats<-list(pmat=pmat, pmat.tilde=pmat.tilde)
  sv.dmats<-list(d=dmat, d.t=d.t, d.t.wt=(d.t/pmat.tilde)/n^2)
  #d.t.wt.pi.pi<-(d.t/pmat.tilde)*(p%*%t(p))
  
  #########
  print("GET ALL BETAS AND PREDICTIONS")
  #########
  res.all.regs<-get.all.regs(pperms, xx, yy, hc.types=c("HC0", "HC2", "HC3"))
  
  #############
  print("GET MSE ESTIMATOR MATRIX, M")
  #############
  sv.mmats<-list()
  cb.sep<-as.vector(colSums(xx.sep)/n)
  cb<-as.vector(colSums(xx)/n)
  m.res<-lapply.reorg(as.data.frame(pperms), function(pp) 
    get.yMy(pp, x=xx, yy, cb, as.matrix(Diagonal(2*n))))
  
  rm.mres<-rowMeans(m.res$m)
  M<-m.res$m%*%t(m.res$m)/ncol(m.res$m) -rm.mres%*%t(rm.mres)
  sv.mmats$M<-M
   
  if(FALSE){ 
    sv.mmats$M.hc0<-Reduce('+', with(m.res, lapply(c(1:ncol(m.res$m)), function(i)
      matrix(i_x.xx.inv.x[,i], ncol=2*n)%*%Diagonal(x=m[,i]^2)%*%matrix(i_x.xx.inv.x[,i], ncol=2*n)
    )
    )
    )/ncol(m.res$m)
    sv.mmats$M.hc0.wt<-sv.mmats$M.hc0/pmat.tilde
  }
 
  
  if(FALSE){ 
    ########
    #Matrixes
    ########
    res.bread.sand<-m.res$xx.inv.x[c(1:(2*n))*3-1,i]
    #with(res.noStrata,  
    #                    sapply(c(1:ncol(pperms)), function(i) {
    #                      -m.res$xx.inv.x[c(1:(2*n))*3-1,i], n*2)
    
    sv.VvMats<-list()
    sv.VvMats$HC0 <- j.cov.n(t(res.bread.sand))  
    #for(i in 1:ncol(res.noStrata$pperms)){
    #  tmp<-(tmp+(res.bread.sand[,i]%*%t(res.bread.sand[,i]))/i)
    #}
    tmp<-diag(rowMeans(res.bread.sand))
    sv.mmats$HC0.resid<-tmp*pmat
    sv.mmats$HC0.resid.wt<-tmp
  
    res.bread.sand.HC2<-m.res$leverage*res.bread.sand
    sv.VvMats$HC2<-j.cov.n(t(res.bread.sand.HC2))  
    #for(i in 1:ncol(res.noStrata$pperms)){
    #  tmp<-(tmp+(res.bread.sand[,i]%*%t(res.bread.sand[,i]))/i)
    #}
    tmp<-diag(rowMeans(res.bread.sand.HC2))
    sv.mmats$HC2.resid<-tmp*pmat
    sv.mmats$HC2.resid.wt<-tmp
  }
  
  #######
  print("GET A BOUND FOR M")
  #######
  m01.m10<-(gimmie.quad(M, 1,2)+gimmie.quad(M, 2,1))/2
  M.t.init<- as.matrix(bdiag(gimmie.quad(M,1,1)-m01.m10
                        , gimmie.quad(M,2,2)-m01.m10))
  ones.mat<-matrix(rep(1, n^2), ncol=n)
  force.ind<-1-bdiag(ones.mat, ones.mat)
  tmat <- as.matrix(M.t.init-M)
  force.val<-force.ind*tmat
  eig.t<-eigen(as.matrix(tmat))
  for(i in 1:10){
    #print("eig.t!2")
    if(min(eig.t$values)< -.00001){
      update.res<-j.t.update(tmat, force.indicator = force.ind, force.value= force.val)
      tmat<-update.res$t
      eig.t<-update.res$eig.res
    }
  }
  sv.mmats$M.t  <- M+tmat
  sv.mmats$M.t.wt  <- (M+tmat)/pmat.tilde
  tmat.AS<-    -(as.matrix(sv.dmats$d)==-1)*M + diag(colSums((as.matrix(sv.dmats$d)==-1)*abs(M)))
  sv.mmats$M.tAS<- M +tmat.AS
  sv.mmats$M.tAS.wt<- (M +tmat.AS)/pmat.tilde
  
  
  ##########
  print("UNBIAS M FOR USE WITH U-HAT (RESIDUALS)") 
  ##########
  
  res.Mt.resid<-get.unbias.uMu(perm.mat = pperms, xx, yy, sv.mmats$M.t
                             , sv.pmats$pmat.tilde, recurs)
  sv.mmats$M.t.resid<-(res.Mt.resid$M.star.mat)
  sv.mmats$M.t.resid.wt<-(res.Mt.resid$M.star.mat/pmat.tilde)

  res.Mt.resid.AS<-get.unbias.uMu(perm.mat = pperms, xx, yy, sv.mmats$M.tAS
                             , sv.pmats$pmat.tilde, recurs)
  sv.mmats$M.tAS.resid<-(res.Mt.resid.AS$M.star.mat)
  sv.mmats$M.tAS.resid.wt<-(res.Mt.resid.AS$M.star.mat/pmat.tilde)
  
  ##HC0 INVERTED (M0)
  res.i_x.d.i_x<-sapply(c(1:ncol(pperms)), function(i) {
    get.i_x.d.i_x(matrix( m.res$xx.inv.x[,i]  , ncol=n*2)
                  , matrix(m.res$i_x.xx.inv.x[,i]  , ncol=n*2)
                  , p=rep(1,2*n)
                  , dmat=diag(2*n)
                  , n=1)
  }
  )
  
  #get matrix for M0
  sv.mmats$M0 <- matrix(rowMeans(res.i_x.d.i_x), ncol=n*2)
  sv.mmats$M0.wt<- sv.mmats$M0/pmat.tilde
  sv.mmats$M0.resid <-get.unbias.uMu(pperms
                                 , x=xx
                                 , y=yy, Mm=sv.mmats$M0
                                 , pmat.tilde=pmat.tilde
                                 , recurs=recurs)$M.star.mat
  sv.mmats$M0.resid.wt<-sv.mmats$M0.resid /pmat.tilde
  
  ##HC2 INVERTED (M2)
  
  res.i_x.d.i_x_M2<-sapply(c(1:ncol(pperms)), function(i) {
    get.i_x.d.i_x(matrix( m.res$xx.inv.x[,i]  , ncol=n*2)
                  , matrix(m.res$i_x.xx.inv.x[,i]  , ncol=n*2)
                  , p=rep(1,2*n)
                  , dmat=diag(2*n)
                  , n=1
                  , lev=m.res$leverage[,i])
  }
  )
  sv.mmats$M2 <- matrix(rowMeans(res.i_x.d.i_x_M2), ncol=n*2)
  sv.mmats$M2.wt<- sv.mmats$M2/pmat.tilde
  sv.mmats$M2.resid <-get.unbias.uMu(pperms
                                     , x=xx
                                     , y=yy, Mm=sv.mmats$M2
                                     , pmat.tilde=pmat.tilde
                                     , recurs=recurs)$M.star.mat
  sv.mmats$M2.resid.wt<-sv.mmats$M2.resid /pmat.tilde
  

  ##########
  print("GET THE VARIANCE ESTIMATORS FOR EACH RANDOMIZATION")
  ##########
  
  res.vhats<-as.data.frame(sapply(c(sv.mmats[c("M.t.wt", "M.tAS.wt", "M0.wt")])
                                      , function(mat) 
                                        get.all.uMu.uGiven(res.all.regs$y.mat*pperms
                                                           , as.matrix(mat)) 
  )
  )
  res.vhats<-cbind(res.vhats, 
                   as.data.frame(
                     sapply(res.all.regs$modelinfo[grep("HC",names(res.all.regs$modelinfo))]
                            , function(m) m[2,])))
  
  
  res.vhats<-cbind(res.vhats, 
                   as.data.frame(sapply(c(sv.mmats[grep("resid.wt", names(sv.mmats))])
                                  , function(mat) 
                                    get.all.uMu.uGiven(res.all.regs$u.mat*pperms
                                                       , as.matrix(mat)) 
  )
  ))
  

  res.vhats<- cbind(res.vhats,HC0=sapply(c(1:ncol(pperms)), function(i) get.MEHW(u=res.all.regs$u.mat[,i]*pperms[,i]
                                                                   , bread=matrix(m.res$xx.inv.x[,i], nrow=ncol(xx))
                                                                   , p=p
                                                                   , sv.dmats$d.t.wt
                                                                   , n
  )[2]))
  res.vhats<- cbind(res.vhats,
                    HC2=sapply(c(1:ncol(pperms)), function(i) get.MEHW(u=res.all.regs$u.mat[,i]*pperms[,i]
                                                                   , bread=matrix(m.res$xx.inv.x[,i], nrow=ncol(xx))
                                                                   , p=p
                                                                   , sv.dmats$d.t.wt
                                                                   , n
                                                                   , lev=(m.res$leverage[,i])
                    )[2]))
  
  
  
  ####################
  #
  ####################
  #estimated vars
  funcs<- c(var=j.var.n, mean=j.mean, min=min, max=max, pcnt.neg=function(b){mean(b<0)}, num.gr.0=function(b){sum(b>0)})
  summary.all   <-t(sapply(as.data.frame(res.vhats), function(ests) sapply(funcs, function(ff) ff(ests))))
  #summary.uMu   <-t(sapply(as.data.frame(res.vhats.uMu), function(ests) sapply(funcs, function(ff) ff(ests))))
  #summary.yMy   <-t(sapply(as.data.frame(res.vhats.yMy), function(ests) sapply(funcs, function(ff) ff(ests))))
  #summary.HC    <-t(sapply(as.data.frame(res.vhats.HC), function(ests) sapply(funcs, function(ff) ff(ests))))
  #summary.MEHW    <-t(sapply(as.data.frame(gEHW<-as.data.frame(cbind(gEHW0=res.vhats.MEHW0, gEHW2=res.vhats.MEHW2
  #                                                                   ,gEHW.avg= (res.vhats.MEHW0+res.vhats.MEHW2)/2)
  #)), function(ests) sapply(funcs, function(ff) ff(ests))))
  summary.b<- j.mat.summary(t(res.all.regs$modelinfo$betas))
  trueV=rep(j.var.n(res.all.regs$modelinfo$betas[cb==1,]), length(res.all.regs$modelinfo$betas[cb==1,]))
  
  list(summary= rbind(true=c(NA, trueV[1],rep(NA, ncol(summary.all)-2))
                     
       , summary=summary.all)
       , summary.b=summary.b
       , all.regs = res.all.regs 
       , b.ests= res.all.regs$modelinfo$betas
       , v.hats=res.vhats
       , sv.mats= list(mmats=sv.mmats
                       #, mHCmats=sv.mHCmats
                       , dmats=sv.dmats, pmats=sv.pmats
                       #, umats=sv.umats
       )
       , pperms=pperms
       , data.2n= list(xx.mats= xx.mats, yy= yy)
       , M.compute.details=m.res
       #, bread.sandwichs=list(HC0=res.bread.sand
        #                      , HC2=res.bread.sand)
       #, VvMats=sv.VvMats
  )
}

#####
#set sim parameters
#############

gen.sample<-function(sm){
  n<<-11
  k<<-1
  
  dgp.coef<<-rep(.7, k)
  set.seed(11111*sm)
  x.vals.all<<-matrix(rnorm(n*k), ncol=k)
  y0<<-sqrt(12)*(runif(n)-.5)*.7*abs((x.vals.all[,1]-min(x.vals.all[,1])))+x.vals.all%*%dgp.coef
  y1<<- y0+.5
  
  ordr<<-order(x.vals.all[,1])
  x.vals<<-x.vals.all[ordr,2:k]
  x.vals.all<<-x.vals.all[ordr,]
  y0<<-y0[ordr]
  y1<<-y1[ordr]
  
  summary(lm(y0~x.vals))
  strat<<-(c(1:n)>6)+1
  treat<<- rep(0, n) 
  #strat.samp.sizes<<- c(2,2)
  tb.strat<<-table(strat)
  for(str.nm in names(tb.strat)){
    treat[sample(c(1:n)[strat==as.numeric(str.nm)]
                 , tb.strat[str.nm]/3)
          ]<<-1
  }
  #x.vals<<-cbind(x.vals, strat)
  cbind(y0, y1, x.vals, strat, treat)
}

funcs<- c(var=j.var.n, mean=j.mean
          , min=min, max=max)

gen.sample.nostrat<-function(sm){
  n<<-9
  k<<-1
  
  dgp.coef<<-rep(.7, k)
  set.seed(11111*sm)
  x.vals<<-matrix(rnorm(n*k), ncol=k)
  y0<<-sqrt(12)*(runif(n)-.5)*.7*abs((x.vals[,1]-min(x.vals[,1])))+x.vals%*%dgp.coef
  y1<<- y0+.5
  ordr<<-order(x.vals[,1])
  x.vals<<-x.vals[ordr,1:k]
  y0<<-y0[ordr]
  y1<<-y1[ordr]
  print(summary(lm(y0~x.vals)))
  treat<<- rep(0, n) 
  treat[1:(round(n/2.5+.4999))]<<-1
  cbind(y0, y1, x.vals, treat)
}

#############
#run analysis
#############

##################
####no strata
#################
dat<-gen.sample.nostrat(3)
cbind(dat, x.vals)
summary(lm(dat[,"y0"]~x.vals))$r.squared
summary(lm(dat[,"y1"]~x.vals))$r.squared
apply(dat, 2, sd)

#debugonce(sim.var.ests)
res.noStrata<-sim.var.ests(y0, y1, x.vals, treat, strat=rep(1,n)
                           , recurs=12, interact=F)

res.noStrata$summary






#############
#############










dig<-c(0,3,3,3,3,0,0)
library(xtable)
xtable(print(res.noStrata$summary), digits=dig)





##################
#sims with strata 
##################

dat<-gen.sample(3)
res.noInt<-sim.var.ests(y0, y1, x.vals, treat, strat
                        , recurs=8, interact=F)

res.noInt.strat<-sim.var.ests(y0, y1, cbind(x.vals, strat), treat, strat
                              , recurs=8, interact=F)

res.noInt.stratOnly<-sim.var.ests(y0, y1, strat, treat, strat
                                  , recurs=8, interact=F)

dig<-c(0,3,3,3,3,0,0)

#xtable(res.noInt.strat.x1$summary)

round(res.noInt$summary,3)
res.noInt$summary

library(xtable)
xtable(res.noInt$summary)
tab.noInt<-cbind(y0=y0, y1=y1, x.vals=x.vals, treat=treat, strat)
colnames(tab.noInt)<-c("y0", "y1","x.vals", "treat", "strat")
xtable(tab.noInt)

###########################
#x interacts with treatment
############################
res<-sim.var.ests(y0, y1, x.vals, treat, strat
                  , recurs=8, interact=T)

#options(scipen=999)
round(res$summary,3)

t.for.hists<-cbind(M.t.resid.wt=res$t.ratos.MeanCentered$uMu[,"M.t.resid.wt"]
                   , M.t.wt=res$t.ratos.MeanCentered$yMy[,"M.t.wt"]
                   , res$t.ratos.MeanCentered$HC[,c(2:3)])
t.for.hists[t.for.hists< -5]<- -5
t.for.hists[t.for.hists>5]<- 5
j.hists(t.for.hists, tdf=n-2-k*2)




