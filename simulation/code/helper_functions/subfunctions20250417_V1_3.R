rbern = function(n,prob)
{
  x = runif(n,min=0,max=1)
  x.bern = ifelse(x<=prob,1,0)
  return(x.bern)
}
get.y = function(betas, design.x, ICCy, M.psu, psu.size){
  dim(betas) = c(length(betas), 1)   #make it matrix of ncol=1
  odds = exp( design.x %*% betas)
  N = nrow(design.x)
  p = odds/(1+odds) #Pr(y=1|x)
  if (ICCy==0) {
    y = ifelse(runif(N)<=p,1,0)
  }
  if (ICCy!=0){
    ei0 = rep(rnorm(M.psu),psu.size)
    eij = rnorm(N)
    Uij = rbern(N,sqrt(ICCy))
    thetaij = qnorm(p)
    threshold = Uij*ei0 + (1-Uij)*eij
    y = ifelse(threshold <= thetaij,1,0)
  }
  return(y)
}
#################################################################################
#-FUNCTION sam.pps is to select a sample with pps sampling                      #
# INPUT:  popul - the population including response and covariates              #
#         MSize - the size of measure to compute the selection probabilities(>0)#
#         n -     the sample size                                               #
# OUTPUT: sample selected with pps sampling of size n,                          #
#         including response, covariates, and sample weights                    #
#################################################################################
sam.pps=function(popul,Msize, n){
  prob.s = Msize/sum(Msize)
  N=nrow(popul)
  pps.samID=sample(N,n,replace=F,prob=Msize)   #F but make sure N is large!
  if (dim(popul)[2] == 1){
    sam.pps.data = as.data.frame(popul[pps.samID,])
    names(sam.pps.data) = names(popul)
  }else{sam.pps.data = popul[pps.samID,]}
  sam.pps.data$wt = sum(Msize)/n/Msize[pps.samID]               
  return(sam.pps = sam.pps.data)
}
#################################################################################
samp.slct = function(seed, fnt.pop, n, psu.name=NULL, m.psu=NULL, dsgn, size = NULL, size.I = NULL){
  #fnt.pop = pop; seed = seed.sim[1]; psu.name="psu"; 
  #dsgn = "pps"; 
  #size = "odds.s"; 
  #dsgn = "srs-pps"; 
  #dsgn = "pps-pps"; 
  #size.I = aggregate(fnt.pop[,size], list(fnt.pop[,psu.name]), sum)[,2]
  set.seed(seed)
  N = nrow(fnt.pop)
  # one-ste sample design
  if(dsgn=="pps"){
    samp = sam.pps(fnt.pop,fnt.pop[,size], n)
  }
  # two-stage cluster sampling design (informative design at the second stage)
  if(dsgn == "srs-pps"){
    #-- first stage: select clusters by srs
    size.psu = as.numeric(table(fnt.pop[,psu.name]))
    M.psu = length(size.psu)
    # m.psu = 25
    # n = 1000
    index.psuI = sample(1:M.psu, m.psu, replace = F)
    sample.I = fnt.pop[fnt.pop[,psu.name] %in% index.psuI,]
    sample.I$wt.I = M.psu/m.psu
    #-- second stage: select subj.samp within selected psus by pps 
    # calcualte the size for the second stage
    samp=NULL
    for (i in 1: m.psu){
      popn.psu.i= sample.I[sample.I[,psu.name]==index.psuI[i],]
      samp.i = sam.pps(
        popn.psu.i, 
        popn.psu.i[,size], 
        min(n / m.psu, nrow(popn.psu.i)))
      samp.i$wt = samp.i$wt*samp.i$wt.I
      samp = rbind(samp,samp.i)
    }#sum(samp.cntl$wt);nrow(fnt.cntl)
  }
  if(dsgn == "pps-pps"){
    #-- first stage: select clusters by pps
    size.psu = as.numeric(table(fnt.pop[,psu.name]))
    M.psu = length(size.psu)
    # Clt.samp = 25
    # n = 1000
    index.psuI = sam.pps(matrix(1:M.psu,,1),size.I, m.psu)
    index.psuI = index.psuI[order(index.psuI[,1]),]  #sort selected psus
    sample.I = fnt.pop[fnt.pop[,psu.name]%in% index.psuI[,1],]
    sample.I = sample.I[order(sample.I[,psu.name]),]
    sample.I$wt.I = rep(index.psuI[,'wt'], size.psu[index.psuI[,1]])
    #-- second stage: select subj.samp within selected psus by pps 
    # calcualte the size for the second stage
    samp=NULL
    for (i in 1: m.psu){
      popn.psu.i = sample.I[sample.I$psu==index.psuI[i,1],    ]
      samp.i = sam.pps(popn.psu.i,popn.psu.i[,size], n/m.psu)
      samp.i$wt = samp.i$wt*samp.i$wt.I
      samp = rbind(samp,samp.i)
    }#sum(samp.cntl$wt);nrow(fnt.cntl)
  }
  rownames(samp) = as.character(1:dim(samp)[1])
  return(samp)      
}



greg.f = function(samp, wt0, N.hat, aux.mtx=NULL, aux.tot, f_w=T){
  n=length(samp[,wt0]);n
  if(!is.null(aux.mtx)) samp = cbind(samp, aux.mtx)
  ds.wt = svydesign(ids=~1, data=samp, weights=as.formula(paste0("~",wt0)))
  calib.fm = as.formula(paste0("~", paste0(names(aux.tot)[-1],collapse="+")))
  V.hat = c(sum(samp[,wt0])/N.hat, svytotal(calib.fm,ds.wt))
  V.hat[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]=V.hat[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]/aux.tot[1]
  V = aux.tot
  V[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]=V[!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]/aux.tot[1]
  v.mtx = model.matrix(calib.fm,samp)
  v.mtx[,1] = aux.mtx[,1]
  v.mtx[,!grepl( "Delta_beta", names(V.hat), fixed = TRUE)] = v.mtx[,!grepl( "Delta_beta", names(V.hat), fixed = TRUE)]/aux.tot[1]
  
  vWv_inv = solve(t(samp[,wt0]*v.mtx)%*%v.mtx)
  f = c(t(1+(V-V.hat)%*%vWv_inv%*%t(v.mtx)))
  #f[f<0]=1e-10
  if(f_w){
    f_w1  = (V-V.hat)%*%vWv_inv
    vWv_w = lapply(1:n, function(i) outer(v.mtx[i,], v.mtx[i,]))
    f_w2  = vWv_inv%*%t(v.mtx)
    f_w = -sapply(1:n, function(i) f_w1%*%vWv_w[[i]]%*%f_w2)-
      (v.mtx)%*%vWv_inv%*%t(v.mtx)
    return(list(f = c(f), f_w = f_w))
  }else{return(list(f=c(f)))}
  
}
################################################################################################
beta.est.var = function(ds, fit, fit.i, sub, sub.wt, structure="combine"){
  #ds = ds.all; fit = as.formula(fit.y); fit.i = as.formula(gsub("x4", "x4.i", fit.y))
  #cyc = "cycle"; sub.cyl=1; method="linear"
  # design variables
  samp = ds$variables
  samp$psu = as.numeric(as.character(unlist(ds$cluster)))
  samp$strata = as.numeric(unlist(ds$strata))
  samp$wt = as.numeric(unlist(1/ds$allprob))
  samp0 = samp[samp[,sub]==1,]
  samp0$wt = samp0[,sub.wt]
  rm(ds)
  # outcome variable
  fit = as.formula(fit); fit.i = as.formula(fit.i)
  y = all.vars(fit)[1]
  # sampling design for the combined sample
  ds.all = svydesign(ids=~psu, strata = ~strata, weights=samp$wt, data=samp)
  # outcome model using surrogate variables
  imp.mdl = svyglm(fit.i, ds.all, family="binomial")
  p.tilde = imp.mdl$fitted.values
  x.tilde = model.matrix(imp.mdl)
  beta.tilde = imp.mdl$coefficients
  # influence functions
  U.tilde_beta.tilde=-t(c(samp$wt*p.tilde*(1-p.tilde))*x.tilde)%*%x.tilde
  v.mtx = t(solve(U.tilde_beta.tilde)%*%
          t(c(samp[,y]-p.tilde)*x.tilde))
  colnames(v.mtx)[1] = "x0" 
  colnames(v.mtx) = paste0("delta.", gsub(":", ".", colnames(v.mtx))) 

  v.mtx0 = v.mtx[samp[,sub]==1,]
  # subsample estimate of influence total 
  #Vs0.hat = c(samp0$wt)%*%v.mtx0
  # combined sample estimate of influence total
  VS.hat  = c(samp$wt)%*%v.mtx
  # GREG adjustment factor
    greg.out = greg_f(wt0=samp0$wt, v.mtx0=v.mtx0, VS.hat=VS.hat#, Vs0.hat=Vs0.hat
                      )
    f = greg.out$f
    f[f<0]=0
    eta = greg.out$eta
  # GREG weights for the subsample
  samp0$wt.c = samp0$wt*f
  # sampling design for the subsample, with GREG weights
  ds.c   = svydesign(ids=~psu, strata = ~strata, weights=~wt.c, data=samp0)
  # final outcome model
  mdl = svyglm(fit, ds.c, family="binomial")
  p = mdl$fitted.values
  X = model.matrix(mdl)
  u=c(samp0[,y]-p)*X
  beta.est = mdl$coefficient
  B = solve(t(samp0$wt*v.mtx0)%*%v.mtx0)%*%(t(samp0$wt*v.mtx0)%*%u)
  ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c*p*(1-p))*X)%*%X
  if(structure=="combine"){
    #partial derivative of the estimating equation system, w.r.t. the sample weights
    U_w = matrix(0,nrow(samp), ncol(X))
    U_w[samp[,sub]==1,] = samp0$wt.c*c(samp0[,y]-p)*X
    S1_w = matrix(0,nrow(v.mtx), ncol(v.mtx))
    S1_w[samp[,sub]==1,] = samp0$wt.c*v.mtx0
    S_w  = S1_w-samp$wt*v.mtx
    U.tilde_w = samp$wt*c(samp[,y]-p.tilde)*x.tilde
    Phi_w = cbind(U_w, S_w, U.tilde_w)
    #partial derivative of the estimating equation system, w.r.t. the parameters
    #U_beta.fun = function(x){
    #  odds.p = exp(X%*%x)
    #  p1 = odds.p/(1+odds.p)
    #  c(samp0$wt.c*c(samp0[,y]-p1))%*%X
    #}
    #U_beta.fun(x=beta.est)
    #U_beta = jacobian(func=U_beta.fun, x = beta.est)
    #U_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  t(c(f1*samp0$wt*c(samp0[,y]-p)))%*%X
    #}
    #U_eta.fun(x=eta)
    #U_eta = jacobian(func=U_eta.fun, x = eta)
    U_eta = t(t(c(samp0$wt*c(samp0[,y]-p))*v.mtx0)%*%X)
    
    U_beta.tilde.fun = function(x){
      odds.p.tilde = exp(x.tilde%*%x)
      p.tilde1 = odds.p.tilde/(1+odds.p.tilde)
      U.tilde_beta = -t(c(samp$wt*p.tilde1*(1-p.tilde1))*x.tilde)%*%x.tilde
      v.mtx1 = t(solve(U.tilde_beta)%*%
                   t(c(samp[,y]-p.tilde1)*x.tilde))
      v.mtx10 = v.mtx1[samp[,sub]==1,]
      f1 = c(t(1+eta%*%t(v.mtx10)))
      t(c(f1*samp0$wt*c(samp0[,y]-p)))%*%X
    }
    U_beta.tilde.fun(x=beta.tilde)
    U_beta.tilde = jacobian(func=U_beta.tilde.fun, x = beta.tilde)
    
    S_beta = matrix(0,length(beta.est),length(beta.est))
    #S_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  c(f1*samp0$wt)%*%v.mtx0-VS.hat
    #}
    #S_eta.fun(x=eta)
    #S_eta = jacobian(func=S_eta.fun, x = eta)
    S_eta = t(c(samp0$wt)*v.mtx0)%*%v.mtx0
    S_beta.tilde.fun = function(x){
      odds.p.tilde = exp(x.tilde%*%x)
      p.tilde1 = odds.p.tilde/(1+odds.p.tilde)
      U.tilde_beta = -t(c(samp$wt*p.tilde1*(1-p.tilde1))*x.tilde)%*%x.tilde
      v.mtx1 = t(solve(U.tilde_beta)%*%
                   t(c(samp[,y]-p.tilde1)*x.tilde))
      v.mtx10 = v.mtx1[samp[,sub]==1,]
      f1 = c(t(1+eta%*%t(v.mtx10)))
      c(f1*samp0$wt)%*%v.mtx10-c(samp$wt)%*%v.mtx1
    }
    #S_beta.tilde.fun(x=beta.tilde)
    S_beta.tilde = jacobian(func=S_beta.tilde.fun, x = beta.tilde)
    
    U.tilde_beta = matrix(0,length(beta.tilde),length(beta.est))
    U.tilde_eta  = matrix(0,length(beta.tilde),length(beta.tilde))
    #U.tilde_beta.tilde.fun = function(x){
    #  odds.p.tilde = exp(x.tilde%*%x)
    #  p.tlde1 = odds.p.tilde/(1+odds.p.tilde)
    #  c(samp$wt*c(samp[,y]-p.tlde1))%*%x.tilde
    #}
    #U.tilde_beta.tilde.fun(x=beta.tilde)
    #U.tilde_beta.tilde = jacobian(func=U.tilde_beta.tilde.fun, x = beta.tilde)
    #-t(c(samp$wt*p.tilde*(1-p.tilde))*x.tilde)%*%x.tilde
    
    Phi_theta = rbind(cbind(U_beta,       U_eta,       U_beta.tilde),
                      cbind(S_beta,       S_eta,       S_beta.tilde),
                      cbind(U.tilde_beta, U.tilde_eta, U.tilde_beta.tilde))
    US.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%U_eta%*%solve(S_eta)),
                   cbind(S_beta, solve(S_eta)))
    Phi_theta.inv = rbind(cbind(US.inv, -US.inv%*%rbind(U_beta.tilde, S_beta.tilde)%*%solve(U.tilde_beta.tilde)),
                          cbind(U.tilde_beta, U.tilde_eta, solve(U.tilde_beta.tilde))
    )
    #Phi_theta.inv1 = solve(Phi_theta)
    #t(solve(U.tilde_beta.tilde)%*%
    #    t(c(samp[,y]-p.tilde)*x.tilde))
    beta_w = t(Phi_theta.inv%*%t(Phi_w))
    TD.dat = as.data.frame(cbind(1, samp$psu, samp$strata, beta_w))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.all = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD)))
  }
  if(structure=="subsample"){
    samp0$wt0 = samp0$wt/samp$wt[samp[,sub]==1]
    U_w.1 = matrix(0,nrow(samp), ncol(X))
    U_w.1[samp[,sub]==1,] = samp0$wt0*(u-v.mtx0%*%B)
    U_w.1 = v.mtx%*%B+U_w.1
    TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, U_w.1))
    #U_w.1 = c(samp[,y]-predict(mdl, samp, type="response"))*model.matrix(lm(as.formula(gsub(y, "wt", fit.y)),data=samp))
    #colSums(samp0$wt.c*c(samp0[,y]-p)*X)
    #colSums(samp$wt*c(samp[,y]-predict(mdl, samp, type="response"))*as.matrix(cbind(1,samp$x1,samp$x3,samp$x4,samp$x3*samp$x4)))
    #colSums((samp$wt*v.mtx)%*%B)+colSums(samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B))
    #e = matrix(0,nrow(samp), ncol(X))
    #e[samp[,sub]==1,] = samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B)
    #var((samp$wt*v.mtx)%*%B+e)
    #S_w.1 = matrix(0,nrow(v.mtx), ncol(v.mtx))
    #S_w.1 = c(v.mtx%*%eta)*v.mtx
    #Phi_w.1 = cbind(U_w.1, S_w.1)
    #TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, Phi_w.1))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.Phi.1 = vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD))
    #var.Phi.2 = t(c((samp0$wt0-1)*samp$wt[samp[,sub]==1]*samp0$wt)*cbind(u,v.mtx0))%*%cbind(u,v.mtx0)
    #var.all = diag(Phi_theta.inv%*%(var.Phi.1+var.Phi.2)%*%t(Phi_theta.inv))
    var.all = diag(solve(U_beta)%*%(var.Phi.1)%*%t(solve(U_beta)))
  }
  return(list(beta.est = beta.est, beta.var = vcov(mdl),
              beta.tilde = beta.tilde, beta.tilde.var = vcov(imp.mdl),
              eta = eta,
              var.all = var.all))
}

################################################################################################
# ds: survey design of the full sample
# fit: outcome model
# v.mtx:matrix of auxiliary variables in the full sample
# sub: variable name of the indicator for subsample
# sub.wt: weight variable name for the sub sample
# structure: data structure "combine" for combining multiple cycles; "subsample" for a random sub sample 
beta.est.varPop = function(ds, fit, v.mtx0, X.pop){
  #ds = ds.all; fit = fit.y; v.mtx =v.mtx.i;  
  #sub = "cycle"; method="linear", sub.wt = "wt1", structure="comb"
  # design variables
  samp0 = ds$variables
  samp0$psu = as.numeric(as.character(unlist(ds$cluster)))
  samp0$strata = as.numeric(unlist(ds$strata))
  samp0$wt = as.numeric(unlist(1/ds$allprob))
  rm(ds)
  # outcome variable
  fit = as.formula(fit); 
  y = all.vars(fit)[1]
  # sample estimate of totals
  X.hat  = c(c(samp0$wt)%*%v.mtx0)
  # GREG adjustment factor
  greg.out = greg_f(wt0=samp0$wt, v.mtx0=v.mtx0, VS.hat=X.pop#, Vs0.hat=X.hat
                    )
  f = greg.out$f
  f[f<0]=0
  eta = greg.out$eta

  # GREG weights for the subsample
  samp0$wt.c = samp0$wt*f
  # sampling design for the subsample, with GREG weights
  ds.c   = svydesign(ids=~psu, strata = ~strata, weights=~wt.c, data=samp0)
  # final outcome model
  mdl = svyglm(fit, ds.c, family="binomial")
  p = mdl$fitted.values
  X = model.matrix(mdl)
  u=c(samp0[,y]-p)*X
  beta.est = mdl$coefficient
  B = solve(t(samp0$wt*v.mtx0)%*%v.mtx0)%*%(t(samp0$wt*v.mtx0)%*%u)
  ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c*p*(1-p))*X)%*%X
  U_eta = t(samp0$wt*u)%*%v.mtx0
  S_beta = matrix(0,ncol(v.mtx0),length(beta.est))
  S_eta = t(c(samp0$wt)*v.mtx0)%*%v.mtx0
  Phi_theta = rbind(cbind(U_beta, U_eta),
                    cbind(S_beta, S_eta))
  Phi_theta.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%t(B)),
                        cbind(S_beta, solve(S_eta)))
  U_w = samp0$wt.c*u
  S_w = samp0$wt.c*v.mtx0
  Phi_w = cbind(U_w, S_w)
  beta_w = -t(Phi_theta.inv%*%t(Phi_w))
  TD.dat = as.data.frame(cbind(1, samp0$psu, samp0$strata, beta_w))
  names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
  ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
  var.all = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD)))

  return(list(beta.est = beta.est, beta.var = vcov(mdl),
              eta = eta,
              var.all = var.all))
}


beta.est.var0 = function(ds, fit, v.mtx, sub, sub.wt, structure){
  #ds = ds.all; fit = fit.y; v.mtx =v.mtx.i;  
  #sub = "cycle"; method="linear", sub.wt = "wt1", structure="comb"
  # design variables
  samp = ds$variables
  samp$psu = as.numeric(as.character(unlist(ds$cluster)))
  samp$strata = as.numeric(unlist(ds$strata))
  samp$wt = as.numeric(unlist(1/ds$allprob))
  samp0 = samp[samp[,sub]==1,]
  samp0$wt = samp0[,sub.wt]
  rm(ds)
  # outcome variable
  fit = as.formula(fit); 
  y = all.vars(fit)[1]
  v.mtx0 = v.mtx[samp[,sub]==1,]
  # subsample estimate of influence total 
  #Vs0.hat = c(c(samp0$wt)%*%v.mtx0)
  # combined sample estimate of influence total
  VS.hat  = c(c(samp$wt)%*%v.mtx)
  # GREG adjustment factor
  #if(method=="linear"){
  greg.out = greg_f(wt0=samp0$wt, v.mtx0=v.mtx0, VS.hat=VS.hat#, Vs0.hat=Vs0.hat
                    )
  f = greg.out$f
  f[f<0]=0
  eta = greg.out$eta
  #c_eta = v.mtx0
  #c_v = eta
  #}
  # GREG weights for the subsample
  samp0$wt.c = samp0$wt*f
  # sampling design for the subsample, with GREG weights
  ds.c   = svydesign(ids=~psu, strata = ~strata, weights=~wt.c, data=samp0)
  # final outcome model
  mdl = svyglm(fit, ds.c, family="binomial")
  p = mdl$fitted.values
  X = model.matrix(mdl)
  u=c(samp0[,y]-p)*X
  beta.est = mdl$coefficient
  B = solve(t(samp0$wt*v.mtx0)%*%v.mtx0)%*%(t(samp0$wt*v.mtx0)%*%u)
  ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c*p*(1-p))*X)%*%X
  if(structure=="combine"){
    #partial derivative of the estimating equation system, w.r.t. the parameters
    #U_beta.fun = function(x){
    #  odds.p = exp(X%*%x)
    #  p1 = odds.p/(1+odds.p)
    #  c(samp0$wt.c*c(samp0[,y]-p1))%*%X
    #}
    #U_beta.fun(x=beta.est)
    #U_beta = jacobian(func=U_beta.fun, x = beta.est)
    #U_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  t(c(f1*samp0$wt*c(samp0[,y]-p)))%*%X
    #}
    #U_eta.fun(x=eta)
    #U_eta = jacobian(func=U_eta.fun, x = eta)
    U_eta = t(samp0$wt*u)%*%v.mtx0
    S_beta = matrix(0,length(beta.est),length(beta.est))
    #S_eta.fun = function(x){
    #  f1 = c(t(1+x%*%t(v.mtx0)))
    #  c(f1*samp0$wt)%*%v.mtx0-VS.hat
    #}
    #S_eta.fun(x=eta)
    #S_eta = jacobian(func=S_eta.fun, x = eta)
    S_eta = t(c(samp0$wt)*v.mtx0)%*%v.mtx0
    Phi_theta = rbind(cbind(U_beta, U_eta),
                      cbind(S_beta, S_eta))
    #t(solve(U.tilde_beta.tilde)%*%
    #    t(c(samp[,y]-p.tilde)*x.tilde))
    Phi_theta.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%t(B)),
                          cbind(S_beta, solve(S_eta)))
    #Phi_theta.inv = solve(Phi_theta)
    #round(Phi_theta.inv,5)==round(solve(Phi_theta),5)
    #partial derivative of the estimating equation system, w.r.t. the sample weights
    U_w = matrix(0,nrow(samp), ncol(X))
    U_w[samp[,sub]==1,] = samp0$wt.c*u
    S1_w = matrix(0,nrow(v.mtx), ncol(v.mtx))
    S1_w[samp[,sub]==1,] = samp0$wt.c*v.mtx0
    S_w  = S1_w-samp$wt*v.mtx
    Phi_w = cbind(U_w, S_w)
    beta_w = -t(Phi_theta.inv%*%t(Phi_w))
    #eta_w = -t(solve(S_eta)%*%t(S_w))
    #beta_w.c = -t(solve(U_beta)%*%t(u))
    #U.dat = as.data.frame(cbind(samp0$wt.c, samp0$psu, samp0$strata, beta_w.c))
    #names(U.dat)=c("wt.c", paste0("V", 1:(ncol(U.dat)-1)))
    #ds.U = svydesign(ids=~psu, strata =~strat, weights=~wt.c, data=U.dat)
    #var.beta = diag(vcov(svytotal(as.formula(paste0("~", paste(names(U.dat)[-1], collapse="+"))), ds.U)))
    #diag(vcov(mdl))
    #S.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, eta_w))
    #names(S.dat)=c("wt", "cycle", paste0("V", 1:(ncol(S.dat)-2)))
    #ds.S = svydesign(ids=~psu, strata =~strata, weights=~wt, data=S.dat)
    #var.eta = diag(vcov(svytotal(as.formula(paste0("~", paste(names(S.dat)[-c(1,2)], collapse="+"))), ds.S)))
    
    #Phi.dat = as.data.frame(cbind(samp$wt, samp[,cyc], Phi_w))
    #names(Phi.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(Phi.dat)-3)))
    #ds.Phi = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=Phi.dat)
    #var.Phi = vcov(svytotal(as.formula(paste0("~", paste(names(Phi.dat)[-c(1,2)], collapse="+"))), ds.Phi))
    #diag((Phi_theta.inv%*%var.Phi)%*%t(Phi_theta.inv))
    TD.dat = as.data.frame(cbind(1, samp$psu, samp$strata, beta_w))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.all = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD)))
  }
  if(structure=="subsample"){
    samp0$wt0 = samp0$wt/samp$wt[samp[,sub]==1]
    U_w.1 = matrix(0,nrow(samp), ncol(X))
    U_w.1[samp[,sub]==1,] = samp0$wt0*(u-v.mtx0%*%B)
    U_w.1 = v.mtx%*%B+U_w.1
    TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, U_w.1))
    #U_w.1 = c(samp[,y]-predict(mdl, samp, type="response"))*model.matrix(lm(as.formula(gsub(y, "wt", fit.y)),data=samp))
    #colSums(samp0$wt.c*c(samp0[,y]-p)*X)
    #colSums(samp$wt*c(samp[,y]-predict(mdl, samp, type="response"))*as.matrix(cbind(1,samp$x1,samp$x3,samp$x4,samp$x3*samp$x4)))
    #colSums((samp$wt*v.mtx)%*%B)+colSums(samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B))
    #e = matrix(0,nrow(samp), ncol(X))
    #e[samp[,sub]==1,] = samp0$wt*(c(samp0[,y]-p)*X-v.mtx0%*%B)
    #var((samp$wt*v.mtx)%*%B+e)
    #S_w.1 = matrix(0,nrow(v.mtx), ncol(v.mtx))
    #S_w.1 = c(v.mtx%*%eta)*v.mtx
    #Phi_w.1 = cbind(U_w.1, S_w.1)
    #TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, Phi_w.1))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.Phi.1 = vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD))
    #var.Phi.2 = t(c((samp0$wt0-1)*samp$wt[samp[,sub]==1]*samp0$wt)*cbind(u,v.mtx0))%*%cbind(u,v.mtx0)
    #var.all = diag(Phi_theta.inv%*%(var.Phi.1+var.Phi.2)%*%t(Phi_theta.inv))
    var.all = diag(solve(U_beta)%*%(var.Phi.1)%*%t(solve(U_beta)))
  }
  #c(var.beta, var.eta)
  #var.all[1:5]
  return(list(beta.est = beta.est, beta.var = vcov(mdl),
              eta = eta,
              var.all = var.all))
}

greg_f = function(wt0, v.mtx0, VS.hat#, Vs0.hat=NULL
                  ){
  #wt0 = samp0$wt
  vWv_inv = solve(t(wt0*v.mtx0)%*%v.mtx0)
  #if(is.null(Vs0.hat)) 
  Vs0.hat=c(c(wt0)%*%v.mtx0)
  eta = c(c(VS.hat-Vs0.hat)%*%vWv_inv)
  f = c(t(1+eta%*%t(v.mtx0)))
  return(list(eta = eta, f = f))    
}



raking_f = function(wt0, v.mtx0, VS.hat, N.S, Vs0.hat, eta0){
  wt0 = samp0$wt
  N.S=sum(samp$wt)
  dat = as.data.frame(cbind(wt=wt0, v.mtx0))
  names(dat)[-1] = paste0("V", 1:(ncol(dat)-1))
  ds = svydesign(ids=~1, data=dat, weights=~wt)
  names(VS.hat)=names(dat)[-1]
  f=weights(calibrate(ds, as.formula(paste0("~", paste(names(dat)[-1], collapse = "+"))), 
            c("(Intercept)"=N.S, VS.hat), calfun="linear"))/wt0
  solve(t(v.mtx0)%*%v.mtx0, t(v.mtx0)%*%f)
  
}

#weights(calibrate(ds.wt, as.formula(paste0("~", paste0(names(aux.tot)[-c(1,length(aux.tot))],collapse="+"))), 
#                  aux.tot[-c(length(aux.tot))], calfun="raking"))


winsorize <- function(x, quant = 0.99) {
  q <- quantile(x, probs = c(1 - quant, quant), na.rm = TRUE)
  x[x > q[2]] <- q[2]
  x[x < q[1]] <- q[1]
  return(x)
}



#Function to calculate the V(S)
# ds: if the individual level is accessible then ds=ds.S, if not we use the ds1 to estimate
# formula: fit.y[i, fm.idx]
compute_covariance_S <- function(ds, formula, structure = "combined",var_beta_star = NULL, data = NULL, beta = NULL,ds2=NULL, n1=NULL,n2=NULL,matched_idx=NULL) {
  library(survey)  # Ensure survey package is loaded
  
  weights_S <- weights(ds)
  
  if (structure == "combined") {
    if(is.null(ds2)||is.null(n1)||is.null(n2)||is.null(matched_idx)){
      stop("Error: ds2,n1,n2 must be provided for structure = 'combined'.") 
    }
    # Fit the regression model
    fit <- svyglm(as.formula(formula), ds, family = "binomial")
    
    # Extract fitted values and model matrix
    p_hat_S <- fit$fitted.values
    X_S <- model.matrix(fit)
    response_var <- all.vars(formula)[1]
    y <- ds$variables[[response_var]]
    u_S <- (y[matched_idx] - p_hat_S[matched_idx]) * X_S[matched_idx, ]
    
    #u_S_win <- apply(u_S, 2, winsorize)
    
    # Create column names dynamically
    colnames(u_S) <- paste0("u_S_", seq_len(ncol(u_S)))
    
    # Update the survey design with u_S as additional variables (avoiding for-loop)
    ds2$variables <- cbind(ds2$variables, as.data.frame(u_S))
    
    # Construct formula dynamically for svytotal()
    u_formula <- as.formula(paste("~", paste(colnames(u_S), collapse = " + ")))
    
    
    
    N2 <- sum(weights(ds2))  # approximate population size based on weights
    fpc_factor <- 1 - (n2 / N2)
    
    # Compute variance using svytotal()
    var_S <- (n2^2 / (n1 + n2)^2)*vcov(svytotal(u_formula, ds2))*fpc_factor
    
  } else if (structure == "model-based") {
    # Model-based structure: require var_beta_star
    if (is.null(var_beta_star)) {
      stop("Error: var_beta_star must be provided for structure = 'model-based1'.")
    }
    
    # Compute model-based p_hat
    X_S <- model.matrix(lm(as.formula(formula), data = data))
    odds.hat <- exp(X_S %*% beta)
    p_hat_S <- odds.hat / (1 + odds.hat)
    w_S <- as.numeric(weights_S * p_hat_S * (1 - p_hat_S))  # Vector of weights
    
    U_beta_S <- -crossprod(X_S, w_S * X_S)
    
    # Compute variance using var_beta_star
    var_S <- U_beta_S %*% var_beta_star %*% t(U_beta_S)
  }  else if (structure == "model-based2") {
    # Model-based structure: require var_beta_star
    if (is.null(var_beta_star)) {
      stop("Error: var_beta_star must be provided for structure = 'model-based1'.")
    }
    
    # Compute variance using var_beta_star
    var_S <- solve(var_beta_star)
  }else {
    stop("Invalid structure specified. Choose 'combined' or 'model-based'.")
  }
  
  return(var_S)
}

#Function to calculate the V(s) and V11, V12,V22
# ds:ds1.cwt.i the design after the calibrate weight
# u_star: ui
# formula_i:fit.y[i, 1]
compute_covariance_s <- function(ds, formula_i, structure,wt0,u_star,ds1=NULL,n1=NULL,n2=NULL) {
  library(survey)  # Load survey package
  if (structure == "model-based"||structure == "model-based1") {
  # Extract survey weights
  weights <- weights(ds)
  
  # Fit the survey-weighted regression model
  fit_i <- svyglm(as.formula(formula_i), ds, family = "binomial")
  
  # Extract fitted values and model matrices
  p_hat <- fit_i$fitted.values
  X <- model.matrix(fit_i)
  
  
  response_var <- all.vars(formula_i)[1]
  y <- ds$variables[[response_var]]
  # Compute influence functions
  u <- (y - p_hat) * X  # Influence function for model
  
  # Convert influence functions to a data frame
  u_df <- as.data.frame(u)
  colnames(u_df) <- paste0("u_", seq_len(ncol(u_df)))  # Name columns dynamically
  
  u_star_df <- as.data.frame(u_star)
  colnames(u_star_df) <- paste0("u_star_", seq_len(ncol(u_star_df)))  # Name columns dynamically
  
  # Combine u and u_star into one dataset
  combined_df <- cbind(u_df, u_star_df)
  
  # Update survey design with combined influence function variables
  ds$variables <- cbind(ds$variables, combined_df)
  
  # Construct dynamic formula for svytotal()
  combined_formula <- as.formula(paste("~", paste(colnames(combined_df), collapse = " + ")))
  
  # Compute full variance-covariance matrix using svytotal()
  combined_cov_matrix <- vcov(svytotal(combined_formula, ds))
  
  # Extract matrix dimensions
  d1 <- ncol(u_df)  # Number of columns in u
  d2 <- ncol(u_star_df)  # Number of columns in u_star
  
  # Extract submatrices
  V11 <- combined_cov_matrix[1:d1, 1:d1]  # Variance of u
  V22 <- combined_cov_matrix[(d1 + 1):(d1 + d2), (d1 + 1):(d1 + d2)]  # Variance of u_star
  V12 <- combined_cov_matrix[1:d1, (d1 + 1):(d1 + d2)]  # Covariance between u and u_star
  } else if (structure == "combined") {
    if(is.null(ds1)||is.null(n1)||is.null(n2)){
      stop("Error: ds1,n1,n2 must be provided for structure = 'combined'.") 
    }
    # Extract survey weights
    weights <- weights(ds)
    
    # Fit the survey-weighted regression model
    fit_i <- svyglm(as.formula(formula_i), ds, family = "binomial")
    
    # Extract fitted values and model matrices
    p_hat <- fit_i$fitted.values
    X <- model.matrix(fit_i)
    
    response_var <- all.vars(formula_i)[1]
    y <- ds$variables[[response_var]]
    
    # Compute influence functions
    u <- (y - p_hat) * X  # Influence function for model
    
    F_i <-weights/wt0
    
    r <-n1/(n1+n2)
    
    u_df <- as.data.frame(F_i*u)
    
    colnames(u_df) <- paste0("u_", seq_len(ncol(u_df)))  # Name columns dynamically
    
    u_star_df <- as.data.frame((F_i-r)*u_star)
    colnames(u_star_df) <- paste0("u_star_", seq_len(ncol(u_star_df)))  # Name columns dynamically
    
    
    # Combine u and u_star into one dataset
    combined_df <- cbind(u_df, u_star_df)
    
    # Update survey design with combined influence function variables
    ds1$variables <- cbind(ds1$variables, combined_df)
    
    # Construct dynamic formula for svytotal()
    combined_formula <- as.formula(paste("~", paste(colnames(combined_df), collapse = " + ")))
    
    # Compute full variance-covariance matrix using svytotal()
    combined_cov_matrix <- vcov(svytotal(combined_formula, ds1))
    
    # Extract matrix dimensions
    d1 <- ncol(u_df)  # Number of columns in u
    d2 <- ncol(u_star_df)  # Number of columns in u_star
    
    # Extract submatrices
    V11 <- combined_cov_matrix[1:d1, 1:d1]  # Variance of u
    V22 <- combined_cov_matrix[(d1 + 1):(d1 + d2), (d1 + 1):(d1 + d2)]  # Variance of u_star
    V12 <- combined_cov_matrix[1:d1, (d1 + 1):(d1 + d2)]  # Covariance between u and u_star
    
  }
  
  return(list(combined_cov_matrix = combined_cov_matrix, V11 = V11, V12 = V12, V22 = V22))
}





# Custom variance estimation function using the first derivative (influence function)
custom_variance <- function(design, formula) {
  # Fit the survey-weighted logistic regression model
  fit <- svyglm(as.formula(formula), design, family = "binomial")
  
  # Extract design matrix, weights, and fitted probabilities
  X <- model.matrix(fit)  # Predictor variables (design matrix)
  weights <- weights(design)  # Extract survey weights
  p_hat <- fit$fitted.values  # Predicted probabilities
  
  # Compute the influence function (first derivative)
  response_var <- all.vars(formula)[1]
  y <- design$variables[[response_var]]
  u <- (y - p_hat) * X  # Influence function (score function)
  
  var_UB <-t(weights * u)%*%(weights * u)  # Outer product of influence functions
  
  
  # Final variance estimate: directly use the meat matrix
  
  W <- as.numeric(weights * p_hat * (1 - p_hat))  # Vector of weights
  
  U_beta <- -t(X)%*% (W * X)
  
  
  var_beta_1 <- solve(U_beta)%*%var_UB%*%solve(U_beta)
  
  return(list(variance_matrix = var_beta_1,U_beta=U_beta,u=u))
}





# Combined function to calculate the analytical variance
calibrate_beta_variance <- function(ds_s,ds_S,formula,formula_i,u_star,wt0,structure = "combined",var_beta_star = NULL,data=NULL,beta=NULL,ds2=NULL,n1=NULL,n2=NULL,ds1=NULL,matched_idx=NULL){
  
  # Using variance function to calculate U_beta
  result <- custom_variance(ds_s,formula_i)
  u <-  result$u
  U_beta <- result$U_beta
  
  # Calculate Q_eta
  Q_eta <- -t(u_star) %*% (as.numeric(wt0) * u_star)
  
  U_eta <- -t(u_star) %*% (as.numeric(wt0) * u)
  
  
  # Q_eta inverse
  Q_eta_inv <- solve(Q_eta)
  
  
  U_beta_inv <- solve(U_beta)
  
  # Calculate b
  b <- -U_beta_inv %*% t(U_eta) %*% solve(Q_eta)
 
  
   
  # Calculate V(S)
  var_S1 <- compute_covariance_S(ds_S,formula,structure,var_beta_star,data,beta,ds2,n1,n2,matched_idx)
  
  # Calculate V(s)
  stacked_cov_matrix <- compute_covariance_s(ds_s, formula_i, structure,wt0,u_star,ds1,n1,n2)
  V11 <- stacked_cov_matrix$V11
  V12 <- stacked_cov_matrix$V12
  V22 <- stacked_cov_matrix$V22
  
  var_beta_calib <- U_beta_inv%*%V11%*%t(U_beta_inv) + b%*%t(V12)%*%t(U_beta_inv) + U_beta_inv%*%V12%*%t(b) + b%*%V22%*%t(b)+b%*%var_S1%*%t(b)
  
  
  return(list(var_beta_calib=var_beta_calib, var_S1=var_S1,V11=V11,V12=V12,V22=V22,U_beta=U_beta,Q_eta=Q_eta,Q_eta_inv=Q_eta_inv,U_beta_inv=U_beta_inv,b=b,U_eta=U_eta))
}



calib_if_logit <- function(wt0, X, y, p_hat,u=NULL,VS.hat = NULL, ridge = 1e-8) {
  # wt0  : n-vector internal base weights
  # X    : n x p model matrix (same as used in score method)
  # y    : n-vector binary response
  # p_hat: n-vector probabilities computed from external/surrogate beta
  # VS.hat: p-vector target total, default 0
  # ridge: small stabilizer for matrix inversion
  
  # Score contributions u_i = (y - p_hat) x_i
  if(is.null(u)){
  u <- c(y - p_hat) * X  # n x p
  }
  # Plug-in weighted Jacobian/information using internal weights + p_hat
  W <- as.numeric(wt0 * p_hat * (1 - p_hat))  # length n
  U_beta_hat <- -t(X) %*% (W * X)            # p x p
  
  # Stabilize in case of near singularity
  # U_beta_hat <- U_beta_hat - diag(ridge, nrow(U_beta_hat))
  
  # Estimated influence function contributions: phi_i = U^{-1} u_i
  # (n x p)
  v.mtx0 <- t(solve(U_beta_hat, t(u)))
  
  if (is.null(VS.hat)) VS.hat <- rep(0, ncol(X))
  

  return(list(U_beta_hat = U_beta_hat,v.mtx0=v.mtx0))
}
