
beta.est.var = function(ds, fit, fit.i, sub, sub.wt, structure = "combine"){
  # --- design variables
  samp = ds$variables
  samp$psu    = as.numeric(as.character(unlist(ds$cluster)))
  samp$strata = as.numeric(unlist(ds$strata))
  samp$wt     = as.numeric(unlist(1/ds$allprob))
  samp0       = samp[samp[, sub] == 1, ]
  samp0$wt    = samp0[, sub.wt]
  rm(ds)
  
  # --- outcomes & formulas
  fit    = as.formula(fit)
  fit.i  = as.formula(fit.i)
  y      = all.vars(fit)[1]
  
  # --- combined-sample design
  ds.all = svydesign(ids = ~psu, strata = ~strata, weights = ~wt, data = samp)
  
  # --- outcome model using surrogate variables
  imp.mdl     = svyglm(fit.i, ds.all, family = "binomial")
  p.tilde     = imp.mdl$fitted.values
  x.tilde     = model.matrix(imp.mdl)
  beta.tilde  = imp.mdl$coefficients
  
  # --- influence functions (tilde-model)
  U.tilde_beta.tilde = -t(c(samp$wt * p.tilde * (1 - p.tilde)) * x.tilde) %*% x.tilde
  v.mtx = t(solve(U.tilde_beta.tilde) %*% t(c(samp[, y] - p.tilde) * x.tilde))
  
  # --- make names safe for formulas (one-pass, robust)
  colnames(v.mtx)[1] = "x0"
  colnames(v.mtx)    = paste0("delta.", gsub(":", ".", colnames(v.mtx)))
  colnames(v.mtx)    = make.names(colnames(v.mtx), unique = TRUE)
  
  v.mtx0 = v.mtx[samp[, sub] == 1, , drop = FALSE]
  
  # --- combined sample estimate of influence total (NAME IT!)
  VS.hat  = as.numeric(c(samp$wt) %*% v.mtx)
  names(VS.hat) = colnames(v.mtx)
  
  # --- calibrate on subsample with SAFE formula
  ds.calib0 <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt,
                         data = cbind(samp0, v.mtx0))
  form.cal  <- reformulate(termlabels = colnames(v.mtx0), intercept = FALSE)
  
  calib <- survey::calibrate(
    design     = ds.calib0,
    formula    = form.cal,
    population = VS.hat,
    calfun     = "linear"
  )
  
  # optional: clip negatives if any
  wts <- weights(calib)
  wts[wts < 0] <- 0
  calib <- svydesign(ids = ~psu, strata = ~strata, weights = ~w,
                     data = transform(calib$variables, w = wts))
  
  # --- final outcome model on calibrated design
  mdl              = svyglm(fit, calib, family = "binomial")
  svyclib.beta.est = coef(mdl)
  svyclib.beta.var = diag(vcov(mdl))
  
  # --- GREG adjustment factor (your helper must exist)
  greg.out = greg_f(wt0 = samp0$wt, v.mtx0 = v.mtx0, VS.hat = VS.hat)
  f   = greg.out$f;  f[f < 0] = 0
  eta = greg.out$eta
  
  # --- GREG weights for subsample & design
  samp0$wt.c = samp0$wt * f
  ds.c       = svydesign(ids = ~psu, strata = ~strata, weights = ~wt.c, data = samp0)
  mdl        = svyglm(fit, ds.c, family = "binomial")
  
  # pieces used below
  p  = mdl$fitted.values
  X  = model.matrix(mdl)
  u  = c(samp0[, y] - p) * X
  svybeta.est = coef(mdl)
  
  samp$sub1 = (samp[, sub] == 1)
  samp      = cbind(samp, v.mtx)
  
  # regression of u on v.mtx0 (B)
  B = solve(t(samp0$wt * v.mtx0) %*% v.mtx0) %*% (t(samp0$wt * v.mtx0) %*% u)
  
  ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c * p * (1 - p)) * X) %*% X
  
  if (structure == "combine") {
    # partial derivative wrt weights
    U_w         = matrix(0, nrow(samp), ncol(X))
    U_w[samp[, sub] == 1, ] = samp0$wt.c * u
    
    S1_w        = matrix(0, nrow(v.mtx), ncol(v.mtx))
    S1_w[samp[, sub] == 1, ] = samp0$wt.c * v.mtx0
    S_w         = S1_w - samp$wt * v.mtx
    
    U.tilde_w   = samp$wt * c(samp[, y] - p.tilde) * x.tilde
    Phi_w       = cbind(U_w, S_w, U.tilde_w)
    
    U_eta = t(t(c(samp0$wt * c(samp0[, y] - p)) * v.mtx0) %*% X)
    
    U_beta.tilde.fun = function(x){
      odds.p.tilde = exp(x.tilde %*% x)
      p.tilde1     = odds.p.tilde / (1 + odds.p.tilde)
      U.tilde_beta = -t(c(samp$wt * p.tilde1 * (1 - p.tilde1)) * x.tilde) %*% x.tilde
      v.mtx1       = t(solve(U.tilde_beta) %*% t(c(samp[, y] - p.tilde1) * x.tilde))
      v.mtx10      = v.mtx1[samp[, sub] == 1, , drop = FALSE]
      f1           = c(t(1 + eta %*% t(v.mtx10)))
      t(c(f1 * samp0$wt * c(samp0[, y] - p))) %*% X
    }
    U_beta.tilde.fun(x = beta.tilde)
    U_beta.tilde = jacobian(func = U_beta.tilde.fun, x = beta.tilde)
    
    S_beta = matrix(0, length(svybeta.est), length(svybeta.est))
    S_eta  = t(c(samp0$wt) * v.mtx0) %*% v.mtx0
    
    S_beta.tilde.fun = function(x){
      odds.p.tilde = exp(x.tilde %*% x)
      p.tilde1     = odds.p.tilde / (1 + odds.p.tilde)
      U.tilde_beta = -t(c(samp$wt * p.tilde1 * (1 - p.tilde1)) * x.tilde) %*% x.tilde
      v.mtx1       = t(solve(U.tilde_beta) %*% t(c(samp[, y] - p.tilde1) * x.tilde))
      v.mtx10      = v.mtx1[samp[, sub] == 1, , drop = FALSE]
      f1           = c(t(1 + eta %*% t(v.mtx10)))
      c(f1 * samp0$wt) %*% v.mtx10 - c(samp$wt) %*% v.mtx1
    }
    S_beta.tilde = jacobian(func = S_beta.tilde.fun, x = beta.tilde)
    
    U.tilde_beta = matrix(0, length(beta.tilde), length(svybeta.est))
    U.tilde_eta  = matrix(0, length(beta.tilde), length(beta.tilde))
    
    US.inv = rbind(
      cbind(solve(U_beta), -solve(U_beta) %*% U_eta %*% solve(S_eta)),
      cbind(S_beta,        solve(S_eta))
    )
    Phi_theta.inv = rbind(
      cbind(US.inv, -US.inv %*% rbind(U_beta.tilde, S_beta.tilde) %*% solve(U.tilde_beta.tilde)),
      cbind(U.tilde_beta, U.tilde_eta, solve(U.tilde_beta.tilde))
    )
    
    beta_w  = t(Phi_theta.inv %*% t(Phi_w))
    US_w    = cbind(U_w, S_w)
    beta_w1 = -t(US.inv %*% t(US_w))
    
    TD.dat = as.data.frame(cbind(1, samp$psu, samp$strata, beta_w, beta_w1, US_w))
    names(TD.dat) = c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat) - 3)))
    ds.TD = svydesign(ids = ~psu, strata = ~strata, weights = ~wt, data = TD.dat)
    
    var.all  = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[1:ncol(beta_w)+3], collapse = "+"))), ds.TD)))
    var.all1 = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[1:ncol(beta_w1)+3+ncol(beta_w)], collapse = "+"))), ds.TD)))
    var.US_w = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[1:ncol(US_w)+3+ncol(beta_w)+ncol(beta_w1)], collapse = "+"))), ds.TD)))
  }
  
  if (structure == "subsample") {
    samp0$wt0 = samp0$wt / samp$wt[samp[, sub] == 1]
    U_w.1     = matrix(0, nrow(samp), ncol(X))
    U_w.1[samp[, sub] == 1, ] = samp0$wt0 * (u - v.mtx0 %*% B)
    U_w.1     = v.mtx %*% B + U_w.1
    
    TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, U_w.1))
    names(TD.dat) = c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat) - 3)))
    ds.TD   = svydesign(ids = ~psu, strata = ~strata, weights = ~wt, data = TD.dat)
    var.Phi.1 = vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse = "+"))), ds.TD))
    var.all   = diag(solve(U_beta) %*% (var.Phi.1) %*% t(solve(U_beta)))
    var.all1  = NA
  }
  
  return(list(
    beta.tilde       = imp.mdl$coefficients,
    beta.tilde.var   = vcov(imp.mdl),
    beta.est         = svybeta.est,
    beta.est.svyclib = svyclib.beta.est,
    beta.var.svyclib = svyclib.beta.var,
    var.all          = var.all,
    var.all1         = var.all1
  ))
}


######################################################################################################
beta.est.varA = function(ds, fit, fit.i, sub, sub.wt, structure="combine", beta.wt=F){
  samp = ds$variables
  samp$psu = as.numeric(as.character(unlist(ds$cluster)))
  samp$strata = as.numeric(unlist(ds$strata))
  samp$wt = as.numeric(unlist(1/ds$allprob))
  samp0 = samp[samp[,sub]==1,]
  samp0$wt = samp0[,sub.wt] # cycle 1 weight (final)
  samp0$wt0 = samp0$wt/samp$wt[samp[,sub]==1]
  
  rm(ds)
  # outcome variable
  fit = as.formula(fit); fit.i = as.formula(fit.i)
  y = all.vars(fit)[1]
  # outcome model using surrogate variables
  if(beta.wt){
    ds.all = svydesign(ids=~psu, strata = ~strata, weights=samp$wt, data=samp)
    # outcome model using surrogate variables
    imp.mdl = svyglm(fit.i, ds.all, family="binomial")
  }else{
    imp.mdl = glm(fit.i, family="binomial", data = samp)
  }
  p.tilde = imp.mdl$fitted.values
  x.tilde = model.matrix(imp.mdl)
  beta.tilde = imp.mdl$coefficients
  
  # influence functions
  U.tilde_w = c(c(samp[,y]-p.tilde))*x.tilde
  U.tilde_beta.tilde=-t(c(p.tilde*(1-p.tilde))*x.tilde)%*%x.tilde
  v.mtx = U.tilde_w
  colnames(v.mtx)[1] = "x0" 
  colnames(v.mtx) = paste0("delta.", gsub(":", ".", colnames(v.mtx))) 
  v.mtx0 = v.mtx[samp[,sub]==1,]
  # subsample estimate of influence total 
  #Vs0.hat = c(samp0$wt)%*%v.mtx0
  # combined sample estimate of influence total
  VS.hat  = colSums(v.mtx)
  # GREG adjustment factor
  greg.out = greg_f(wt0=samp0$wt0, v.mtx0=v.mtx0, VS.hat=VS.hat#, Vs0.hat=Vs0.hat
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
  B = solve(t(samp0$wt0*v.mtx0)%*%v.mtx0)%*%(t(samp0$wt0*v.mtx0)%*%u)
    ################### Variance calculation ###################
  U_beta = -t(c(samp0$wt.c*p*(1-p))*X)%*%X
  if(structure=="combine"){
    #partial derivative of the estimating equation system, w.r.t. the sample weights
    U_w = matrix(0,nrow(samp), ncol(X))
    U_w[samp[,sub]==1,] = samp0$wt.c*u
    S1_w = matrix(0,nrow(v.mtx), ncol(v.mtx))
    S1_w[samp[,sub]==1,] = samp0$wt0*f*v.mtx0
    S_w  = S1_w-v.mtx
    Phi_w = cbind(U_w, S_w, U.tilde_w)
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
    S_eta = t(c(samp0$wt0)*v.mtx0)%*%v.mtx0
    S_beta.tilde.fun = function(x){
      odds.p.tilde = exp(x.tilde%*%x)
      p.tilde1 = odds.p.tilde/(1+odds.p.tilde)
      U.tilde_beta = -t(c(p.tilde1*(1-p.tilde1))*x.tilde)%*%x.tilde
      v.mtx1 = t(solve(U.tilde_beta)%*%
                   t(c(samp[,y]-p.tilde1)*x.tilde))
      v.mtx10 = v.mtx1[samp[,sub]==1,]
      f1 = c(t(1+eta%*%t(v.mtx10)))
      c(f1*samp0$wt0)%*%v.mtx10-colSums(v.mtx1)
    }
    #S_beta.tilde.fun(x=beta.tilde)
    S_beta.tilde = jacobian(func=S_beta.tilde.fun, x = beta.tilde)
    
    U.tilde_beta = matrix(0,length(beta.tilde),length(beta.est))
    U.tilde_eta  = matrix(0,length(beta.tilde),length(beta.tilde))
    
    US.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%U_eta%*%solve(S_eta)),
                   cbind(S_beta, solve(S_eta)))
    Phi_theta.inv = rbind(cbind(US.inv, -US.inv%*%rbind(U_beta.tilde, S_beta.tilde)%*%solve(U.tilde_beta.tilde)),
                          cbind(U.tilde_beta, U.tilde_eta, solve(U.tilde_beta.tilde))
    )
    beta_w = -t(Phi_theta.inv%*%t(Phi_w))
    
    US_w = cbind(U_w, S_w)
    beta_w1 = -t(US.inv%*%t(US_w))
    
    TD.dat = as.data.frame(cbind(1, samp$psu, samp$strata, beta_w, beta_w1))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.all = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[1:ncol(beta_w)+3], collapse="+"))), ds.TD)))
    var.all1 = diag(vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[1:ncol(beta_w1)+3+ncol(beta_w)], collapse="+"))), ds.TD)))
  }
  if(structure=="subsample"){
    samp0$wt0 = samp0$wt/samp$wt[samp[,sub]==1]
    U_w.1 = matrix(0,nrow(samp), ncol(X))
    U_w.1[samp[,sub]==1,] = samp0$wt0*(u-v.mtx0%*%B)
    U_w.1 = v.mtx%*%B+U_w.1
    TD.dat = as.data.frame(cbind(samp$wt, samp$psu, samp$strata, U_w.1))
    names(TD.dat)=c("wt", "psu", "strata", paste0("V", 1:(ncol(TD.dat)-3)))
    ds.TD = svydesign(ids=~psu, strata = ~strata, weights=~wt, data=TD.dat)
    var.Phi.1 = vcov(svytotal(as.formula(paste0("~", paste(names(TD.dat)[-c(1:3)], collapse="+"))), ds.TD))
    #var.Phi.2 = t(c((samp0$wt0-1)*samp$wt[samp[,sub]==1]*samp0$wt)*cbind(u,v.mtx0))%*%cbind(u,v.mtx0)
    #var.all = diag(Phi_theta.inv%*%(var.Phi.1+var.Phi.2)%*%t(Phi_theta.inv))
    var.all = diag(solve(U_beta)%*%(var.Phi.1)%*%t(solve(U_beta)))
    var.all1=NA
  }
  
  return(list(beta.est = beta.est, 
              var.all  = var.all,
              var.all1 = var.all1))
}





############################################################################################
beta.est.var.ex = function(ds.1, fit.y, fit.i, beta.star, beta.star.var, b=1, c=NULL, nest = T){
  vars = all.vars(as.formula(fit.i))
  y = vars[1]
  samp = ds.1$variables
  samp$wt = weights(ds.1)
  x.mtx = model.matrix(lm(fit.i, samp))
  n.x = ncol(x.mtx)
  colnames(x.mtx)[1] = "x0"
  p.tilde = 1/(1+exp(-x.mtx%*%c(beta.star)))
  v.mtx = c(samp[,y]-p.tilde)*x.mtx
  colnames(v.mtx) = colnames(x.mtx)
  ds.1 = update(ds.1, v.mtx = v.mtx)
  greg.out = greg_f(wt0=samp$wt, v.mtx=v.mtx, VS.hat=rep(0, ncol(v.mtx)))#, Vs0.hat=Vs0.hat)
  f = greg.out$f
  f[f<0]=0
  samp$wt.c = samp$wt*f
  ds.c   = svydesign(ids=~1, strata = NULL, weights=~wt.c, data=samp)
  #calib = survey::calibrate(design = ds.1,
  #  formula = ~v.mtx-1,
  #                          population=c(0, ncol(v.mtx)), calfun = "linear")
  # final outcome model
  mdl = svyglm(fit.y, ds.c, family="binomial")
  svybeta.est = mdl$coefficients
  X = model.matrix(mdl)
  n.X = ncol(X)
  
  p = mdl$fitted.values
  U_beta = -t(c(samp$wt.c*p*(1-p))*X)%*%X
  U_eta = t(t(c(samp$wt*c(samp[,y]-p))*v.mtx)%*%X)
  S_beta = matrix(0,length(svybeta.est),length(svybeta.est))
  S_eta = t(c(samp$wt)*v.mtx)%*%v.mtx
  US.inv = rbind(cbind(solve(U_beta), -solve(U_beta)%*%U_eta%*%solve(S_eta)),
                 cbind(S_beta, solve(S_eta)))
  
  Ustar_beta.star = -t(c(samp$wt*(p.tilde*(1-p.tilde)))*x.mtx)%*%x.mtx
  var.U.star = Ustar_beta.star%*%beta.star.var%*%t(Ustar_beta.star)
  
  U_w = f*c(samp[,y]-p)*X
  
  if(nest){
    S_w = (f-1/c)*v.mtx
    var.U.star1 = vcov(svytotal(~v.mtx, ds.1))/c^2
    var.U.star = var.U.star-var.U.star1
  }else{
    S_w = f*v.mtx
  }
  ds.1 = update(ds.1,U_w=U_w, S_w = S_w)
  V_US = vcov(svytotal(~U_w+S_w, ds.1))+rbind(matrix(0,n.X, (n.x+n.X)),
                                              cbind(matrix(0, n.x, n.X), var.U.star))
  v.all = diag(US.inv%*%V_US%*%t(US.inv))
  return(list(beta.est = svybeta.est,
              var.all = v.all))
  
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

