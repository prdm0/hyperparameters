# ==============================================================================
#        "family.elliptical" set of functions used by 'elliptical'
# ==============================================================================
library(Matrix)
elliptical.deriv<-structure(
  .Data=list(
    g0=function(z,...) log(1/(sqrt(2*pi))*exp(-0.5*z^2)),
    g1=function(z,...) rep(-0.5,length(z)),
    g2=function(z,...) 1/4,
    g3=function(z,...) 3/4,
    g4=function(z,...) 1,
    g5=function(z,...) rep(0,length(z)),		   
    g0=function(z,...) log((1/pi)*(1+z^2)^(-1)),
    g1=function(z,...) -1/(1+z^2),
    g2=function(z,...) 1/8,
    g3=function(z,...) 3/8,
    g4=function(z,...) 1, ##non exist
    g5=function(z,...) 1/((1+z^2)^2),		   
    g0=function(z,df,...) log(((gamma(0.5*(df+1)))/((pi*df)^0.5*
                                                      gamma(0.5*df)))*(1+z^2/df)^(-0.5*(df+1))),
    g1=function(z,df,...) (df+1)/(-2*(df+z^2)),
    g2=function(z,df,...) (df+1)/(4*(df+3)),
    g3=function(z,df,...)(3*(df+1))/(4*(df+3)),
    g4=function(z,df,...) df/(df-2),
    g5=function(z,df,...) (df+1)/(2*(df+z^2)^2),
    g0=function(z,s,r,...) log((s^(r/2)*gamma(0.5+r/2))/(gamma(0.5)*
                                                           gamma(r/2))*(s+z^2)^(-0.5*(r+1))),
    g1=function(z,s,r,...) (r+1)/(-2*(s+z^2)),
    g2=function(z,s,r,...) (r*(r+1))/(4*s*(r+3)),
    g3=function(z,s,r,...) (3*(r+1))/(4*(r+3)),
    g4=function(z,s,r,...) s/(r-2),
    g5=function(z,s,r,...) (r+1)/(2*(s+z^2)^2),
    g0=function(z,...) log(1.4843300029*exp(-z^2)/(1+exp(-z^2))^2),
    g1=function(z,...) -tanh(z^2/2),
    g2=function(z,...) 1.477240176/4,
    g3=function(z,...) 4.013783934/4,
    g4=function(z,...) 0.79569,
    g5=function(z,...) -0.5+0.5*tanh(0.5*z^2)^2,
    g0=function(z,...) log(exp(z)/(1+exp(z))^2),
    g1=function(z,...) (exp(z)-1)/(-2*(z*(1+exp(z)))),
    g2=function(z,...) 1/12,
    g3=function(z,...) 2.42996/4,
    g4=function(z,...) pi^2/3,
    g5=function(z,...) (2*exp(z)*z-exp(2*z)+1)/(4*z^3*(1+exp(z)^2)),
    g0=function(z,alpha,mp,...) log((alpha*gamma(mp+mp))/(gamma(mp)*gamma(mp))*
                                      (exp(alpha*z)/(1+exp(alpha*z))^2)^mp),
    g1=function(z,alpha,mp,...) alpha*mp*(exp(alpha*z)-1)/(-2*(z*(1+exp(alpha*z)))), 
    g2=function(z,alpha,mp,...) (alpha^2*mp^2)/(4*(2*mp+1)),
    g3=function(z,alpha,mp,...) (2*m)*(2+mp^2*trigamma(mp))/(4*(2*mp+1)),		   
    g4=function(z,alpha,mp,...) 2*trigamma(mp),
    g5=function(z,alpha,mp,...) alpha*mp*(2*alpha*exp(alpha*z)*z-exp(2*alpha*z)+1)/(4*z^3*(1+exp(alpha*z)^2)),
    g0=function(z,epsi,sigmap,...) log( (1-epsi)*1/(sqrt(2*pi))*exp(-0.5*z^2)+
                                          epsi*1/(sqrt(2*pi)*sigmap)*exp(-0.5*z^2/sigmap^2)),
    g1=function(z,epsi,sigmap,...)((1-epsi)*exp(-z^2/2)+
                                     (epsi*(sigmap^2)^(-1.5)*exp(-z^2/(2*sigmap^2))))/
      ((-2)*((1-epsi)*exp(-z^2/2)+
               (epsi*(sigmap^2)^(-0.5)*exp(-z^2/(2*sigmap^2))))),
    g2=function(z,epsi,sigmap,...)
    {
      NULL
    },
    g3=function(z,epsi,sigmap,...) NULL, 
    g4=function(z,epsi,sigmap,...) 1+epsi*(sigmap^2-1),
    g5=function(z,epsi,sigmap,...)
    {
      NULL
    },
    g0=function(z,k,...) log(1/(gamma(1+((1+k)/2))*2^(1+(1+k)/2))*exp(-0.5*(abs(z)^(2/(1+k))))),		   
    g1=function(z,k,...) 1/(-2*(1+k)*(z^2)^(k/(1+k))),
    g2=function(z,k,...) (gamma((3-k)/2))/(4*(2^(k-1)*(1+k)^2*gamma((k+1)/2))),
    g3=function(z,k,...)(k+3)/(4*(k+1)) ,
    g4=function(z,k,...) 2^(1+k)*(gamma(1.5*(k+1))/(gamma((k+1)/2))),
    g5=function(z,k,...) k/(2*(z^2)^((2*k+1)/(1+k))*((1+k)^2))
  ),
  .Dim=c(6,9),
  .Dimnames=list(c("g0","g1","g2","g3","g4","g5"),
                 c("Normal","Cauchy","Student","Gstudent",
                   "LogisI","LogisII","Glogis",
                   "Cnormal","Powerexp")))


Student <- function(df=stop("no df argument"))
{
  if(df<0)
    stop(paste("allowed values for degrees of freedom positive"))
  make.family.elliptical("Student", arg=df)
}

Normal <- function()
{
  make.family.elliptical("Normal")
}
Cauchy <- function()
{
  make.family.elliptical("Cauchy")
}
Gstudent <- function(parm=stop("no s or r argument"))
{
  if((parm[1]<=0)||(parm[2]<=0))
    stop(paste("s and r must be positive"))
  
  make.family.elliptical("Gstudent",arg=list(s=parm[1],r=parm[2]))
}


LogisI <- function()
{
  make.family.elliptical("LogisI")
}
LogisII <- function()
{
  make.family.elliptical("LogisII")
}
Glogis <- function(parma=stop("no alpha=alpha(m) or m argument"))
{
  if((parma[1]<=0)||(parma[2]<=0))
    stop(paste("alpha=alpha(m) and m must be positive"))
  make.family.elliptical("Glogis",arg=list(alpha=parma[1],mp=parma[2]))
}
Cnormal <- function(parmt=stop("no alpha or epsi argument"))
{
  stop(paste("not implement yet"))
  if((parmt[1]<0)||(parmt[1]>1)||(parmt[2]<=0))
    stop(paste("0<=epsilon<=1 and sigma must be positive"))
  
  make.family.elliptical("Cnormal",arg=list(epsi=parmt[1],sigmap=parmt[2]))
}

Powerexp <- function(k=stop("no k argument"))
{
  if(abs(k)>1)
    stop(paste("k must be (-1,1)"))
  make.family.elliptical("Powerexp", arg=k)
}

make.family.elliptical <- function(name, arg, ...)
{
  if(is.character(name) && charmatch(name, dimnames(elliptical.deriv)[[2]], F))
  { 
    g0 <- elliptical.deriv[["g0",name]]
    g1 <- elliptical.deriv[["g1",name]]
    g2 <- elliptical.deriv[["g2",name]]
    g3 <- elliptical.deriv[["g3",name]]
    g4 <- elliptical.deriv[["g4",name]]
    g5 <- elliptical.deriv[["g5",name]]
    
  }
  else 
  {
    obj.deriv <- eval(parse(text=paste(name, ".deriv", sep="")))
    g0 <- obj.deriv[["g0",name]]
    g1 <- obj.deriv[["g1",name]]
    g2 <- obj.deriv[["g2",name]]
    g3 <- obj.deriv[["g3",name]]
    g4 <- obj.deriv[["g4",name]]
    g5 <- obj.deriv[["g5",name]]
  }
  family <- list(g0=g0,g1=g1,g2=g2,g3=g3,g4=g4,g5=g5,
                 df=if(charmatch(name, "Student", F)) arg,
                 s=if(charmatch(name, "Gstudent", F)) arg$s,
                 r=if(charmatch(name, "Gstudent", F)) arg$r,
                 alpha=if(charmatch(name, "Glogis", F)) arg$alpha ,
                 mp=if(charmatch(name, "Glogis", F)) arg$m,
                 epsi=if(charmatch(name, "Cnormal", F)) arg$epsi,
                 sigmap=if(charmatch(name, "Cnormal", F)) arg$sigmap,
                 k=if(charmatch(name, "Powerexp", F)) arg) 
  names(family) <- c("g0","g1","g2","g3","g4","g5",
                     "df","s","r","alpha","mp","epsi","sigmap","k")  
  structure(.Data = c(list(family=name), family), class=c("family.elliptical","family"))
}


family.elliptical <- function(elliptical.object)  
{
  if( length(elliptical.object$call$family) > 1 )
    eval(elliptical.object$call$family)
  else
    eval(parse(text=paste(elliptical.object$call$family,"()")))
}


# ==============================================================================
#                 Fitting routines for regression-elliptical models
# ==============================================================================


elliptical <- function( formula = formula(data),
                        DerB=NULL,
                        parmB=NULL,
                        family = Normal, 
                        data = sys.parent(),
                        dispersion = NULL,   
                        weights,                    
                        subset,
                        na.action = "na.fail", 
                        method = "elliptical.fit",
                        control = glm.control(epsilon=0.0001,maxit=100, trace=F),
                        model = F, x = F, y = T,
                        contrasts = NULL,
                        linear=T,
                        restrict=F,
                        Cres=NULL,
                        sol=NULL,
                        offset,				 
                        ...)
{
  call <- match.call()
  
  dist <- as.character(call$family)[1]
  user.def <- F
  if( charmatch(dist, c("Normal","Cauchy","Student","Gstudent","LogisI","LogisII",
                        "Glogis","Cnormal","Powerexp"), nomatch=F) ) 
    dist <- match.arg(dist, c("Normal","Cauchy","Student","Gstudent","LogisI","LogisII",
                              "Glogis","Cnormal","Powerexp"))
  else user.def <- T
  if(!charmatch(method,c("model.frame","elliptical.fit"),F))
    stop(paste("\n unimplemented method:", method))
  if(is.null(DerB)&& linear==F)	stop(paste("\n necessary derivative matrix:")	)	
  m <- match.call(expand=F)
  m$family <- m$method <- m$control <- m$model <- m$dispersion <- m$x <- 
    m$y <- m$contrasts <-m$linear<-m$restrict<-m$Cres<-m$sol<-m$DerB<-m$offset<-m$parmB<-m$... <- NULL   
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  if(method=="model.frame")      
    return(m)
  
  {
    if(!missing(family) && 
       !charmatch(dist, c("Normal","Cauchy","Student","Gstudent","LogisI","LogisII",
                          "Glogis","Cnormal","Powerexp"),F))
      cat(paste("\n work with user-defined family:", call$family, "\n")) 
    }
  
  if(!missing(dispersion) && is.number(dispersion) && !(dispersion > 0))
    stop("\n no negative values for dispersion parameter")
  
  Terms <- attr(m, "terms")
  Y <- model.extract(m, response) 
  if(!is.numeric(Y))
    stop("\n response must be numeric")
  
  X <- model.matrix(Terms, m, contrasts)   
  if(!is.numeric(X))
    stop("\n model matrix must be numeric")
  offset <- model.extract(m, offset)
  nobs <- nrow(X)
  if(length(offset)==1 && offset==0)
    offset <- rep(0,nobs)
  w <- model.extract(m, weights) ; wzero <- rep(F,nrow(m))
  if(!length(w))
    w <- rep(1,nrow(m))
  else if(any(w < 0))
    stop("\n negative weights not allowed")
  else 
  {
    wzero <- (w == 0)
    Y.org <- Y ; X.org <- X ; offset.org <- offset
    Y <- Y*w ; X <- diag(c(w))%*%X ; offset <- w*offset 
    if(any(wzero))
    {
      wpos <- !wzero
      fitted <- resid <- q1 <- q2 <- Y.org
      Y <- Y[wpos] ; X <- as.matrix(X[wpos,]) ; offset <- offset[wpos] 
      
      
    }        
  }
  
  
  method <- "elliptical.fit"
  
  elliptical.fitter <- get(method)
  
  
  
  
  
  
  offset4fit <- offset
  if(linear)
  {
    if(restrict)
    {
      fit <- elliptical.fitter(X=X, Y=Y,DerB=NULL,parmB=NULL, linear=T,offset=offset4fit, family=family, 
                               dispersion=dispersion,Cres=Cres,sol=sol,restrict=T,
                               maxit=control$maxit, epsilon=control$epsilon, 
                               trace=control$trace, ...)
    }
    else	{
      fit <- elliptical.fitter(X=X, Y=Y, DerB=NULL,parmB=NULL,linear=T,offset=offset4fit, family=family, 
                               dispersion=dispersion, restrict=F,
                               maxit=control$maxit, epsilon=control$epsilon, 
                               trace=control$trace, ...)
    }
    
  }
  else
  {
    if(restrict)
    {
      fit <- elliptical.fitter(X=X, Y=Y,DerB=DerB,parmB=parmB,linear=F,offset=offset4fit, family=family, 
                               dispersion=dispersion,Cres=Cres,sol=sol,restrict=T,
                               maxit=control$maxit, epsilon=control$epsilon, 
                               trace=control$trace, ...)
    }
    else	{
      fit <- elliptical.fitter(X=X, Y=Y,DerB=DerB,parmB=parmB,linear=F,offset=offset4fit, family=family, 
                               dispersion=dispersion,restrict=F,
                               maxit=control$maxit, epsilon=control$epsilon, 
                               trace=control$trace, ...)
    }
  }
  if(any(wzero))
  {
    nas <- is.na(fit$coef)
    
    fitted[wpos] <- fit$fitted.values/w[wpos]
    fitted[wzero] <- X.org[wzero,!nas]%*%as.vector(fit$coef[!nas]) + 
      if(length(offset.org) > 1) 
        offset.org[wzero] else 0
    fit$fitted.values <- fitted
    
    resid[wpos] <- fit$resid
    resid[wzero] <- (Y.org[wzero]-fitted[wzero])/sqrt(fit$dispersion)
    fit$residuals <- resid         
    
    q1[wpos] <- fit$q1 ; q2[wpos] <- fit$q2
    q1[wzero] <- family$g1(resid[wzero],df=family$df,alpha=family$alpha,
                           mp=family$mp,epsi=family$epsi,sigmap=family$sigmap,
                           k=family$k)
    q2[wzero] <- -2*q1[wzero]
    fit$q1 <- q1 ; fit$q2 <- q2
  }
  else 
    fit$fitted.values <- fit$fitted.values/w
  fit$weights <- w
  
  names(fit$fitted.values) <- names(fit$residuals) <- names(fit$q1) <- 
    names(fit$q2) <- NULL
  
  p <- dim(X)[2]  
  rank <- fit$rank
  df.residuals <- length(if(exists("X.org",frame=sys.nframe())) Y.org 
                         else Y)-rank-sum(w==0) + ifelse(restrict,nrow(Cres),0)
  asgn <- attr(if(exists("X.org",frame=sys.nframe())) X.org else X, 
               "assign") 
  if(rank < p)
  {  
    nas <- is.na(fit$coef)
    pasgn <- asgn[!nas]
    if(df.residuals > 0)
      fit$assign.residual <- (rank+1):length(Y)
    fit$R.assign <- pasgn
    fit$x.assign <- asgn
  }
  
  
  fit <- c(fit, list( assign = asgn,
                      df.residuals = df.residuals,
                      family = family,
                      user.def = user.def, 
                      formula = as.vector(attr(Terms,"formula")),
                      terms = Terms,
                      contrasts = attr(X,"contrasts"),
                      control=control,
                      call = call ))
  
  if(y) fit$y <- if(exists("Y.org",frame=sys.nframe())) Y.org else Y
  names(fit$y) <- NULL
  if(x) fit$X <- if(exists("X.org",frame=sys.nframe())) X.org else X
  if(model) fit$model <- m
  
  attr(fit,"class") <- c("elliptical","glm","lm")
  fit$restrict<-restrict
  if(restrict){
    fit$Cres<-Cres
    fit$sol<-sol
    
  }
  fit
}


elliptical.fit <- function(X, Y, DerB,parmB,linear,offset, family, dispersion, maxit, epsilon, trace,restrict,Cres,sol, ...)
{
  n<-nrow(X)
  if(is.null(offset)){offset<-rep(0,n)}
  if(linear)
  {
    p<-ncol(X)
    aux.model <- glm.fit(x=X,y=Y,offset=offset,family=gaussian())
    attr(aux.model,"class") <- c("glm","lm")  
    start <- aux.model$coef
  }
  else
  {
    pn<-length(parmB)
    start<-parmB
    names(start)<- dimnames(DerB(start,X)$gradient)[[2]]
    Xd<-DerB(start,X)$gradient
    aux.model <- glm.fit(x=Xd,y=Y,offset=offset,family=gaussian())
    attr(aux.model,"class") <- c("glm","lm")
  }
  
  is.null.disp <- is.null(dispersion)
  elliptical.disp <- !is.null.disp && !is.number(dispersion)   
  if( is.null.disp )
    dispersion <- (summary(aux.model)$dispersion)
  if( elliptical.disp )
    dispersion <-(summary(aux.model)$dispersion)
  
  args <- resid(aux.model)/sqrt(dispersion)			  
  if(linear){
    if(any(nas <- is.na(start)))    
    {
      names(nas)<- dimnames(X)[[2]]
      X <- X[,!nas]
      aux.model <- glm.fit(x=X,y=Y,offset=offset,family=gaussian())
      attr(aux.model,"class") <- c("glm","lm")
      start <- aux.model$coef
      dispersion<-(summary(aux.model)$dispersion)
    }
  }
  else
  {	 
    if(any(nas <- is.na(start)))  {
      names(nas)<- dimnames(DerB(start,X)$gradient)[[2]]
      Xd<-DerB(start,X)$gradient
      Xd <- Xd[,!nas]
      aux.model <- glm.fit(x=Xd,y=Y,offset=offset,family=gaussian())
      attr(aux.model,"class") <- c("glm","lm")
      start<-parmB
      dispersion<-(summary(aux.model)$dispersion)
    }
    
  }
  iter <- 1
  error2 <- error3 <- 0
  repeat 
  {
    if(trace)
      cat("\n iteration", iter, ":")
    
    
    {
      w.1 <- family$g1(args,df=family$df,r=family$r,
                       s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                       sigmap=family$sigmap,k=family$k)
      dg <- family$g2(args,df=family$df,r=family$r,
                      s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                      sigmap=family$sigmap,k=family$k)          
      fg <- family$g3(args,df=family$df,r=family$r,
                      s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                      sigmap=family$sigmap,k=family$k)
      
      
      if(linear) 
      {
        if(restrict){
          y.aux<-Y-offset
          w.h <- as.vector(-2*w.1)
          aux.model <- glm.fit(x=X,y=y.aux,weights=w.h,family=gaussian())
          attr(aux.model,"class") <- c("glm","lm")
          new.start <-matrix(coef(aux.model),ncol(X),1)+
            solve(as.matrix(t(X)%*%diag(w.h)%*%X))%*%t(Cres)%*%solve(as.matrix(
              Cres%*%solve(as.matrix(t(X)%*%diag(w.h)%*%X))%*%t(Cres)))%*%(
                sol-Cres%*%matrix(coef(aux.model),ncol(X),1))
          new.start<-as.numeric(new.start)
        }
        else
        {
          y.aux<-Y-offset
          w.h <- as.vector(-2*w.1)
          aux.model <- glm.fit(x=X,y=y.aux,weights=w.h,family=gaussian())
          attr(aux.model,"class") <- c("glm","lm")
          new.start <-coef(aux.model)
        }
      }
      else
      {			  
        if(restrict){
          y.aux<-Y-offset
          w.h<--2*w.1
          Xd<-DerB(start,X)$gradient
          mu<-DerB(start,X)$value
          beta<-matrix(start,pn,1)
          attr(aux.model,"class") <- c("glm","lm")
          KBB<-as.matrix(((4*dg)/dispersion)*(t(Xd)%*%(Xd)))
          UB<-(1/dispersion)*t(Xd)%*%diag(w.h)%*%matrix((y.aux-mu),n,1)
          new.start <-solve(KBB)%*%(KBB%*%beta+UB)+
            solve(KBB)%*%t(Cres)%*%solve(Cres%*%solve(KBB)%*%t(Cres))%*%(
              sol-Cres%*%solve(KBB)%*%(KBB%*%beta+UB))
          attr(aux.model,"class") <- c("glm","lm")			    
          new.start<-as.numeric(new.start)
        }
        else
        {
          y.aux<-Y-offset
          w.h<--2*w.1
          Xd<-DerB(start,X)$gradient
          mu<-DerB(start,X)$value
          beta<-matrix(start,pn,1)+(1/(4*dg))*
            solve(as.matrix(t(Xd)%*%Xd))%*%t(Xd)%*%diag(w.h)%*%matrix((y.aux-mu),n,1)
          attr(aux.model,"class") <- c("glm","lm")			    
          new.start<-as.numeric(beta)
        }
      } 
      }      
    error1 <- max(abs((new.start-start)/start))
    start<-new.start
    if(linear){
      abs.res <- Y-X%*%start-offset
    }
    else
    {
      abs.res <- Y-DerB(start,X)$value-offset
    }
    if(is.null.disp)
    {
      aux.dispersion <-dispersion
      new.dispersion <- mean((-2*w.1)*abs.res^2)
      error2 <- abs((new.dispersion-dispersion)/dispersion)
      dispersion <- new.dispersion
    }
    
    old.args <- args
    args <- abs.res/sqrt(dispersion) 
    
    if(trace)
    {
      loglik <- -0.5*length(abs.res)*log((dispersion)) + 
        sum(family$g0(abs.res/sqrt(dispersion), df=family$df, 
                      s=family$s, r=family$r,alpha=family$alpha, mp=family$mp,epsi=family$epsi,
                      sigmap=family$sigmap,k=family$k))
      cat(" log-likelihood =", signif(loglik,6))
    }
    error3 <- sqrt(sum((args-old.args)^2)/max(1e-20,sum(old.args^2)))
    if((iter == maxit) || (max(error1, error2, error3,na.rm = TRUE) < epsilon))
      break
    
    iter <- iter + 1
  }
  
  if(trace) cat("\n")
  
  if(maxit>1 && iter==maxit)
    warning(paste("\n linear convergence not obtained in", maxit,
                  "iterations"))
  
  
  coefs <- rep(NA,length(nas))
  coefs[!nas] <- start
  names(coefs) <- names(nas)	  
  names(dispersion) <- "dispersion"
  if(linear)	{
    fitted <- as.vector(X%*%start+offset)
  }
  else {
    fitted <- as.vector(DerB(start,X)$value+offset)
  }
  residuals <- (Y-fitted)/sqrt(dispersion)
  w.1 <- family$g1(residuals,df=family$df, 
                   s=family$s, r=family$r,alpha=family$alpha, mp=family$mp,epsi=family$epsi,
                   sigmap=family$sigmap,k=family$k)
  w.2 <- -2*w.1
  if( any(w.2<0) ) 
    cat("\n --- negative iterative weights returned! --- \n")
  
  
  if(is.null.disp)
  {
    if(linear) {
      rank <- dim(X)[2]
      Rnames <- dimnames(X)[[2]]
      Xd <- cbind(X,residuals)
    }
    else
    {
      rank<-dim(Xd)[2]
      Rnames <- dimnames(Xd)[[2]]	 
      Xd <- cbind(DerB(start,X)$gradient,residuals)
    }	 
  }
  dimnames(Xd)[[2]] <- c(Rnames,"scale")
  nn <- is.null(Rnames) 			
  Rnames <- list(dimnames(Xd)[[2]],dimnames(Xd)[[2]])
  R <- t(Xd)%*%Xd
  if(is.null.disp)
    R[rank+1,rank+1] <- R[rank+1,rank+1] + length(residuals)
  attributes(R) <- list(dim=dim(R))
  if(!nn) attr(R, "dimnames") <- Rnames
  
  loglik <- -0.5*length(residuals)*log((dispersion)) + 
    sum(family$g0(residuals, df=family$df, 
                  s=family$s, r=family$r,alpha=family$alpha, mp=family$mp,epsi=family$epsi,
                  sigmap=family$sigmap,k=family$k))
  names(loglik)<-NULL
  
  fit <- list(coefficients = coefs,
              dispersion = dispersion,
              fixed = !is.null.disp,
              residuals = residuals,
              fitted.values = fitted,
              loglik = loglik,
              Wg = family$g1(residuals,df=family$df,r=family$r,
                             s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                             sigmap=family$sigmap,k=family$k),
              Wgder = family$g5(residuals,df=family$df,r=family$r,
                                s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                                sigmap=family$sigmap,k=family$k), 
              v = -2*family$g1(residuals,df=family$df,r=family$r,
                               s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                               sigmap=family$sigmap,k=family$k), 
              rank = rank,
              R = as.matrix(R),
              iter = iter-1,
              scale=4*family$g2(residuals,df=family$df,r=family$r,
                                s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                                sigmap=family$sigmap,k=family$k),
              scaledispersion=-1+4*family$g3(args,df=family$df,r=family$r,
                                             s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                                             sigmap=family$sigmap,k=family$k), 
              scalevariance=family$g4(args,df=family$df,r=family$r,
                                      s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                                      sigmap=family$sigmap,k=family$k),
              df=if(charmatch(family$family, "Student", F)) family$df,
              s=if(charmatch(family$family, "Gstudent", F)) family$s,
              r=if(charmatch(family$family, "Gstudent", F)) family$r,
              alpha=if(charmatch(family$family, "Glogis", F))family$alpha ,
              mp=if(charmatch(family$family, "Glogis", F)) family$m,
              epsi=if(charmatch(family$family, "Cnormal", F)) family$epsi,
              sigmap=if(charmatch(family$family, "Cnormal", F)) family$sigmap,
              k=if(charmatch(family$family, "Powerexp", F)) family$k,
              Xmodel=matrix(Xd[,(1:rank)],nrow(Xd),rank),
              linear=linear,
              DerBB=if(linear) {NULL} else {DerB(start,X)$hessian	} ,
              nfunc=DerB
  )        
  
  fit
}
# ==============================================================================
#                              Methods for `summary'
# ==============================================================================


# for `elliptical' objects
# -----------------


summary.elliptical <- function(object, correlation = T)
{
  coef <- object$coef
  disp <- object$dispersion
  scale<-object$scale
  scaledispersion<-object$scaledispersion
  fixed <- object$fixed
  resid <- object$residuals   
  wt <- object$weights
  nas <- is.na(coef)
  n <- length(resid)	  - sum(wt==0)
  p <- object$rank
  Cres<-object$Cres
  restrict<-object$restrict    
  if(is.null(p))
    p <- sum(!nas)
  if(!p) 
  {
    warning("\n This model has zero rank --- no summary is provided")
    return(object)
  }
  rdf <- object$df.resid
  if(is.null(rdf))
    rdf <- n - p - sum(wt==0)+ifelse(object$restrict,nrow(object$Cres),0)
  
  R <- object$R[(1:p),(1:p)]
  Rnames <- dimnames(R)
  covun <- solve(qr(R))
  if(restrict){
    covun<-covun%*%(diag(p)-t(Cres)%*%solve(Cres%*%covun%*%t(Cres))%*%Cres%*%covun)}
  dimnames(covun) <- Rnames         
  rowlen <- sqrt(diag(covun))
  cnames <- names(coef[!nas])
  coef <- matrix(rep(coef[!nas], 4), ncol = 4)
  dimnames(coef) <- list(cnames, c("Value", "Std. Error", "z-value", "p-value"))
  coef[, 2] <- rowlen[1:p] %o% sqrt(disp/scale) 
  coef[, 3] <- coef[, 1]/coef[, 2]
  coef[, 4] <- 2 * pnorm( - abs(coef[, 3]) )
  
  if(!fixed)
  {
    disp <- matrix(c(disp,sqrt((4*disp^2)/(n*scaledispersion))),ncol=2)
    dimnames(disp) <- list("dispersion", c("Value", "Std. Error"))
  }
  
  if(correlation) 
  {
    correl <- covun * outer(1/rowlen, 1/rowlen)
    dimnames(correl) <- Rnames
  }
  else correl <- NULL
  
  summary <- list(coefficients = coef,
                  dispersion = disp,
                  fixed = fixed,
                  residuals = resid,
                  cov.unscaled = covun[(1:p),(1:p)],
                  correlation = correl[(1:p),(1:p)],
                  family = object$family,
                  loglik = object$loglik,
                  terms = object$terms,
                  df = c(p, rdf,n),
                  iter = object$iter, 
                  nas = nas,
                  call = object$call,
                  scale=scale,
                  scaledispersion=scaledispersion
  )
  attr(summary,"class") <- c("summary.elliptical")
  summary
}


# ==============================================================================
#                              Methods for `print'
# ==============================================================================


# for `family.elliptical' objects
# ------------------------


print.family.elliptical <- function(x, ...)
{
  cat("\n",x$family, "family\n")
  cat("\n density  : ", as.character(as.list(x[["g0"]])),"\n")
  cat("\n Wg: ", as.character(as.list(x[["g1"]])), "\n")
  cat("\n scale: ",  as.character(as.list(x[["g2"]])), "\n")
  cat("\n scale dispersion: ",  as.character(as.list(x[["g3"]])), "\n")
  cat("\n scale variance: ",  as.character(as.list(x[["g4"]])), "\n")
  cat("\n Wg': ",  as.character(as.list(x[["g5"]])), "\n")
  if(charmatch(x$family,"Gstudent",F)) 
    cat("\n r :", x$r, "\n","\n s :", x$s, "\n")
  if(charmatch(x$family,"Glogis",F)) 
    cat("\n alpha :", x$alpha, "\n","\n m :", x$mp, "\n")
  if(charmatch(x$family,"Cnormal",F)) 
    cat("\n epsilon :", x$epsi, "\n","\n sigma :", x$sigmap, "\n")
  if(charmatch(x$family,"Student",F)) 
    cat("\n df :", x$df, "\n")
  if(charmatch(x$family,"Powerexp",F)) 
    cat("\n k  :", x$k, "\n")
  return(invisible(0))
}


# for 'elliptical' objects
# -----------------


print.elliptical <- function(x, digits=6, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  coef <- x$coef
  if(any(nas <- is.na(coef))) {
    if(is.null(names(coef))) names(coef) <- paste("b", 1:length(
      coef), sep = "")        
    coef <- coef[!nas]
    cat("\nCoefficients: (", sum(nas), 
        " not defined because of singularities)\n", sep = "")
  }
  else cat("\nCoefficients:\n")
  print(coef, digits=digits, ...)
  cat("\nScale parameter: ", format(x$dispersion, digits=digits), 
      if(x$fixed) " (fixed)\n" else "\n")
  cat("\nError distribution: ", x$family[[1]], "\n")
  rank <- x$rank
  if(is.null(rank))
    rank <- sum(!nas)
  nobs <- length(x$residuals) - sum(x$weights==0)
  rdf <- x$df.resid
  if(is.null(rdf))
    rdf <- nobs - rank  +ifelse(x$restrict,nrow(x$Cres),0)
  cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")
  cat("-2*Log-Likelihood", format(-2 * x$loglik, digits=digits), "\n")
  invisible(x)
}


# for `summary.elliptical' objects
# -------------------------


print.summary.elliptical <- function(x, digits = 6, quote = T, prefix = "")
{
  nas <- x$nas
  p <- sum(!nas)
  coef <- x$coef
  correl <- x$correl
  if(any(nas)) {
    nc <- length(nas)
    cnames <- names(nas)
    coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
    coef1[!nas,  ] <- coef
    coef <- coef1
    if(!is.null(correl)) {
      correl1 <- matrix(NA, nc, nc, dimnames = list(cnames, 
                                                    cnames))
      correl1[!nas, !nas] <- correl[1:p,1:p]
      correl <- correl1
    }
  }
  if(is.null(digits))
    digits <- options()$digits
  else {
    old.digits <- options(digits = digits)
    on.exit(options(old.digits))
  }
  cat("Call: ")
  dput(x$call)
  if(any(nas))
    cat("\nCoefficients: (", sum(nas), 
        " not defined because of singularities)\n", sep = "")
  else cat("\nCoefficients:\n")
  print(coef, digits = digits)
  cat(paste("\nScale parameter for", x$family$family, ": "))
  cat(signif(x$dispersion[1], digits=digits), " (", 
      if(x$fixed) "fixed" else signif(x$dispersion[2], digits=digits), 
      ")\n")
  int <- attr(x$terms, "intercept")
  if(is.null(int))
    int <- 1
  df <- x$df
  nobs <- df[3] 
  cat("\nDegrees of Freedom:", nobs, "Total;", x$df[2], "Residual\n")
  cat("-2*Log-Likelihood", format(-2 * x$loglik), "\n")
  cat("\nNumber  Iterations:", format(trunc(x$iter)), 
      "\n")
  if(!is.null(correl)) {
    p <- dim(correl)[2]
    if(p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1,  - p, drop = F], quote = F, digits = 
              digits)
    }
  }
  invisible(x)
}
# ==============================================================================
#                     ANOVA for regression-elliptic models
# ==============================================================================

anova.elliptical <- function(object, ... , dispersion = NULL,test = c("Chisq"))
{
  test <- match.arg(test)
  
  margs <- function(...)  nargs()
  if(margs(...))
    return(anova.ellipticallist(list(object, ...), test = test))
  
  Terms <- object$terms
  term.labels <- attr(Terms, "term.labels")
  nt <- length(term.labels)
  m <- model.frame(object)
  family <- family(object)[[1]]
  y <- model.extract(m, "response")
  
  loglik <- double(nt + 1)
  df.res <- loglik
  if(nt) 
  {
    loglik[nt + 1] <- -2 * object$loglik
    df.res[nt + 1] <- object$df.residual
    fit <- object
    for(iterm in seq(from = nt, to = 1, by = -1)) 
    {
      ww <- fit$weights
      argslist <- list(object = fit,
                       formula = eval(parse(text = 
                                              paste("~ . -", term.labels[iterm]))))
      fit <- do.call("update", argslist)
      loglik[iterm] <- -2 * fit$loglik
      df.res[iterm] <- fit$df.residual
    }
    dev <- c(NA,  - diff(loglik))
    df <- c(NA,  - diff(df.res))
  }
  else 
  {
    loglik[1] <- -2 * object$loglik
    df.res[1] <- dim(y)[1] - attr(Terms, "intercept")
    dev <- df <- as.numeric(NA)     
  }
  
  heading <- c("Analysis of Deviance Table\n", 
               paste("Error distribution:", family), 
               paste("Response:", as.character(formula(object))[2], "\n", 
                     sep = ""), "Terms added sequentially (first to last)")
  if (is.null(dispersion)) {
    dispersion <- 1
    df.dispersion <- if (dispersion == 1) 
      Inf
    else object$df.residual
  }
  else df.scale <- Inf
  
  aod <- data.frame(Df = df, Deviance = dev, "Resid. Df" = df.res,
                    "-2*LL" = loglik, row.names = c("NULL", term.labels), 
                    check.names = F)
  attr(aod, "heading") <- heading
  class(aod) <- c("anova", "data.frame")
  if (is.null(test)) 
    return(aod)
  else
    aod<-stat.anova(aod, test, scale =1, df.scale = df.dispersion, n = nrow(y))
  structure(aod, heading = heading, class = c("anova", "data.frame"))
}

anova.ellipticallist <- function(object, ... , test = c("Chisq"))
{
  diff.term <- function(term.labels, i)
  {
    t1 <- term.labels[[1]]
    t2 <- term.labels[[2]]
    m1 <- match(t1, t2, F)
    m2 <- match(t2, t1, F)
    if(all(m1)) 
    {
      if(all(m2))
        return("=")
      else return(paste(c("", t2[ - m1]), collapse = "+"))
    }
    else 
    {
      if(all(m2))
        return(paste(c("", t1[ - m2]), collapse = "-"))
      else return(paste(i - 1, i, sep = " vs. "))
    }
  }
  test <- match.arg(test)
  rt <- length(object)
  if(rt == 1) 
  {
    object <- object[[1]]
    UseMethod("anova")
  }
  forms <- sapply(object, function(x)
    as.character(formula(x)))
  subs <- as.logical(match(forms[2,  ], forms[2, 1], F))
  if(!all(subs))
    warning("Some fit objects deleted because response differs from the first model" )
  if(sum(subs) == 1)
    stop("The first model has a different response from the rest")
  forms <- forms[, subs]
  object <- object[subs]
  dfres <- sapply(object, "[[", "df.resid")
  m2loglik <- -2 * sapply(object, "[[", "loglik")
  tl <- lapply(object, labels)
  rt <- length(m2loglik)
  effects <- character(rt)
  for(i in 2:rt)
    effects[i] <- diff.term(tl[c(i - 1, i)], i)
  dm2loglik <-  - diff(m2loglik)
  ddf <-  - diff(dfres)
  family <- family(object[[1]])[[1]]
  fixed <- object[[1]]$fixed
  heading <- c("Analysis of Deviance Table\n", 
               paste("Error distribution:", family), 
               paste("Response:", forms[2, 1], "\n"))
  nmodels <- length(object)
  varirest <- lapply(object, function(x) paste(deparse(x$call$restrict), 
                                               collapse = "\n"))
  
  variables <- lapply(object, function(x) paste(deparse(formula(x)), 
                                                collapse = "\n"))
  topnote <- paste("Model ", format(1:nmodels), ": ", variables,"     restrict =", varirest, 
                   sep = "", collapse = "\n")
  aod <- data.frame(Terms = forms[3,  ], "Resid. Df" = dfres, "-2*LL" = 
                      m2loglik, Test = effects, Df = c(NA, abs(ddf)), Deviance = c(NA, 
                                                                                   abs(dm2loglik)), check.names = F)
  attr(aod, "heading") <- heading
  
  
  if(!is.null(test)) 
  {
    n <- length(object[[1]]$residuals)
    o <- order(dfres)
    aod<- stat.anova(aod, test, 1, 
                     dfres[o[1]], n)
  }
  else aod
  structure(aod, heading = c(heading,topnote), class = c("anova", 
                                                         "data.frame"))
}

# ==============================================================================
#                    residuals for regression-elliptical  models
# ==============================================================================


residuals.elliptical<-
  function(object, type = c("stand", "pearson",  "response"))
  {
    type <- match.arg(type)
    rr <- switch(type,
                 pearson = object$resid/sqrt(object$scalevariance),
                 stand = {
                   Xd <- as.matrix(object$Xmodel)
                   Xdi <- solve(t(Xd) %*% Xd)
                   H <- Xd %*% Xdi %*% t(Xd)
                   H1 <- (1/(object$scalevariance * object$scale)) * H
                   varr <- object$scalevariance * object$dispersion * (1 - diag(H1))
                   ri <- object$y - object$fitted
                   ri/sqrt(varr)
                 },
                 response = object$y - object$fitted)
    if(is.null(object$na.action))
      rr
    else naresid(object$na.action, rr)
  }

# ==============================================================================
#                     Diagnostics for regression-elliptical  models
# ==============================================================================


elliptical.diag <- function(ellipticalfit, weighting = "observed")
{
  #
  #  Calculate diagnostics for objects of class "elliptical".  The diagnostics
  #  calculated are various types of residuals.
  if(is.null(ellipticalfit$DerBB)&& ellipticalfit$linear==F)	stop(paste("\n necessary second derivative matrix:")	)
  
  scalevariance<-ellipticalfit$scalevariance
  scale<-ellipticalfit$scale
  family <- ellipticalfit$family
  user.def <- ellipticalfit$user.def
  f.name <- family[[1]]
  dispersion <- ellipticalfit$dispersion
  w <- if(is.null(ellipticalfit$weights)) rep(1, length(ellipticalfit$residuals))
  else ellipticalfit$weights
  wzero <- (w == 0)
  resid <- ellipticalfit$residuals[!wzero]              # response residuals
  Xd<-diag(c(w[!wzero])) %*%ellipticalfit$Xmodel[!wzero,]
  dev <- 2 * ( family$g0(resid,df=family$df,r=family$r,
                         s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                         sigmap=family$sigmap,k=family$k) - 
                 family$g0(0,df=family$df,r=family$r,
                           s=family$s,alpha=family$alpha,mp=family$mp,epsi=family$epsi,
                           sigmap=family$sigmap,k=family$k)  )   
  p <- ellipticalfit$rank
  H<-Xd%*%solve(t(Xd)%*%Xd)%*%t(Xd)
  h<-diag(H)/(scalevariance*scale)
  rs <- resid/ sqrt(scalevariance*(1-h))    # standardized resid.
  ro <- ellipticalfit$y[!wzero]-ellipticalfit$fitted[!wzero]# ordinal residuals
  n<-length(resid)
  u<-ro^2/dispersion
  ct<-ellipticalfit$Wgder[!wzero]
  op<-ellipticalfit$Wg[!wzero]*ro
  a<-ellipticalfit$v[!wzero]-4*ct*u
  b<-op+(u*ct*ro)
  db<-diag(b)
  b<-matrix(b,n,1)
  u<-matrix(u,n,1)
  da<-diag(a)
  dai<-diag(1/a)
  ro<-matrix(ro,n,1)
  som<-matrix(0,p,p)
  if(!ellipticalfit$linear){
    if(p==1){ DerBB<-array(ellipticalfit$DerBB[,,!wzero],c(p,p,length(!wzero)))}
    else {	 DerBB<-ellipticalfit$DerBB[,,!wzero]  }
    for(i in 1:n)
    {
      som<-som+2*op[i]*DerBB[,,i]
    }
  } 
  M<-as.matrix(som+t(Xd)%*%da%*%Xd)
  lbb<-(-1/dispersion)*M
  lby<-(1/dispersion)*t(Xd)%*%da
  lphiy<-(-2/(dispersion^2))*t(b)
  lbphi<-(2/dispersion^2)*t(Xd)%*%b
  lphib<-t(lbphi)
  lphi<-(1/(dispersion^2))*((n/2)+t(u)%*%diag(ct)%*%u-(1/dispersion)*t(ro)%*%diag(ellipticalfit$v[!wzero])%*%ro)
  lbb1<-solve(lbb)
  lc1<-lbb1
  E<-as.vector(lphi-lphib%*%lbb1%*%lbphi)
  Fi<- -lbb1%*%lbphi
  GLbeta<-Xd%*%(-lbb1)%*%lby
  R<-Xd%*%(Fi%*%t(Fi)%*%lby+Fi%*%lphiy)
  GLphi<-(-1/E)*R
  G<-GLphir<-matrix(0,n,n)
  if(ellipticalfit$restrict){
    Cres<-ellipticalfit$Cres
    K<--lbb1%*%t(Cres)%*%solve(Cres%*%lbb1%*%t(Cres))%*%Cres%*%lbb1
    E<-E-lphib%*%K%*%lbphi
    G<-Xd%*%(-K)%*%lby
    GLphi<- (-1/E)*R
    Zr0<-lbb1%*%lbphi%*%lphib%*%t(K)
    Zr1<-K%*%lbphi%*%lphib%*%lbb1
    Zr2<-K%*%lbphi%*%lphib%*%t(K)
    Zr3<-K%*%lbphi%*%lphiy
    Zr<-Xd%*%((Zr0+Zr1+Zr2)%*%lby-Zr3)
    GLphir<-(-1/E)*Zr
    lc1<-lbb1+K
  }
  Om<-rbind(cbind(lbb,lbphi),cbind(lphib,lphi))
  Fr<- -lc1%*%lbphi
  lc1<-lc1+(1/E)*Fr%*%t(Fr)
  lc2<-(1/E)*Fr
  lc3<-t(lc2)
  lc4<-matrix(1/E,1,1)
  lc<-cbind(rbind(lc1,lc3),rbind(lc2,lc4))
  GLbeta<-diag(GLbeta)
  GLphi<-diag(GLphi)
  GLphir<-diag(GLphir)
  G<-diag(G)
  GL<-GLbeta+G+GLphi+GLphir
  Bi<-(a/dispersion)*GL
  ############pertubacao  wiLi
  
  deltab<-matrix(0,n,p)
  deltad<-matrix(0,n,1)
  deltab<-(1/dispersion)*diag((ellipticalfit$y[!wzero]-ellipticalfit$fitted[!wzero])*ellipticalfit$v[!wzero])%*%Xd
  deltad<-matrix(-(0.5/dispersion)*(1-ellipticalfit$v[!wzero]*u),n,1)
  delta<-t(cbind(deltab,deltad))
  b11<-cbind(matrix(0,p,p),matrix(0,p,1))
  b12<-cbind(matrix(0,1,p),1/E)
  b1<-rbind(b11,b12)
  b211<-cbind(lbb1,matrix(0,p,1))
  b212<-cbind(matrix(0,1,p),matrix(0,1,1))
  b2<-	rbind(b211,b212)
  Cic<--t(delta)%*%(lc)%*%delta
  Cic<-2*diag(Cic)	
  A<-as.matrix(t(delta)%*%(lc)%*%delta)
  decA<-eigen(A)
  Lmax<-decA$val[1]
  dmax<-decA$vec[,1]
  dmax<-dmax/sqrt(Lmax)
  dmaxc<-abs(dmax)
  ############pertubacao na escala (LD) modelo heterocedastico
  deltab<-(-2/dispersion)*t(Xd)%*%db
  deltad<-(-1/(dispersion^2))*t(ro)%*%db
  delta<-rbind(deltab,deltad)	   
  b11<-cbind(matrix(0,p,p),matrix(0,p,1))
  b12<-cbind(matrix(0,1,p),1/lphi)
  b1<-rbind(b11,b12)
  b211<-cbind(lbb1,matrix(0,p,1))
  b212<-cbind(matrix(0,1,p),0)
  b2<-	rbind(b211,b212)
  Ci<--t(delta)%*%(lc)%*%delta
  Cih<-2*diag(Ci)	
  A<-as.matrix(t(delta)%*%(lc)%*%delta)
  decA<-eigen(A)
  Lmax<-decA$val[1]
  dmax<-decA$vec[,1]
  dmax<-dmax/sqrt(Lmax)
  dmax<-abs(dmax)
  ############pertubacao na resposta(yi+wisi) baseado L-1 equivalente a GL
  deltab<-(1/dispersion)*t(Xd)%*%da
  deltad<-(-2/(dispersion^2))*t(b)
  delta<-rbind(deltab,deltad)	   
  b11<-cbind(matrix(0,p,p),matrix(0,p,1))
  b12<-cbind(matrix(0,1,p),1/lphi)
  b1<-rbind(b11,b12)
  b211<-cbind(lbb1,matrix(0,p,1))
  b212<-cbind(matrix(0,1,p),0)
  b2<-	rbind(b211,b212)
  Ci<--t(delta)%*%(lc)%*%delta
  Ci<-2*diag(Ci)	
  #########################################
  #################Pertubacao predicao na resposta
  ds<-diag(sqrt(dispersion),n)
  deltai<-(1/dispersion)*t(Xd)%*%da%*%ds
  Lmax<-NULL
  A<-matrix(0,n,n)
  for( i in 1:n){
    A[,i]<-as.matrix(t(deltai)%*%solve(t(Xd)%*%da%*%Xd)%*%matrix(Xd[i,],p,1))
  }
  Lmax<-abs(diag(A))
  ######################### Pertubacao predicao na var explicativa
  Cmax<-matrix(0,p,n)
  for( j in 1:p){
    Ff<-matrix(0,p,n)
    Ff[j,]<-rep(1,n)
    De<-diag(as.vector(ro),n)
    Dv<-diag(ellipticalfit$v[!wzero])
    st<-sqrt(var(Xd[,j]))
    for( i in 1:n){
      A[,i]<-st*(1/dispersion)*t(Ff%*%De%*%Dv-ellipticalfit$coef[j]*t(Xd)%*%da)%*%solve(t(Xd)%*%da%*%Xd)%*%matrix(Xd[i,],p,1)
      Cmax[j,i]<-2*abs(t(A[,i])%*%A[,i])
    }
  }
  
  list(resid = resid,rs=rs, dispersion = dispersion, GL=GL,GLbeta=GLbeta,GLphi=GLphi,
       G=G,GLphir=GLphir,dmax=dmax,Ci=Ci,delta=delta,Bi=Bi,Om=Om,IOm=lc,a=a,b=as.vector(b),c=ct,Cmax=Cmax,Lmax=Lmax,
       Cic=Cic,dmaxc=dmaxc,Cih=Cih,h=diag(H))
}

elliptical.diag.plots <- function(ellipticalfit, ellipticaldiag = NULL, weighting, which,
                                  subset = NULL, iden = F, labels = NULL, ret = F ,...)
{
  #  Diagnostic plots for objects of class "elliptical"
  
  if(is.null(ellipticaldiag))
  { 
    if(missing(weighting)) 
    { 
      family <- ellipticalfit$family
      user.def <- ellipticalfit$user.def
      f.name <- family[[1]]
      weighting <- "observed"
    }
    ellipticaldiag <- elliptical.diag(ellipticalfit, weighting = weighting)
  }
  if(is.null(subset))
    subset <- c(1:length(ellipticaldiag$GL))
  else if(is.logical(subset))
    subset <- (1:length(subset))[subset]
  else if(is.numeric(subset) && all(subset < 0))
    subset <- (1:(length(subset) + length(ellipticaldiag$GL)))[subset]
  else if(is.character(subset)) 
  {
    if(is.null(labels))
      labels <- subset
    subset <- seq(along = subset)
  }
  w <- if(is.null(ellipticalfit$weights)) rep(1, length(ellipticalfit$residuals))
  else ellipticalfit$weights
  wzero <- (w == 0)
  
  choices <- c("All", 
               "Response residual against fitted values",
               "Response residual against index",
               "Standardized residual against fitted values",
               "Standardized residual against index",
               "QQ-plot of  response residuals",
               "QQ-plot of  Standardized residuals",
               "Generalized Leverage ",
               "Ci against index ",
               "|Lmax| against index (local influence on coefficients)",
               "Bii against index\n")
  tmenu <- paste("plot:", choices)
  
  if( missing(which) )
    pick <- menu(tmenu, title = "\n Make a plot selection (or 0 to exit)\n")
  else if( !match(which, 2:11, nomatch=F) )
    stop("choice not valid") 
  else pick <- which
  
  if(pick == 0)
    stop(" no graph required ! ")
  
  #  old.par <- par()
  #    on.exit(par(old.par))
  
  repeat
  {
    switch(pick, 
           "1" = { par(pty="s", mfrow=c(1,1))
             
             close.screen(all = T)
             split.screen(c(2, 2))
             
             screen(1)       
             #  Plot the response residuals against the fitted values
             x1 <- ellipticalfit$fitted.values[!wzero]
             y1 <- ellipticaldiag$resid
             plot(x1,y1, xlab = "Fitted values", 
                  ylab = "Response residual", ...)
             
             screen(2)       #
             #  Plot the response residuals against the index values
             x2 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y2<-	 ellipticaldiag$resid
             plot(x2, y2, xlab = "Index", 
                  ylab = "Response residual", ...)
             
             screen(3)       #
             #  Plot the standardized residuals against the fitted values
             x3 <- ellipticalfit$fitted.values[!wzero]
             y3 <- ellipticaldiag$rs
             plot(x3, y3, xlab = "Fitted values", 
                  ylab = "Standardized residual", ...)
             
             screen(4)       #
             #  Plot the standardized residuals against the index values
             x4 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y4 <- ellipticaldiag$rs
             plot( x4,y4, xlab = "Index", 
                   ylab = "Standardized residual", ...)
             
             xx <- list(x1,x2,x3,x4)
             yy <- list(y1,y2,y3,y4)
             if(is.null(labels))
               labels <- names(model.extract(model.frame(ellipticalfit),
                                             "response"))
             
             yes <- iden
             while(yes) 
             {
               #  If interaction with the plots is required then ask the user which plot
               #  they wish to interact with and then run identify() on that plot.
               #  When the user terminates identify(), reprompt until no further interaction
               #  is required and the user inputs a 0.
               cat("****************************************************\n")
               cat("Please Input a screen number  (1,2,3, or 4)\n")
               cat("0 will terminate the function \n")
               num <- scan(n = 1)
               if((length(num) > 0) && ((num == 1) || (num == 2) || (num == 3) ||
                                        (num == 4))) {
                 cat(paste("Interactive Identification for screen", num, 
                           "\n"))
                 cat("left button = Identify, center button = Exit\n")
                 screen(num, new = F)
                 identify(xx[[num]], yy[[num]], labels, ...)
               }
               else yes <- F
             }
             close.screen(all=T)
             par(ask=T)
             split.screen(fig=c(2,2))
             
             screen(1)       #
             #  Plot a Normal QQ-plot of the response residuals
             y5 <- ellipticaldiag$resid
             x5 <- qnorm(ppoints(length(y5)))
             x5 <- x5[rank(y5)]
             .lim <- c(min(x5, y5), max(x5,y5))
             plot(x5, y5, xlab = paste("Quantiles of standard normal"), 
                  ylab = "Ordered response residual", 
                  xlim = .lim, ylim = .lim, ...)
             abline(0, 1, lty = 2)
             
             screen(2)       #
             #  Plot a Normal QQ-plot of the standardized residuals
             y6 <- ellipticaldiag$rs
             x6 <- qnorm(ppoints(length(y6)))	
             x6 <- x6[rank(y6)]
             .lim <- c(min(x6, y6), max(x6,y6))
             plot(x6, y6, xlab = paste("Quantiles of standard normal"), 
                  ylab = "Ordered standardized residual", 
                  xlim = .lim, ylim = .lim, ...)
             abline(0, 1, lty = 2)
             
             screen(3)       #
             
             #  Plot a leverage generalized									  
             y7 <- ellipticaldiag$GL
             x7 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x7, y7, xlab = "Index", 
                  ylab = "Generalized leverage ", ...)
             
             screen(4)       #
             #  Plot a Ci									  
             y8 <- ellipticaldiag$Ci
             x8 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x8, y8, xlab = "Index", 
                  ylab = "Ci", ...)
             
             xx <- list(x5, x6,x7,x8)
             yy <- list(y5, y6,y7,y8)
             if(is.null(labels))
               labels <- names(model.extract(model.frame(ellipticalfit),
                                             "response"))
             yes <- iden
             while(yes) 
             {
               #  If interaction with the plots is required then ask the user which plot
               #  they wish to interact with and then run identify() on that plot.
               #  When the user terminates identify(), reprompt until no further interaction
               #  is required and the user inputs a 0.
               cat("****************************************************\n")
               cat("Please Input a screen number  (1,2,3, or 4)\n")
               cat("0 will terminate the function \n")
               num <- scan(n = 1)
               if((length(num) > 0) && ((num == 1) || (num == 2) || (num == 3) ||
                                        (num == 4))) {
                 cat(paste("Interactive Identification for screen", num, 
                           "\n"))
                 cat("left button = Identify, center button = Exit\n")
                 screen(num, new = F)
                 identify(xx[[num]], yy[[num]], labels, ...)
               }
               else yes <- F
             }
             close.screen(all=T)
             par(ask=T) 
             split.screen(fig=c(2,2))
             
             screen(1)       #
             #  Plot a Dmax									  
             y9 <- ellipticaldiag$dmax
             x9 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x9, y9, xlab = "Index", 
                  ylab = "|Lmax|", ...)
             screen(2)       #
             #  Plot a Bi									  
             y10 <- ellipticaldiag$Bi
             x10 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x10, y10, xlab = "Index", 
                  ylab = "Bii", ...)
             
             #
             xx <- list(x9,x10)
             yy <- list(y9,y10)
             if(is.null(labels))
               labels <- names(model.extract(model.frame(ellipticalfit),
                                             "response"))
             
             yes <- iden
             while(yes) 
             {
               #  If interaction with the plots is required then ask the user which plot
               #  they wish to interact with and then run identify() on that plot.
               #  When the user terminates identify(), reprompt until no further interaction
               #  is required and the user inputs a 0.
               cat("****************************************************\n")
               cat("Please Input a screen number (1, 2 , 3 or 4)\n")
               cat("0 will terminate the function \n")
               num <- scan(n = 1)
               if((length(num) > 0) && ((num == 1) || (num == 2) || (num == 3)|| (num == 4))) {
                 cat(paste("Interactive Identification for screen", num, 
                           "\n"))
                 cat("left button = Identify, center button = Exit\n")
                 screen(num, new = F)
                 identify(xx[[num]], yy[[num]], labels, ...)
               }
               else yes <- F
             }                     
             close.screen(all = T) 
             par(ask=F)
           },
           "2" = { par(pty="s")
             #  Plot the response residuals against the fitted values
             x2 <- ellipticalfit$fitted.values[!wzero]
             y2 <-  ellipticaldiag$resid
             plot(x2,y2, xlab = "Fitted values", 
                  ylab = "Response residual", ...)
             xx <- list(x2)
             yy <- list(y2)
           },
           "3" = { par(pty="s")
             #  Plot the response residuals against the index
             x3 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y3 <-  ellipticaldiag$resid
             plot(x3,y3, xlab = "Index", 
                  ylab = "Response residual", ...)
             xx <- list(x3)
             yy <- list(y3)
           },					 
           "4" = { par(pty="s")
             #  Plot the Standardized residuals against the fitted values
             x4 <- ellipticalfit$fitted.values[!wzero]
             y4 <- ellipticaldiag$rs
             plot(x4, y4, xlab = "Fitted values", 
                  ylab = "Standardized residual", ...)
             xx <- list(x4)
             yy <- list(y4)
           },
           "5" = { par(pty="s")
             #  Plot the standardized residuals against the index
             x5 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y5 <-  ellipticaldiag$rs
             plot(x5,y5, xlab = "Index", 
                  ylab = "Standardized residual", ...)
             xx <- list(x5)
             yy <- list(y5)
           },
           "6" = { par(pty="s")
             #  Plot a Normal QQ-plot of the response residuals
             y6 <- ellipticaldiag$resid
             x6 <- qnorm(ppoints(length(y6)))
             x6 <- x6[rank(y6)] 
             .lim <- c(min(x6, y6), max(x6,y6))
             plot(x6, y6, xlab = "Quantiles of standard normal", 
                  ylab = "Ordered response residual", 
                  xlim = .lim, ylim = .lim, ...)
             abline(0, 1, lty = 2)
             xx <- list(x6)
             yy <- list(y6)
           },
           "7" = { 
             par(pty="s")
             #  Plot a Normal QQ-plot of the standardized residuals
             y7 <- ellipticaldiag$rs
             x7 <- qnorm(ppoints(length(y7)))[rank(y7)]
             .lim <- c(min(x7, y7), max(x7,y7))
             plot(x7, y7, xlab = paste("Quantiles of standard normal"), 
                  ylab = "Ordered standardized residual", 
                  xlim = .lim, ylim = .lim, ...)
             abline(0, 1, lty = 2)
             xx <- list(x7)
             yy <- list(y7)
           },
           "8"={
             par(pty="s")
             #  Plot a leverage generalized									  
             y8 <- ellipticaldiag$GL
             x8 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x8, y8, xlab = "Index", 
                  ylab = "Generalized leverage ", ...)
             xx <- list(x8)
             yy <- list(y8)
           },
           
           "9"={
             #  Plot a Ci
             par(pty="s")
             y9 <- ellipticaldiag$Ci
             x9 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x9, y9, xlab = "Index", 
                  ylab = "Ci", ...)
             xx <- list(x9)
             yy <- list(y9)
           },
           
           "10"={
             #  Plot a Dmax
             par(pty="s")
             y10 <- ellipticaldiag$dmax
             x10 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x10, y10, xlab = "Index", 
                  ylab = "|Lmax|", ...)
             xx <- list(x10)
             yy <- list(y10)
           },
           
           "11"={
             #  Plot a Bii
             par(pty="s")
             y11 <- ellipticaldiag$Bi
             x11 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x11, y11, xlab = "Index", 
                  ylab = "Bii", ...)
             xx <- list(x11)
             yy <- list(y11)
           }
           
    )
    
    if( (!(pick == 1)) )
    {           
      if(is.null(labels))
        labels <- names(model.extract(model.frame(ellipticalfit),"response"))
      yes <- iden
      while(yes) 
      {
        #  If interaction with the plots is required then ask the user which plot
        #  they wish to interact with and then run identify() on that plot.
        #  When the user terminates identify(), reprompt until no further interaction
        #  is required and the user inputs a 0.
        cat("****************************************************\n")
        cat("Interactive Identification\n")
        cat("left button = Identify, center button = Exit\n")
        identify(xx[[1]], yy[[1]], labels, ...)
        yes <- F
      }
    }
    
    if( missing(which) )
      pick <- menu(tmenu, 
                   title = "\n Make a plot selection (or 0 to exit)\n")
    if( (pick == 0) || !missing(which) )
    {
      invisible(close.screen(all=T))    
      break
    }
  }       
  
  if(ret)
    ellipticaldiag
  else invisible()
}

elliptical.diag.rest.plots <- function(ellipticalfit, ellipticaldiag = NULL, weighting, which,
                                       subset = NULL, iden = F, labels = NULL, ret = F ,...)
{
  #  Diagnostic plots for objects of class "elliptical"
  
  if(is.null(ellipticaldiag))
  { 
    if(missing(weighting)) 
    { 
      family <- ellipticalfit$family
      user.def <- ellipticalfit$user.def
      f.name <- family[[1]]
      weighting <- "observed"
    }
    ellipticaldiag <- elliptical.diag(ellipticalfit, weighting = weighting)
  }
  if(is.null(subset))
    subset <- c(1:length(ellipticaldiag$GL))
  else if(is.logical(subset))
    subset <- (1:length(subset))[subset]
  else if(is.numeric(subset) && all(subset < 0))
    subset <- (1:(length(subset) + length(ellipticaldiag$GL)))[subset]
  else if(is.character(subset)) 
  {
    if(is.null(labels))
      labels <- subset
    subset <- seq(along = subset)
  }
  w <- if(is.null(ellipticalfit$weights)) rep(1, length(ellipticalfit$residuals))
  else ellipticalfit$weights
  wzero <- (w == 0)
  
  choices <- c("All",
               "Leverage ",
               "H-Leverage",
               "G-Leverage",
               "Phi-Leverage",
               "Phirest-Leverage",
               "M-Leverage",
               "M-phi-Leverage",
               " ",
               " ",
               " ")
  tmenu <- paste("plot:", choices)
  
  if( missing(which) )
    pick <- menu(tmenu, title = "\n Make a plot selection (or 0 to exit)\n")
  else if( !match(which, 2:11, nomatch=F) )
    stop("choice not valid") 
  else pick <- which
  
  if(pick == 0)
    stop(" no graph required ! ")
  
  old.par <- par()
  on.exit(par(old.par))
  
  repeat
  {
    switch(pick, 
           "1" = { par(pty="s", mfrow=c(1,1))
             #
             close.screen(all = T)
             split.screen(c(2, 2))
             #
             screen(1)       #
             #  Plot dos leverage
             x1 <-1:length(ellipticalfit$fitted.values[!wzero])
             y1 <- ellipticaldiag$GL
             plot(x1,y1, xlab = "Index", 
                  ylab = "Generalized Leverage", ...)
             #
             screen(2)       #
             #  Plot the H-Leverage
             x2 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y2<-	 ellipticaldiag$GLbeta
             plot(x2, y2, xlab = "Index", 
                  ylab = "H-Leverage", ...)
             #
             screen(3)       #
             #  Plot the -G-Leverage
             x3 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y3 <- -ellipticaldiag$G
             plot(x3, y3, xlab = "Index", 
                  ylab = "G-Leverage", ...)
             #
             screen(4)       #
             #  Plot the Phi-leverage
             x4 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y4 <- ellipticaldiag$GLphi
             plot( x4,y4, xlab = "Index", 
                   ylab = "Phi-Leverage", ...)
             #
             xx <- list(x1,x2,x3,x4)
             yy <- list(y1,y2,y3,y4)
             if(is.null(labels))
               labels <- names(model.extract(model.frame(ellipticalfit),
                                             "response"))
             #
             yes <- iden
             while(yes) 
             {
               #  If interaction with the plots is required then ask the user which plot
               #  they wish to interact with and then run identify() on that plot.
               #  When the user terminates identify(), reprompt until no further interaction
               #  is required and the user inputs a 0.
               cat("****************************************************\n")
               cat("Please Input a screen number  (1,2,3, or 4)\n")
               cat("0 will terminate the function \n")
               num <- scan(n = 1)
               if((length(num) > 0) && ((num == 1) || (num == 2) || (num == 3) ||
                                        (num == 4))) {
                 cat(paste("Interactive Identification for screen", num, 
                           "\n"))
                 cat("left button = Identify, center button = Exit\n")
                 screen(num, new = F)
                 identify(xx[[num]], yy[[num]], labels, ...)
               }
               else yes <- F
             }
             close.screen(all=T)
             par(ask=T)
             split.screen(fig=c(2,2))
             
             screen(1)       #
             #  Plot -Phirest-Leverage
             x5 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y5 <- -ellipticaldiag$GLphir
             plot( x5,y5, xlab = "Index", 
                   ylab = "Phirestricted-Leverage", ...)
             #
             screen(2)       #
             #  Plot M-Leverage
             x6 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y6 <- ellipticaldiag$GLbeta+ellipticaldiag$G
             plot( x4,y4, xlab = "Index", 
                   ylab = "M-Leverage", ...)
             #
             screen(3)       #
             #
             #  Plot a Phi-leverage generalized									  
             y7 <- ellipticaldiag$GLphi+ellipticaldiag$Gphir
             x7 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x7, y7, xlab = "Index", 
                  ylab = "Phi-Generalized leverage ", ...)
             
             screen(4)       #
             #  Plot a Ci									  
             y8 <- ellipticaldiag$Ci
             x8 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x8, y8, xlab = "Index", 
                  ylab = "Ci", ...)
             
             xx <- list(x5, x6,x7,x8)
             yy <- list(y5, y6,y7,y8)
             if(is.null(labels))
               labels <- names(model.extract(model.frame(ellipticalfit),
                                             "response"))
             yes <- iden
             while(yes) 
             {
               #  If interaction with the plots is required then ask the user which plot
               #  they wish to interact with and then run identify() on that plot.
               #  When the user terminates identify(), reprompt until no further interaction
               #  is required and the user inputs a 0.
               cat("****************************************************\n")
               cat("Please Input a screen number  (1,2,3, or 4)\n")
               cat("0 will terminate the function \n")
               num <- scan(n = 1)
               if((length(num) > 0) && ((num == 1) || (num == 2) || (num == 3) ||
                                        (num == 4))) {
                 cat(paste("Interactive Identification for screen", num, 
                           "\n"))
                 cat("left button = Identify, center button = Exit\n")
                 screen(num, new = F)
                 identify(xx[[num]], yy[[num]], labels, ...)
               }
               else yes <- F
             }
             close.screen(all=T)
             par(ask=T)
             split.screen(fig=c(2,2))
             
             screen(1)       #
             #  Plot a Dmax									  
             y9 <- ellipticaldiag$dmax
             x9 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x9, y9, xlab = "Index", 
                  ylab = "|Lmax|", ...)
             screen(2)       #
             #  Plot a Bi									  
             y10 <- ellipticaldiag$Bi
             x10 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x10, y10, xlab = "Index", 
                  ylab = "Bii", ...)
             
             #
             xx <- list(x9,x10)
             yy <- list(y9,y10)
             if(is.null(labels))
               labels <- names(model.extract(model.frame(ellipticalfit),
                                             "response"))
             
             #
             yes <- iden
             while(yes) 
             {
               #  If interaction with the plots is required then ask the user which plot
               #  they wish to interact with and then run identify() on that plot.
               #  When the user terminates identify(), reprompt until no further interaction
               #  is required and the user inputs a 0.
               cat("****************************************************\n")
               cat("Please Input a screen number (1, 2 , 3 or 4)\n")
               cat("0 will terminate the function \n")
               num <- scan(n = 1)
               if((length(num) > 0) && ((num == 1) || (num == 2) || (num == 3)|| (num == 4))) {
                 cat(paste("Interactive Identification for screen", num, 
                           "\n"))
                 cat("left button = Identify, center button = Exit\n")
                 screen(num, new = F)
                 identify(xx[[num]], yy[[num]], labels, ...)
               }
               else yes <- F
             }
             #                       
             close.screen(all = T) 
             par(ask=F)
           },
           "2" = { par(pty="s")
             #  Plot the response residuals against the fitted values
             x2 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y2 <-  ellipticaldiag$GL
             plot(x2,y2, xlab = "Index", 
                  ylab = "Generalized Leverage", ...)
             xx <- list(x2)
             yy <- list(y2)
           },
           "3" = { par(pty="s")
             #  Plot the H-Leverage
             x3 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y3 <-  ellipticaldiag$GLbeta
             plot(x3,y3, xlab = "Index", 
                  ylab = "H-Leverage", ...)
             xx <- list(x3)
             yy <- list(y3)
           },					 
           "4" = { par(pty="s")
             #  Plot the -G-Leverage
             x4 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y4 <- -ellipticaldiag$G
             plot(x4, y4, xlab = "Index", 
                  ylab = "G-Leverage", ...)
             xx <- list(x4)
             yy <- list(y4)
           },
           "5" = { par(pty="s")
             #  Plot the Phi-leverage
             x5 <- 1:length(ellipticalfit$fitted.values[!wzero])
             y5 <-  ellipticaldiag$GLphi
             plot(x5,y5, xlab = "Index", 
                  ylab = "Phi-Leverage", ...)
             xx <- list(x5)
             yy <- list(y5)
           },
           "6" = { par(pty="s")
             #  Plot -Phirest-Leverage
             y6 <- -ellipticaldiag$GLphir
             x6 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot( x6,y6, xlab = "Index", 
                   ylab = "Phirestricted-Leverage", ...)
             xx <- list(x6)
             yy <- list(y6)
           },
           "7" = { 
             par(pty="s")
             #  Plot M-Leverage
             y7 <-ellipticaldiag$GLbeta+ellipticaldiag$G
             x7 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot( x7,y7, xlab = "Index", 
                   ylab = "M-Leverage", ...)
             #                       xx <- list(x7)
             yy <- list(y7)
           },
           "8"={
             par(pty="s")
             #  Plot a Mphi-Leverage									  
             y8 <- ellipticaldiag$GLphi+ellipticaldiag$GLphir
             x8 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x8, y8, xlab = "Index", 
                  ylab = "Mphi-leverage ", ...)
             xx <- list(x8)
             yy <- list(y8)
           },
           
           "9"={
             #  Plot a Ci
             par(pty="s")
             y9 <- ellipticaldiag$Ci
             x9 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x9, y9, xlab = "Index", 
                  ylab = "Ci", ...)
             xx <- list(x9)
             yy <- list(y9)
           },
           
           "10"={
             #  Plot a Dmax
             par(pty="s")
             y10 <- ellipticaldiag$dmax
             x10 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x10, y10, xlab = "Index", 
                  ylab = "|Lmax|", ...)
             xx <- list(x10)
             yy <- list(y10)
           },
           
           "11"={
             #  Plot a Bii
             par(pty="s")
             y11 <- ellipticaldiag$Bi
             x11 <- 1:length(ellipticalfit$fitted.values[!wzero])
             plot(x11, y11, xlab = "Index", 
                  ylab = "Bii", ...)
             xx <- list(x11)
             yy <- list(y11)
           }
           
    )
    
    if( (!(pick == 1)) )
    {           
      if(is.null(labels))
        labels <- names(model.extract(model.frame(ellipticalfit),"response"))
      yes <- iden
      while(yes) 
      {
        #  If interaction with the plots is required then ask the user which plot
        #  they wish to interact with and then run identify() on that plot.
        #  When the user terminates identify(), reprompt until no further interaction
        #  is required and the user inputs a 0.
        cat("****************************************************\n")
        cat("Interactive Identification\n")
        cat("left button = Identify, center button = Exit\n")
        identify(xx[[1]], yy[[1]], labels, ...)
        yes <- F
      }
    }
    
    if( missing(which) )
      pick <- menu(tmenu, 
                   title = "\n Make a plot selection (or 0 to exit)\n")
    if( (pick == 0) || !missing(which) )
    {
      invisible(close.screen(all=T))    
      break
    }
  }       
  
  if(ret)
    ellipticaldiag
  else invisible()
}
rms.curv.elliptical<-
  function(object, fit.val = get.fit.val(object, data), data = object$call$data,alpha=0.05,...)
  {
    get.fit.val <- function(object, data)
    {
      if(is.null(data))
        data <- sys.parent(2)
      if(is.numeric(data))
        data <- sys.frame(data)
      if(is.name(data))
        data <- get(data)
      class(data) <- NULL
      pp <- as.list(b <- coef(object))
      np <- names(b)
      data[np] <- pp
      eval(as.expression(object$formula[3]), local = data)
    }
    #	v <- attr(fit.val, "gradient")
    dist<-object$family[[1]]
    v<-object$Xmodel
    if(is.null(v))
      stop("gradient attribute missing")
    #	a <- attr(fit.val, "hessian")
    p <- ncol(v)
    n <- nrow(v)
    b<-object$DerBB
    a<-array(0,c(n,p,p))
    for(j in 1:p){
      for(k in 1:p){
        for(i in 1:n){		
          a[i,j,k] <- b[j,k,i]}}}
    if(is.null(a))
      stop("hessian attribute missing")
    if(charmatch(dist,"normal",F))
    {
      rho<- sqrt((sum((object$y-object$fitted)^2)/(n-p))*p)
    }
    else{
      rho <- sqrt(object$dispersion/object$scale)
    }
    D <- v
    for(j in 1:p)
      D <- cbind(D, a[, 1:j, j])
    qrd <- qr(D)
    Q <- qr.Q(qrd)
    rnk <- qrd$rank
    if(rnk <= p)
      warning("regression apparently linear")
    Q1 <- Q[, 1:rnk]
    C <- array(0, c(rnk, p, p))
    for(j in 1:p)
      C[,  , j] <- crossprod(Q1, a[,  , j])
    C <- aperm(C, c(2, 3, 1))
    r11i <- solve(qr.R(qrd)[1:p, 1:p])
    ct <- 0
    for(j in 1:p) {
      C[,  , j] <- crossprod(r11i, C[,  , j]) %*% r11i * rho
      ct <- ct + 2 * sum(C[,  , j]^2) + sum(diag(C[,  , j]))^2
    }
    ci <- 0
    for(j in (p + 1):rnk) {
      C[,  , j] <- crossprod(r11i, C[,  , j]) %*% r11i * rho
      ci <- ci + 2 * sum(C[,  , j]^2) + sum(diag(C[,  , j]))^2
    }
    ct <- sqrt(ct/(p * (p + 2)))
    ci <- sqrt(ci/(p * (p + 2)))
    if(charmatch(dist,"Normal",F))
    {
      pe <- ct * sqrt(qf((1-alpha), p, n - p))
      ic <- ci * sqrt(qf((1-alpha), p, n - p))
    } else
    {
      pe <- ct * sqrt(qchisq((1-alpha), p))
      ic <- ci * sqrt(qchisq((1-alpha), p))
    }
    val <- list(pe = pe, ic = ic, ct = ct, ci = ci, C = C,D=D,dist=dist)
    class(val) <- "rms.curv.elliptical"
    val
  }
print.rms.curv.elliptical<-function(x, ...)
{
  if(charmatch(x$dist,"normal",F))
  {
    
    cat("Parameter effects: c^theta x sqrt(F) =", round(x$pe, 4), "\n", 
        "      Intrinsic: c^iota  x sqrt(F) =", round(x$ic, 4), "\n", ...)
  }
  else
  {
    cat("Parameter effects: c^theta x sqrt(chisq) =", round(x$pe, 4), "\n", 
        "      Intrinsic: c^iota  x sqrt(chisq) =", round(x$ic, 4), "\n", ...)
  }
  invisible(x)
}

A.curv.elliptical<-function(object,...){
  gradient<-object$Xmodel
  if(is.null(gradient))
    stop("gradient attribute missing")
  hessian<-object$DerBB	
  n<-nrow(gradient)
  p<-ncol(gradient)
  dist<-object$family[[1]]
  if(charmatch(dist,"normal",F))
  {
    rho<- sqrt((sum((object$y-object$fitted)^2)/(n-p))*p)
  }
  else{
    rho <- sqrt(object$dispersion/object$scale)
  }	
  gradient.scaled<-gradient/rho
  qrstr<-qr(as.matrix(gradient.scaled))
  R<-qr.R(qrstr)
  Q<-qr.Q(qrstr,complete=T)
  K<-solve(R)
  hessian.scaled<-hessian/rho
  hessian.red<-hessian.scaled
  G<-array(0,dim=c(p,p,n))  
  A<-array(0,dim=c(p,p,n))
  for(i in 1:n)
  {
    G[,,i]<-t(K)%*%hessian.red[,,i]%*%K
  }
  aux<-matrix(0,n,1)
  for(i in 1:p)
  {
    for(j in 1:p){
      for(k in 1:n){
        aux[k,1]<-G[i,j,k]
      }
      A[i,j,]<-t(Q)%*%aux
    } 
  }
  A.pe<-A[,,1:p]
  A.ic<-A[,,(p+1):n]
  if(p==1) { A.pe<-array(A[,,1:p],dim=c(1,1,1))}
  if((p==1)&&((n-p)==1)){ A.ic<-array(A[,,(p+1):n],dim=c(1,1,1))}
  eta<-matrix(0,n,1)
  Xin<-solve(t(gradient)%*%gradient)
  for(i in 1:n){
    eta[i,1]<-(-0.5*(rho^2))*sum(diag(Xin%*%hessian[,,i]))
  }
  if(charmatch(dist,"normal",F))
  {
    eta<-eta/p
  }
  bias<-Xin%*%t(gradient)%*%eta
  pbias<-(100*as.vector(bias))/coef(object)
  curv<-list(A=A,A.pe=A.pe,A.ic=A.ic,K=K,bias=bias,pbias=pbias,eta=eta)
  curv
}
max.curv.elliptical<-function(object,maxit=100,...){
  gradient<-object$Xmodel
  n<-nrow(gradient)
  p<-ncol(gradient)
  dist<-object$family[[1]]
  if(charmatch(dist,"normal",F))
  {
    rho<- sqrt((sum((object$y-object$fitted)^2)/(n-p))*p)
  }
  else{
    rho <- sqrt(object$dispersion/object$scale)
  }
  
  A.curvt<-A.curv.elliptical(object)
  A.pe<-A.curvt$A.pe
  A.ic<-A.curvt$A.ic
  di<-matrix(0,p,1)
  di[p,1]<-1
  j<-0
  iter<-0
  eps<-0.0001
  diff<-NULL
  gi<-gif<-matrix(0,p,1)
  curv.pe<-NULL
  repeat {
    for(i in 1:p){
      gi<-gi+as.vector((t(di)%*%A.pe[,,i]%*%di))*(A.pe[,,i]%*%di)
    }
    gi<-4*gi
    gif<-matrix(as.vector(gi)/sqrt(sum(as.vector(gi)^2)),p,1)
    diff<-t(gif)%*%di
    if((iter == maxit) || (diff >= eps))
    {
      curv.pe<-0
      for(i in 1:p){
        curv.pe<-curv.pe+abs(t(di)%*%A.pe[,,i]%*%di)^2
      }
      curv.pe<-sqrt(curv.pe)
      di<-gif
      break 
    }
    iter <- iter + 1
    const<-3*gif+di
    const<-sqrt(sum(const))
    di<-gif/const
  }
  iter<-0
  di<-matrix(0,p,1)
  di[p,1]<-1
  eps<-0.0001
  diff<-NULL
  gi<-gif<-matrix(0,p,1)
  curv.ic<-NULL
  repeat{
    for(i in 1:(n-p)){
      gi<-gi+as.vector((t(di)%*%A.ic[,,i]%*%di))*(A.ic[,,i]%*%di)
    }
    gi<-4*gi
    gif<-matrix(as.vector(gi)/sqrt(sum(as.vector(gi)^2)),p,1)
    diff<-t(gif)%*%di
    if((iter == maxit) || (diff >= eps))
    {
      curv.ic<-0
      for(i in 1:(n-p)){
        curv.ic<-curv.ic+abs(t(di)%*%A.ic[,,i]%*%di)^2
      }
      curv.ic<-sqrt(curv.ic)
      di<-gif
      break 
    }
    iter <- iter + 1
    const<-3*gif+di
    const<-sqrt(sum(const))
    di<-gif/const
  }
  curv.el<-list(curv.pe=curv.pe,curv.ic=curv.ic,bias=A.curvt$bias,pbias=A.curvt$pbias)
  class(curv.el) <- "curv.elliptical"
  curv.el
}
print.curv.elliptical<-function(x, ...)
{
  cat("Parameter effects:  =", round(x$curv.pe, 5), "\n", 
      "Parameter Intrinsic:  =", round(x$curv.ic, 5), "\n", ...)
  cat("bias of location parameter =", round(x$bias,5),"\n")
  cat("% relative bias of location parameter =", round(x$pbias,5),"%","\n")
  
  
  invisible(x)
}


rpowerexp <- function(n,k=0.5)                          
{
  if( is.na(n) )
    return(NA)
  u<-runif(n,-1,1)
  r<-2/(1+k)
  ff<- rgamma(n,(1+1/r),1)
  (2*ff)^(1/r)*u
}

ppowerexp<- function(q,k=0.5)
{
  r<-2/(1+k)
  punif(q^r/2, -1, 1)*pgamma(q,(1+1/r),1)
}

dpowerexp <- function(x,k=0.5,mean=0,sd=1)
{
  z<-(x-mean)/sd
  1/(sd*(gamma(1+((1+k)/2))*2^(1+(1+k)/2)))*exp(-0.5*(abs(z)^(2/(1+k))))		   		 
}

rinvgamma<-function(n,s=1/2,r=1)
{
  if( is.na(n) )
    return(NA)
  1/(rgamma(n,s,r))
}

pinvgamma<-function(x,s=1/2,r=1)
{
  if( is.na(n) )
    return(NA)
  1-pgamma((1/x),s,r)
}

dinvgamma<-function(q,s=1/2,r=1)
{
  if( is.na(n) )
    return(NA)
  dgamma((1/q),s,r)*(1/x^2)
}

rgstudent<- function(n,s=1,r=2)                          
{
  if( is.na(n) )
    return(NA)
  z<-rnorm(n,0,1)
  v<-rinvgamma(n,r/2,s/2)
  v^(-0.5)*z
}

pgstudent<- function(q,s=1,r=2)
{
  (1- pinvgamma(1/sqrt(q), s/2,r/2 ))*pnorm(q,0,1)
}

dgstudent <- function(x,s=1,r=2,mean=0,sd=1)
{
  z<-(x-mean)/sd
  dens<-(1/sd)(s^(r/2)*gamma(0.5+r/2))/(gamma(0.5)*
                                          gamma(r/2))*(s+z^2)^(-0.5*(r+1))
}

dlogisI <- function(x,mean=0,sd=1)
{
  z<-(x-mean)/sd
  ((1.484300029*exp(-(z^2)))/(sd*(1+exp(-(z^2)))^2))
}

rlogisII <- function(n)                          
{
  if( is.na(n) )
    return(NA)
  u<-runif(n,0,1)
  log(u/(1-u))
}

plogisII<- function(q)
{
  punif((exp(q)/(1+exp(q))), 0, 1)
}

dlogisII <- function(x,mean=0,sd=1)
{
  z<-(x-mean)/sd
  (exp(z)/(sd*(1+exp(z))^2))
}

vari<-function(x){
  wnas <- x[!is.na(x)]
  var(x,na.rm=TRUE)*(length(wnas)-1)/length(wnas)
}










skewn<-function(x, na.rm = F, method = "fisher")
{
  method <- char.expand(method, c("fisher", "moment"), stop(
    "argument 'method' must match either \"fisher\" or \"moment\""))
  if(na.rm) {
    wnas <- x[!is.na(x)]
    if(length(wnas))
      x <- wnas
  }
  else if(any(is.na(x[!is.na(x)])))
    return(NA)
  n <- length(x)
  if(method == "fisher" && n < 3)
    return(NA)
  x <- x - mean(x)
  if(method == "moment")
    (sum(x^3)/n)/(sum(x^2)/n)^1.5
  else ((sqrt(n * (n - 1))/(n - 2)) * (sum(x^3)/n))/((sum(x^2)/n)^1.5)
}

kurt <-function(x, na.rm = F, method = "fisher")
{
  method <- char.expand(method, c("fisher", "moment"), stop(
    "argument 'method' must match either \"fisher\" or \"moment\""))
  if(na.rm) {
    wnas <- x[!is.na(x)]
    if(length(wnas))
      x <- wnas
  }
  else if(any(is.na(x[!is.na(x)])))
    return(NA)
  n <- length(x)
  if(method == "fisher" && n < 4)
    return(NA)
  x <- x - mean(x)
  if(method == "moment")
    (sum(x^4)/n)/(sum(x^2)/n)^2 - 3
  else ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 - (3 * (n - 1))/(n + 1)))/((n - 2) * (n - 3))
}

envelope<-function(object,DerB=NULL,B=100,arg=arg,...)
{
  if(object$linear==F) {
    initial<-object$coef
    DerB<-object$nfunc
  } else{ initial<-NULL}
  if(is.null(DerB)&& object$linear==F)	stop(paste("\n necessary derivate matrix:")	)	
  X <- model.matrix(object$terms)
  Xd<-as.matrix(object$Xmodel)
  n <- nrow(Xd)
  p <- ncol(Xd)
  ro <- object$resid
  tdf <- ro/sqrt(object$scalevariance)
  ####residual Cox and Snell
  Xdi<-solve(t(Xd)%*%Xd)
  H<-Xd%*%Xdi%*%t(Xd)
  H1<-(1/(object$scalevariance*object$scale))*H
  varr<-object$scalevariance*object$dispersion*(diag(1,n)-H1)
  varr<-diag(varr)
  ri<-object$y-object$fitted
  tdf<-ri/sqrt(varr)
  ############
  #
  e <- matrix(0,n,B)
  #
  mu<-object$fitted
  phi<-object$dispersion
  resp<-NULL
  for(i in 1:B){
    dist<-object$family[[1]]
    if(charmatch(dist,"Normal",F))
    {
      resp<-rnorm(n,0,1)
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=Normal(),control=glm.control(maxit=1000))
    }
    else if(charmatch(dist,"Cauchy",F))
    {
      resp<-rcauchy(n,0,1)
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=Cauchy(),control=glm.control(maxit=1000))
    }
    else if(charmatch(dist,"Student",F))
    {
      resp<-rt(n,arg)
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=Student(arg),control=glm.control(maxit=1000))
    }
    else if(charmatch(dist,"Gstudent",F))
    {
      resp<-rgstudent(n,arg[1],arg[2])
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=Gstudent(arg),control=glm.control(maxit=1000) )
    }
    else if(charmatch(dist,"LogisI",F))
    {
      stop(paste("not implemented yet"))
      resp<-rlogisI(n,0,1)
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=LogisI(),control=glm.control(maxit=1000))
    }
    else if(charmatch(dist,"LogisII",F))
    {
      resp<-rlogisII(n)
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=LogisII(),control=glm.control(maxit=1000))
    }
    else if(charmatch(dist,"Glogis",F))
    {
      stop(paste("not implement yet"))
      resp<-rglogis(n,arg[1],arg[2])
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=Glogis(arg),control=glm.control(maxit=1000))
      
    }
    else if(charmatch(dist,"Cnormal",F))
    {
      stop(paste("not implemented yet"))
      resp<-rcnormal(n,arg[1],arg[2])
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=Cnormal(arg),control=glm.control(maxit=1000))
    }
    else if(charmatch(dist,"Powerexp",F))
    {
      resp<-rpowerexp(n,arg)
      resp<-mu+sqrt(phi)*resp
      fit<-elliptical(resp~X+(-1),DerB=DerB,linear=object$linear,parmB=initial,family=Powerexp(arg),control=glm.control(maxit=1000))
    }
    ro <- fit$resid
    td <- ro/sqrt(fit$scalevariance)
    ####residual Cox and Snell
    Xd<-as.matrix(fit$Xmodel)
    Xdi<-solve(t(Xd)%*%Xd)
    H<-Xd%*%Xdi%*%t(Xd)
    H1<-(1/(fit$scalevariance*fit$scale))*H
    varr<-fit$scalevariance*fit$dispersion*(diag(1,n)-H1)
    varr<-diag(varr)
    ri<-fit$y-fit$fitted
    td<-ri/sqrt(varr)
    ############
    
    e[,i] <- sort(td)}
  #
  e1 <- numeric(n)
  e2 <- numeric(n)
  e3 <- numeric(n)
  e4 <- numeric(n)
  e5 <- numeric(n)
  e6 <- numeric(n)
  e7 <- numeric(n)
  #
  for(i in 1:n){
    eo <- sort(e[i,])
    e1[i] <- eo[ceiling(B*0.05)]
    e2[i] <- eo[ceiling(B*0.95)]}
  e3<-t(t(apply(e,2,mean)))
  e4<-t(t(apply(e,2,vari)))
  e5<-t(t(apply(e,2,skewn)))
  e6<-t(t(apply(e,2,kurt)))
  e7<-cbind(e3,e4,e5,e6)
  desc<-apply(e7,2,mean)
  med <- apply(e,1,mean)
  faixa <- range(tdf,e1,e2)
  screen(4)
  par(pty="s")
  qqnorm(tdf,xlab="Quantiles of N(0,1)",
         ylab="Standardized residual", ylim=faixa, pch=16)
  par(new=TRUE)
  #
  qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
  par(new=TRUE)
  qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
  par(new=TRUE)
  qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2)
  x<-list(mean = desc[1],var=desc[2], skewness = desc[3],kurtosis=desc[4])
  invisible(x)
}


".First.lib" <-
  function(lib, pkg)
  {
    library.dynam("elliptical", package = pkg, lib.loc = lib)  
    cat("\n")
    cat("------------------------------------------------\n")
    if(is.R()){
      cat(package.description("elliptical", lib = lib, field="Title"))
      cat("\n")
      ver <- package.description("elliptical", lib = lib, field="Version")
      cat(paste("elliptical version", ver,  "is now loaded\n"))
    }
    cat("------------------------------------------------\n")
    cat("\n")
    return(invisible(0))
  }
