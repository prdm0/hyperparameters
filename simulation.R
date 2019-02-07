rm(list = ls())
######################################## LIBRARIES ###########################################

library(mvtnorm) 
library(sm)
library(MASS)
library(iRegression)
library(quantreg)
library(robustbase)

library(foreach)
library(doParallel)

# Altere dir para o diretorio em que foi salvo 
######################################## FUNCTIONS ###########################################
# o codigo C++ de nome fast.cpp.
directory <- "/home/pedro/Downloads/hyperparameters/"
setwd(dir = directory)
source(file = paste(directory,"config_cpp.R", sep = ""))
config_cpp() #Compilando o codigo C++.

#COMBINA??O DE HIPERPARAMETROS
comb = function(x){
  n = length(x)
  s = NULL
  for (i in 1:n) s = c(s, rep(x[i],n))
  tcomb = cbind(s, rep(x,n))
  colnames(tcomb) = c('A','B')
  tcomb
}

#HIPERPARAMETROS
S1 = function(y1, y2, frac = 0.5){
  n = length(y1)
  m = floor(n*frac)
  idx1 = sample(1:n, m, replace = T)
  idx2 = sample(1:n, m, replace = T)
  tmp = (y1[idx1] - y2[idx2])^2
  mean(quantile(tmp[tmp != 0], probs = c(.9, .1)))
}

# S2 = function(y1, y2){
#     D = outer(y1, y2, '-')
#     D = D^2
#     D_no_zero = D[which(!D == 0)]
#     median(D_no_zero)
#  }
# 
# #S2 <- S2_cpp

S3 = function(y1, y2, n, p){
  sum((y1-y2)^2)/(n-p-1)
}

S4 = function(y1,y2){
  h.select(y1,y2,method="aicc")
}

COV = function(y1,yest1,y2,yest2){
  a1 = yest1 - yest2/2
  b1 = yest1 + yest2/2
  
  a2 = y1 - y2/2
  b2 = y1 + y2/2
  
  Xbar1 = sum(a1 + b1)/ (2 * length(y1))
  Xbar2 = sum(a2 + b2)/ (2 * length(y1))
  
  cov = sum( 2*(a1 - Xbar1)*(a2 - Xbar2) + (a1 - Xbar1)*(b2 - Xbar2) + (b1 - Xbar1)*(a2 - Xbar2) + 2*(b1 - Xbar1)*(b2 - Xbar2) ) / (6*length(y1))
  return(cov)
}

#KERNEL GAUSSIANO
gauss.kern = function(a, b, s){
  as.vector(exp(-( 1/s)*(a-b)^2)) #/2 
}

ETKRR = function(y, x, s, tol = 1e-10, maxit = 100, tolsolve = 1e-100) {
  invg = 0
  x = as.matrix(x)
  n = nrow(x)
  p = ncol(x)
  
  x = cbind(1, x)
  
  # Initialization
  
  txkx = t(x)%*%x
  if(det(txkx) == 0 || abs(det(txkx)) < tol) {invtxkx = ginv(txkx); invg = invg + 1} else {invtxkx = solve(txkx, tol = tolsolve)}
  
  
  betahat = invtxkx%*%t(x)%*%y #BETA OLS
  yhat = x%*%betahat #YEST OLS
  
  hparameter = switch(s, 
                      S1 = S1(y, yhat),
                      S2 = S2(y, yhat),
                      S3 = S3(y, yhat, n, p),
                      S4 = S4(y, yhat))
  
  K = gauss.kern(y, yhat, hparameter)
  S = sum(2 - 2*K)
  
  it = 1
  # Model Step
  
  repeat {
    it = it+1
    
    txkx = t(x)%*%diag(K)%*%x
    if(det(txkx) == 0 || abs(det(txkx)) < tol) {invtxkx = ginv(txkx); invg = invg + 1} else {invtxkx = solve(txkx, tol = tolsolve)}
    
    
    betahat = invtxkx%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, hparameter)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coefficients = as.vector(betahat), fitted.values = as.vector(yhat), criterion = S, weigth = K, iter = it, nginv = invg, hp = hparameter))
}

#ETKRR PRODUTO
ETKRRP =  function(y1, x1, y2, x2, s, t = NULL, tol = 1e-10, maxit = 100, tolsolve = 1e-100) {

  #1: CENTRO/MAX
  #2: AMPLITUDE/MIN
  #TOLSOLVE ? A TOLERANCIA DE CADA ELEMENTO DA MATRIZ, FOI COLOCADO ESSA TOLERANCIA POIS FOI OBSERVADO NUMEROS MUITO PEQUENOS COMO 1e-50
  invg = 0
  
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  
  n = nrow(x1)
  p = ncol(x1)
  
  x1 = cbind(1, x1) #CENTRO
  x2 = cbind(1, x2) #AMPLITUDE
  
  #BETAS CENTRO
  txkx1 = t(x1)%*%x1
  if(det(txkx1) == 0 || abs(det(txkx1)) < tol) {invtxkx1 = ginv(txkx1); invg = invg + 1} else {invtxkx1 = solve(txkx1, tol = tolsolve)}
  
  betahat1 = invtxkx1%*%t(x1)%*%y1 
  #Y1 ESTIMADO OLS
  yhat1 = x1%*%betahat1 
  
  #BETAS AMPLITUDE
  
  txkx2 = t(x2)%*%x2
  if(det(txkx2) == 0 || abs(det(txkx2)) < tol) {invtxkx2 = ginv(txkx2); invg = invg + 1} else {invtxkx2 = solve(txkx2, tol = tolsolve)}
  
  betahat2 = invtxkx2%*%t(x2)%*%y2 
  #Y2 ESTIMADO OLS
  yhat2 = x2%*%betahat2 
  
  
  
  #CASO 1, HIPERPARAMETROS IGUAIS
  
  if (!is.null(t)){
    
    if(!is.null(s)){
      hp1 = switch(s[1], S1 = S1(y1, yhat1), S2 = S2(y1, yhat1), S3 = S3(y1, yhat1, n, p), S4 = S4(y1, yhat1))
      hp2 = switch(s[1], S1 = S1(y2, yhat2), S2 = S2(y2, yhat2), S3 = S3(y2, yhat2, n, p), S4 = S4(y2, yhat2))}
    hp = switch(t, MAX = max(hp1, hp2), MIN = min(hp1, hp2), MEAN = (hp1+hp2)/2, COV = COV(y1, yhat1, y2, yhat2))
    hp1 = hp; hp2 = hp;
    
  } 
  #CASO 2, HIPERPARAMETROS DIFERENTES
  else {
    
    hp1 = switch(s[1], S1 = S1(y1, yhat1), S2 = S2(y1, yhat1), S3 = S3(y1, yhat1, n, p), S4 = S4(y1, yhat1))
    hp2 = switch(s[2], S1 = S1(y2, yhat2), S2 = S2(y2, yhat2), S3 = S3(y2, yhat2, n, p), S4 = S4(y2, yhat2))
    
  }
  
  K = gauss.kern(y1, yhat1, hp1) * gauss.kern(y2, yhat2, hp2)
  S = sum(2 - 2*K)
  
  it = 1
  
  repeat {
    it = it+1
    
    txkx1 = t(x1)%*%diag(K)%*%x1
    if(det(txkx1) == 0 || abs(det(txkx1)) < tol) {invtxkx1 = ginv(txkx1); invg = invg + 1} else {invtxkx1 = solve(txkx1, tol = tolsolve)}
    
    betahat1 = invtxkx1 %*% t(x1 )%*%diag(K)%*%y1
    yhat1 = x1%*%betahat1
    
    txkx2 = t(x2)%*%diag(K)%*%x2
    if(det(txkx2) == 0 || abs(det(txkx2)) < tol) {invtxkx2 = ginv(txkx2); invg = invg + 1} else {invtxkx2 = solve(txkx2, tol = tolsolve)}
    
    betahat2 = invtxkx2%*%t(x2)%*%diag(K)%*%y2
    yhat2 = x2%*%betahat2
    
    K = gauss.kern(y1, yhat1, hp1) * gauss.kern(y2, yhat2, hp2) 
    S = c(S, sum(2-2*K))
    
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coefficients = cbind(betahat1, betahat2), fitted.values = cbind(yhat1, yhat2), criterion = S, weigth = K, iter = it, nginv = invg, hp = c(hp1,hp2)))
}

func_metodos = function(amostra, percentual, scenario) {
    
    reg.elliptical = source(file = paste(directory,"ellipticalReg.R", sep = ""))
    
    tol = 1e-10 #TOLERANCIA DO MODELO
    maxit = 100 #MAXIMO DE ITERA??O DO MODELO
    betac = c(11,2)
    betar = c(1,1)
    
    cenario = scenario
    n = amostra
    nout = round(n * percentual)
    data.sum = matrix(0,51,18) #MATRIZ ONDE SER? SALVO A SOMA DAS ESTIMATIVAS 

        
        xc = runif(n, 5, 10)
        xr = runif(n, 1, 2)
        yc = betac[1] + betac[2] * xc + rnorm(n, 0, 0.5)
        yr = betar[1] + betar[2] * xr + rnorm(n, 0, 0.1)
        
        data = cbind(yc,xc,yr,xr);
        
        ################################## CENARIO #################################################
        
        #MATRIZ DE VARIANCIA E COVARIANCIA BIVARIADA
        sc = diag(0.5,2,2) #CENTER
        sr = diag(0.1,2,2) #RANGE
        
        #M?DIA DOS CENARIOS (Table 1 - LimaNetoDeCarvalho_2018)
        
        #CENTER m = mean; ms = max and standard deviation
        m.xc = mean(xc)
        ms.xc = max(xc) + 1.5 * sd(xc)
        m.yc = mean(yc)
        ms.yc = max(yc) + 1.5 * sd(yc)
        
        #RANGE m = mean; ms = max and standard deviation
        m.xr = mean(xr)
        ms.xr = max(xr) + 1.5 * sd(xr)
        m.yr = mean(yr)
        ms.yr = max(yr) + 1.5 * sd(yr)
        
        #VETORES DE M?DIA DOS CENARIOS 
        ucx = c(m.xc, ms.xc, ms.xc, 0, 0, 0, m.xc, ms.xc, ms.xc)
        ucy = c(ms.yc, m.yc, ms.yc, 0, 0, 0, ms.yc, m.yc, ms.yc)
        urx = c(0,0,0, m.xr, ms.xr,ms.xr, m.xr, ms.xr, ms.xr)
        ury = c(0,0,0, ms.yr, m.yr, ms.yr, ms.yr,m.yr,ms.yr )
        
        #M?DIAS DA NORMAL BIVARIADA
        uc = c(ucy[cenario], ucx[cenario]) #CENTER
        ur = c(ury[cenario], urx[cenario]) #RANGE
        
        #OUTLIERS CENTER E RANGE
        outc = rmvnorm(nout, mean = uc , sigma = sc); colnames(outc) = c('yc', 'xc')
        outr = rmvnorm(nout, mean = ur , sigma = sr); colnames(outr) = c('yr', 'xr')
        
        #OUTLIERS USANDO COLUNA DO BANCO CASO OUTLIER SEJA APENAS CENTER OU RANGE, E OU OUTLIERS USANDO COLUNAS OUTC E OUTR
        if (sum(uc) != 0 & sum(ur) == 0  ) out = cbind(outc, data[1:nout,c(3,4)])
        if (sum(uc) == 0 & sum(ur) != 0  ) out = cbind(data[1:nout,c(1,2)], outr)
        if (sum(uc) != 0 & sum(ur) != 0  ) out = cbind(outc, outr)
        
        #obs: OUT S?O OS OUTLIERS GERADOS
        
        #BANCO COM OUTLIERS
        cenario.out = rbind(data[-(1:nout),], out)
        
        y1 = cenario.out[,'yc']
        y2 = cenario.out[,'yr']
        x1 = cenario.out[,'xc']
        x2 = cenario.out[,'xr']
        
        id = 1 #LINHA QUE REPRESENTA CADA M?TODO
        
        data.cenario = matrix(0,51,18)
        
        colname.data =  c('B0c(mean)', 'B0c(bias)', 'B0c(mse)', 'B1c(mean)', 'B1c(bias)', 
                          'B1c(mse)','B0r(mean)', 'B0r(bias)', 'B0r(mse)', 'B1r(mean)', 
                          'B1r(bias)', 'B1r(mse)', 'Iter','Pnr','Prbc', 'Time', 'IterMax', 'Ginv')
        # B0c(mean) B0c(bias) B0c(mse) B1c(mean) B1c(bias) B1c(mse) B0r(mean) B0r(bias) B0r(mse) B1r(mean) B1r(bias) B1r(mse) Iter Pnr Prbc  Time IterMax Ginv
        colnames(data.cenario) = colname.data #NOMES DAS COLUNAS
        chp_names = NULL  #COMBINA??ES DE HIPERPARAMETROS
        
        ########################### ETKRR #######################################################
        
        S = c('S1', 'S2', 'S3', 'S4')
        chp = comb(S);
        
        for (i in 1:nrow(chp)) {
          
          #HIPERPARAMETROS
          hp = chp[i,]
          hp1 = hp[1] ; hp2 = hp[2];
          
          #TEMPO DO AJUSTE DO MODELO
          t0 = proc.time()[3]
          mod.k_c = ETKRR(y1, x1, hp1, tol = tol, maxit = maxit )
          mod.k_r = ETKRR(y2, x2, hp2, tol = tol, maxit = maxit )
          data.cenario[id, 'Time'] = proc.time()[3] - t0
          
          data.cenario[id, c('B0c(mean)','B1c(mean)')] = mod.k_c$coefficients #AVERAGE
          data.cenario[id, c('B0c(bias)','B1c(bias)')] = betac - mod.k_c$coefficients #BIAS
          data.cenario[id, c('B0c(mse)','B1c(mse)')] = (betac - mod.k_c$coefficients)^2 #MSE
          
          data.cenario[id, c('B0r(mean)','B1r(mean)')] = mod.k_r$coefficients
          data.cenario[id, c('B0r(bias)','B1r(bias)')] = betar - mod.k_r$coefficients
          data.cenario[id, c('B0r(mse)','B1r(mse)')] = (betar - mod.k_r$coefficients)^2
          
          data.cenario[id, 'Ginv'] = max(mod.k_c$nginv, mod.k_r$nginv)
          data.cenario[id, 'Iter'] = max(mod.k_c$iter, mod.k_r$iter)
          if (max(mod.k_c$iter, mod.k_r$iter) >= maxit) data.cenario[id, 'IterMax'] = 1 
          
          data.cenario[id,'Pnr'] = sum(mod.k_r$fitted.values<0)/length(y1) 
          data.cenario[id,'Prbc'] = sum(mod.k_r$fitted.values >= mod.k_c$fitted.values)/ length(y1) 
          
          #NOME DAS COMBINA??ES
          chp_name = paste("iETKRR ", "(", hp1, ",", hp2, ")", sep = ''  );
          chp_names = c(chp_names, chp_name)
          id = id + 1
        }
        
        ################################### ETKRR PRODUTO  CASO 1 #################################
        
        t = c('MAX', 'MIN', 'MEAN')
        S = c('S1', 'S2', 'S3', 'S4')
        
        for (i in 1:length(S)) {
          
          for (j in 1:length(t)) { #HIPERPARAMETROS (S1,S2,S3,S4)
            
            t0 = proc.time()[3]
            mod.k_p = ETKRRP(y1,x1,y2,x2, S[i], t[j], tol = tol, maxit = maxit) #SINGULAR S3 E MIN 
            data.cenario[id, 'Time'] = proc.time()[3] - t0
            
            data.cenario[id, c('B0c(mean)','B1c(mean)')] = mod.k_p$coefficients[,1] #AVERAGE
            data.cenario[id, c('B0c(bias)','B1c(bias)')] = betac - mod.k_p$coefficients[,1] #BIAS
            data.cenario[id, c('B0c(mse)','B1c(mse)')] = (betac - mod.k_p$coefficients[,1])^2 #MSE
            
            data.cenario[id, c('B0r(mean)','B1r(mean)')] = mod.k_p$coefficients[,2]
            data.cenario[id, c('B0r(bias)','B1r(bias)')] = betar - mod.k_p$coefficients[,2]
            data.cenario[id, c('B0r(mse)','B1r(mse)')] = (betar - mod.k_p$coefficients[,2])^2
            
            data.cenario[id, 'Ginv'] = mod.k_p$nginv
            data.cenario[id, 'Iter'] = mod.k_p$iter
            if (mod.k_p$iter >= maxit) data.cenario[id, 'IterMax'] = 1
            
            data.cenario[id,'Pnr'] = sum(mod.k_r$fitted.<0)/length(y1) 
            data.cenario[id,'Prbc'] = sum(mod.k_p$fitted[,2] >= mod.k_p$fitted[,1])/ length(y1) 
            
            
            chp_name = paste("iETKRR(P) ",t[j], "(", S[i],")", sep = ''  ); #ROW NAME
            chp_names = c(chp_names, chp_name) #ROWS NAMES
            
            id = id + 1
          }
        }
        
        ############################### ETKRR PRODUTO CASO 1 COVARIANCIA ###########################
        
        t0 = proc.time()[3]
        mod.k_p = ETKRRP(y1,x1,y2,x2, 'Y,y', 'COV', tol = tol, maxit = maxit) #SINGULAR S3 E MIN 
        data.cenario[id, 'Time'] = proc.time()[3] - t0
        
        data.cenario[id, c('B0c(mean)','B1c(mean)')] = mod.k_p$coefficients[,1] #AVERAGE
        data.cenario[id, c('B0c(bias)','B1c(bias)')] = betac - mod.k_p$coefficients[,1] #BIAS
        data.cenario[id, c('B0c(mse)','B1c(mse)')] = (betac - mod.k_p$coefficients[,1])^2 #MSE
        
        data.cenario[id, c('B0r(mean)','B1r(mean)')] = mod.k_p$coefficients[,2]
        data.cenario[id, c('B0r(bias)','B1r(bias)')] = betar - mod.k_p$coefficients[,2]
        data.cenario[id, c('B0r(mse)','B1r(mse)')] = (betar - mod.k_p$coefficients[,2])^2
        
        data.cenario[id, 'Ginv'] = mod.k_p$nginv
        data.cenario[id, 'Iter'] = mod.k_p$iter
        if (mod.k_p$iter >= maxit) data.cenario[id, 'IterMax'] = 1
        
        data.cenario[id,'Pnr'] = sum(mod.k_r$fitted.<0)/length(y1) 
        data.cenario[id,'Prbc'] = sum(mod.k_p$fitted[,2] >= mod.k_p$fitted[,1])/ length(y1) 
        
        
        chp_name = paste("iETKRR(P) ",'COV(Y,y)', sep = ''  ); #ROW NAME
        chp_names = c(chp_names, chp_name) #ROWS NAMES
        
        id = id + 1
        
        ################################## ETKRR PRODUTO CASO 2 ####################################
        
        S = c('S1','S2','S3','S4')
        chp = comb(S);
        
        for (i in 1:nrow(chp)) {
          
          hp = chp[i,]
          
          t0 = proc.time()[3]
          mod.k_p = ETKRRP(y1,x1,y2,x2, hp, tol = tol, maxit = maxit) 
          data.cenario[id, 'Time'] = proc.time()[3] - t0
          
          data.cenario[id, c('B0c(mean)','B1c(mean)')] = mod.k_p$coefficients[,1]
          data.cenario[id, c('B0c(bias)','B1c(bias)')] = betac - mod.k_p$coefficients[,1]
          data.cenario[id, c('B0c(mse)','B1c(mse)')] = (betac - mod.k_p$coefficients[,1])^2
          
          data.cenario[id, c('B0r(mean)','B1r(mean)')] = mod.k_p$coefficients[,2]
          data.cenario[id, c('B0r(bias)','B1r(bias)')] = betar - mod.k_p$coefficients[,2]
          data.cenario[id, c('B0r(mse)','B1r(mse)')] = (betar - mod.k_p$coefficients[,2])^2
          
          data.cenario[id, 'Ginv'] = mod.k_p$nginv
          data.cenario[id, 'Iter'] = mod.k_p$iter
          if (mod.k_p$iter >= maxit) data.cenario[id, 'IterMax'] = 1
          
          data.cenario[id,'Pnr'] = sum(mod.k_r$fitted.<0)/length(y1) 
          data.cenario[id,'Prbc'] = sum(mod.k_p$fitted[,2] >= mod.k_p$fitted[,1])/ length(y1) 
          
          
          chp_name = paste("iETKRR(P) ", "(", hp[1], ",", hp[2], ")", sep = ''  )
          chp_names = c(chp_names, chp_name)
          
          id = id + 1
        }
        
        ########################## IRR - Fagundes et al (2013) ####################################
        
        t0 = proc.time()[3]
        mod.irr_c = lmrob(y1 ~ x1)
        mod.irr_r = lmrob(y2 ~ x2)
        data.cenario[id, 'Time'] = proc.time()[3] - t0
        
        data.cenario[id, c('B0c(mean)','B1c(mean)')] = mod.irr_c$coefficients
        data.cenario[id, c('B0c(bias)','B1c(bias)')] = betac - mod.irr_c$coefficients
        data.cenario[id, c('B0c(mse)','B1c(mse)')] = (betac - mod.irr_c$coefficients)^2
        
        data.cenario[id, c('B0r(mean)','B1r(mean)')] = mod.irr_r$coefficients
        data.cenario[id, c('B0r(bias)','B1r(bias)')] = betar - mod.irr_r$coefficients
        data.cenario[id, c('B0r(mse)','B1r(mse)')] = (betar - mod.irr_r$coefficients)^2
        
        data.cenario[id, 'Ginv'] = 0
        data.cenario[id, 'Iter'] = max( mod.irr_c$iter, mod.irr_r$iter)
        if ( max( mod.irr_c$iter, mod.irr_r$iter) >= maxit) data.cenario[id, 'IterMax'] = 1
        
        data.cenario[id,'Pnr'] = sum(mod.irr_r$fitted<0)/length(y1) 
        data.cenario[id,'Prbc'] = sum(mod.irr_r$fitted >= mod.irr_c$fitted)/ length(y1) 
        
        
        chp_name = paste("IRR" , sep = '')
        chp_names = c(chp_names, chp_name)
        
        id = id + 1
        
        ################################## QIR ################################################
        
        t0 = proc.time()[3]
        mod.qir_c = rq(y1 ~ x1, 0.5)
        mod.qir_r = rq(y2 ~ x2, 0.5)
        data.cenario[id, 'Time'] = proc.time()[3] - t0
        
        data.cenario[id, c('B0c(mean)','B1c(mean)')] = mod.qir_c$coefficients
        data.cenario[id, c('B0c(bias)','B1c(bias)')] = betac - mod.qir_c$coefficients
        data.cenario[id, c('B0c(mse)','B1c(mse)')] = (betac - mod.qir_c$coefficients)^2
        
        data.cenario[id, c('B0r(mean)','B1r(mean)')] = mod.qir_r$coefficients
        data.cenario[id, c('B0r(bias)','B1r(bias)')] = betar - mod.qir_r$coefficients
        data.cenario[id, c('B0r(mse)','B1r(mse)')] = (betar - mod.qir_r$coefficients)^2
        
        data.cenario[id, 'Ginv'] = 0
        data.cenario[id, 'Iter'] = 0
        #if ( max( mod.qir_c$iter, mod.qir_r$iter) >= maxit) 
        data.cenario[id, 'IterMax'] = 0
        
        data.cenario[id,'Pnr'] = sum(mod.qir_r$fitted<0)/length(y1) 
        data.cenario[id,'Prbc'] = sum(mod.qir_r$fitted >= mod.qir_c$fitted)/ length(y1) 
        
        chp_name = paste("QIR" , sep = '')
        chp_names = c(chp_names, chp_name)
        
        id = id + 1
        
        ##################################### SSLR - t-Student ########################################
        
        gl = c(2,6,8,10)
        
        for (i in 1:length(gl)) {
          
          t0 = proc.time()[3]
          mod.sslr_c = elliptical(formula= y1 ~ x1,family=Student(gl[i]))
          mod.sslr_r = elliptical(formula= y2 ~ x2,family=Student(gl[i]))
          data.cenario[id, 'Time'] = proc.time()[3] - t0
          
          data.cenario[id, c('B0c(mean)','B1c(mean)')] = mod.sslr_c$coefficients #AVERAGE
          data.cenario[id, c('B0c(bias)','B1c(bias)')] = betac - mod.sslr_c$coefficients #BIAS
          data.cenario[id, c('B0c(mse)','B1c(mse)')] = (betac - mod.sslr_c$coefficients)^2 #MSE
          
          data.cenario[id, c('B0r(mean)','B1r(mean)')] = mod.sslr_r$coefficients
          data.cenario[id, c('B0r(bias)','B1r(bias)')] = betar - mod.sslr_r$coefficients
          data.cenario[id, c('B0r(mse)','B1r(mse)')] = (betar - mod.sslr_r$coefficients)^2
          
          data.cenario[id, 'Ginv'] = 0
          data.cenario[id, 'Iter'] = max(mod.sslr_c$iter, mod.sslr_r$iter)
          if (max(mod.sslr_c$iter, mod.sslr_r$iter) >= maxit) data.cenario[id, 'IterMax'] = 1 
          
          data.cenario[id,'Pnr'] = sum(mod.sslr_r$fitted.values<0)/length(y1) 
          data.cenario[id,'Prbc'] = sum(mod.sslr_r$fitted.values >= mod.sslr_c$fitted.values)/ length(y1) 
          
          #NOME DAS COMBINA??ES
          chp_name = paste("SSLR ", "(", gl[i], ",", gl[i], ")", sep = ''  );
          chp_names = c(chp_names, chp_name)
          id = id + 1
          
        }
        #rownames(data.cenario) = chp_names
        data.sum = data.sum + data.cenario
        rownames(data.sum) = chp_names
        
        return(list("data" = data.sum))}



        ###################################### FIM DA FUN??O #################################################

func_metodos(1000,0.05,1)
packpages = c('mvtnorm', 'sm', 'MASS', 'iRegression', 'quantreg', 'robustbase') 

profvis::profvis(func_metodos(1000,0.05,1)) 
   


    
          
################################## PEQUENA COMPARA??O ENTRE FOR NORMAL E FOR PARALELO ################################

       
        
#50 - Tamanho da amostra
#0.05 - Percentual de Outliers
#1 - Cenario 
        
#Normal
set.seed(0)
TempoA = system.time({for (i in 1:7) {func_metodos(1000,0.05,1)}})[3]
#S2 <- S2_cpp
TempoB = system.time({for (i in 1:7) {func_metodos(1000,0.05,1)}})[3]
        
#Paralelo

#Checa quantos n?cleos existem
#ncl<-getDoParWorkers()
#Checa quantos n?cleos existem
# ncl = detectCores() - 1
# 
# #Registra os n?cleos a serem utilizados
# cl <- makeCluster(ncl, type = 'PSOCK')
# registerDoParallel(cl)
# 
# #checa quantos nucleos est?o sendo usados
# getDoParWorkers()
# 
# TempoB = system.time({foreach(i = 1:7, .packages=packpages) %dopar% {func_metodos(1000,0.05,1)}})[3]
# 
# stopCluster(cl)  
# 
# comp_func_metodos <- cmpfun(func_metodos)
# 
# TempoC = system.time({foreach(i = 1:100, .packages=packpages) %dopar% {comp_func_metodos(200,0.05,1)}})[3]


#    TempoA   TempoB(4)  TempoC(3) 200 50
#   136         62.41     57.34

#    TempoA   TempoB(4)  TempoC(3) TempoD(2) 200 100
#    261.13    101.42     106.70    127.09