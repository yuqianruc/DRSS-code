#setting c/d/f low dim (adj CI)
require(glmnet);require(itertools);require(foreach);require(doParallel)
ptm=proc.time()
set.seed(123)
#numCores=detectCores()
registerDoParallel(numCores)
N=10000;p=10;pi_N=0.01;tmax=500;K=5;M=floor(N/K)
beta=c(-0.5,1,1,1,rep(0,p-3))
alpha=rep(0,p+1) #O1
#alpha=c(0,1,1,1,rep(0,p-3)) #O2/O3
zeta=rep(0,p+1) #O1/O2
#zeta=c(0,0.2,0.2,0.2,rep(0,p-3)) #O3
#gamma=rep(0,p+1) #P1
gamma=c(0,1,rep(0,p-1)) #P2
#find the gamma(1) corresponds to gamma(-1) and pi_N
var_sig=sum(gamma[-1]^2)
gamma1_lower=-3;gamma1_upper=3
for (i in 1:4){
  gamma1_range=seq(from=gamma1_lower,to=gamma1_upper,by=10^(-i))
  i_max=100000
  i_range=1:(i_max-1)
  values=apply(t(gamma1_range),2,function(x){sum(1/(1+exp(-qnorm(i_range/i_max,mean=x+log(pi_N),sd=var_sig^0.5))))/(i_max-1)})
  gamma1_upper=gamma1_range[min(which(values-pi_N>0))]
  gamma1_lower=gamma1_upper-10^(-i)
  if (i==4){
    gamma1=gamma1_range[which.min(abs(values-pi_N))]
  }
}
gamma[1]=gamma1
theta=beta[1]+sum(alpha)
num_mest=4;num_piest=3;num_totest=num_mest*num_piest+2
grid=expand.grid(est_m=1:num_mest,est_pi=1:num_piest)
w=which(grid$est_m==1)[-(1:2)]
grid_nw=grid[-w,]
grid=rbind(grid_nw[1:(which(grid$est_m==1)[2]),],grid[w,],grid_nw[-(1:(which(grid$est_m==1)[2])),])
results_par <- foreach (t=1:tmax, .combine=rbind, .packages=c("grf","glmnet","KRLS")) %dopar% {
  pred_m=matrix(0,nrow=num_mest,ncol=N);pred_pi=matrix(0,nrow=num_piest,ncol=N)
  X=matrix(rnorm(N*p),nrow=N)
  eps=rnorm(N)
  Y=beta[1]+X%*%beta[-1]+X^2%*%alpha[-1]+X^3%*%zeta[-1]+eps
  pi=pi_N*exp(gamma[1]+X%*%gamma[-1])/(1+pi_N*exp(gamma[1]+X%*%gamma[-1]))
  R=rbinom(N,1,pi)
  pred_m[1,]=Y-eps #oracle
  pred_pi[1,]=pi #oracle
  X=scale(X) #standardize the covariates
  opt_param_RF=rep(0,6*K)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    Y_nk=Y[-index];R_nk=R[-index]
    Y_mtrain=Y_nk[which(R_nk==1)]
    #propensity score models
    pred_pi[2,index]=rep(mean(R[-index]),M) #MCAR
    fit_Log=glm(R_nk~.-pi_N+offset(log(pi_N)),data=data.frame(x=X[-index,],pi_N=rep(mean(R_nk),N-M)),family=binomial)
    pred_pi[3,index]=predict(fit_Log,newdata=data.frame(x=X[index,],pi_N=rep(mean(R_nk),M)),type="response") #offset
    #outcome models
    fit_LS=lm(Y_mtrain~.,data=data.frame(x=(X[-index,])[which(R_nk==1),]))
    pred_m[2,index]=predict(fit_LS,newdata=data.frame(x=X[index,])) #LS
    fit_poly2=lm(Y_mtrain~.,data=data.frame(x=cbind((X[-index,])[which(R_nk==1),],(X[-index,])[which(R_nk==1),]^2)))
    pred_m[3,index]=predict(fit_poly2,newdata=data.frame(x=cbind(X[index,],X[index,]^2))) #polynomial(degree 3, without interactions terms)
    fit_poly3=lm(Y_mtrain~.,data=data.frame(x=cbind((X[-index,])[which(R_nk==1),],(X[-index,])[which(R_nk==1),]^2,(X[-index,])[which(R_nk==1),]^3)))
    pred_m[4,index]=predict(fit_poly3,newdata=data.frame(x=cbind(X[index,],X[index,]^2,X[index,]^3))) #polynomial(degree 3, without interactions terms)
  }
  #check mse of the outcome model
  mse_m_t=colMeans((pred_m[1,]-t(pred_m[-1,]))^2)
  #check mse of the propensity score model
  mse_pi_t=colMeans((1-pred_pi[1,]/t(pred_pi[-1,]))^2)
  #mean estimator
  pred_theta_t=c();var_opt_t=c();var_mod_t=c()
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    pred_theta_t=c(pred_theta_t,mean(pred_m[x[1],]+R*(Y-pred_m[x[1],])/pred_pi[x[2],]))
  }
  pred_theta_t=c(sum(R*Y)/sum(R),mean(Y),pred_theta_t)
  #asymptotic variance estimator
  J3=cbind(t((pred_pi[3,]*(1-pred_pi[3,]))%*%cbind(rep(1,N),X)/N),mapply(function(y){(pred_pi[3,]*(1-pred_pi[3,])*X[,y])%*%cbind(rep(1,N),X)/N},1:p))
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    psi_mu=pred_m[x[1],]+R*(Y-pred_m[x[1],])/pred_pi[x[2],]-pred_theta_t[x[1]+num_mest*(x[2]-1)]
    var_opt_t=c(var_opt_t,mean(psi_mu^2))
    if (x[2]==1){
      var_mod_t=c(var_mod_t,mean(psi_mu^2)) #oracle
    } else if (x[2]==2){
      delta_mu=mean(pred_m[x[1],]-Y)
      var_mod_t=c(var_mod_t,mean((psi_mu+(R/mean(R)-1)*delta_mu)^2)) #MCAR
    } else if (x[2]==3){
      IF_pi=as.numeric(t((1-pred_pi[x[2],])*(pred_m[x[1],]-Y)/N)%*%cbind(rep(1,N),X)%*%solve(J3)%*%t(cbind(rep(1,N),X))*(R-pred_pi[x[2],]))
      var_mod_t=c(var_mod_t,mean((psi_mu+IF_pi)^2)) #offset
    }
  }
  var_opt_t=c(var(Y[which(R==1)])*N/sum(R),var(Y),var_opt_t)
  var_mod_t=c(var(Y[which(R==1)])*N/sum(R),var(Y),var_mod_t)
  c(pred_theta=pred_theta_t,var_opt=var_opt_t,var_mod=var_mod_t,mse_m=mse_m_t,mse_pi=mse_pi_t,opt_param_RF=opt_param_RF)
}
pred_theta=results_par[,1:num_totest]
var_opt=results_par[,(num_totest+1):(2*num_totest)];var_mod=results_par[,(2*num_totest+1):(3*num_totest)]
mse_m=results_par[,(3*num_totest+1):(3*num_totest+num_mest-1)]
mse_pi=results_par[,(3*num_totest+num_mest):(3*num_totest+num_mest+num_piest-2)]
mse_m_mean=colMeans(mse_m)
mse_m_sd=apply(mse_m,2,sd)
mse_pi_mean=colMeans(mse_pi)
mse_pi_sd=apply(mse_pi,2,sd)
bias=colMeans(pred_theta-theta)
RMSE=colMeans((pred_theta-theta)^2)^0.5
AL_opt=colMeans(2*sqrt(var_opt/N)*qnorm(1-0.05/2))
AL_mod=colMeans(2*sqrt(var_mod/N)*qnorm(1-0.05/2))
AC_opt=colMeans(abs(pred_theta-theta)<=sqrt(var_opt/N)*qnorm(1-0.05/2))
AC_mod=colMeans(abs(pred_theta-theta)<=sqrt(var_mod/N)*qnorm(1-0.05/2))
ESD=apply(pred_theta,2,sd)
ASD_opt=colMeans(var_opt^0.5)/sqrt(N)
ASD_mod=colMeans(var_mod^0.5)/sqrt(N)
table_result=function(x){
  paste(format(round(bias[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(RMSE[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(AL_opt[x],digits=3),nsmall=3,scientific=FALSE),"(",format(round(AL_mod[x],digits=3),nsmall=3,scientific=FALSE),")&",
        format(round(AC_opt[x],digits=3),nsmall=3,scientific=FALSE),"(",format(round(AC_mod[x],digits=3),nsmall=3,scientific=FALSE),")&",
        format(round(ESD[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(ASD_opt[x],digits=3),nsmall=3,scientific=FALSE),"(",
        format(round(ASD_mod[x],digits=3),nsmall=3,scientific=FALSE),")",sep="")
}
result=apply(t(1:num_totest),2,table_result)
result=c(result,round(mse_m_mean,digits=3),round(mse_m_sd,digits=3),
         round(mse_pi_mean,digits=3),round(mse_pi_sd,digits=3))
result
proc.time()-ptm

