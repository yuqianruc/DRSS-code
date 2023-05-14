#setting e high dim
require(glmnet);require(itertools);require(foreach);require(doParallel)
ptm=proc.time()
#numCores=detectCores()
registerDoParallel(10)
N=10000;p=500;pi_N=0.01;tmax=500;K=5;M=floor(N/K)
beta=c(-0.5,1,1,1,rep(0,p-3))
alpha=c(0,1,1,1,rep(0,p-3))
theta=beta[1]+sum(alpha^2)
num_mest=3;num_piest=4;num_totest=num_mest*num_piest+2
grid=expand.grid(est_m=1:num_mest,est_pi=1:num_piest)
w=which(grid$est_m==1)[-(1:2)]
grid_nw=grid[-w,]
grid=rbind(grid_nw[1:(which(grid$est_m==1)[2]),],grid[w,],grid_nw[-(1:(which(grid$est_m==1)[2])),])
results_par <- foreach (t=1:tmax, .combine=rbind, .packages=c("glmnet")) %dopar% {
  pred_m=matrix(0,nrow=num_mest,ncol=N);pred_pi=matrix(0,nrow=num_piest,ncol=N)
  X=matrix(rnorm(N*p),nrow=N)
  eps=rnorm(N)
  Y=beta[1]+X%*%beta[-1]+alpha[1]+X^2%*%alpha[-1]+eps
  p_delta=exp(X[,1])/(1+exp(X[,1]))
  delta=rbinom(N,1,p_delta)
  pi_1=0.5*pi_N;pi_0=1.5*pi_N
  pi=pi_1*p_delta+pi_0*(1-p_delta)
  pi_delta=pi_1*delta+pi_0*(1-delta)
  R=rbinom(N,1,pi_delta)
  pred_m[1,]=Y-eps #oracle
  pred_pi[1,]=pi #oracle
  X=scale(X) #standardize the covariates
  pred_p_delta=rep(0,N)
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    Y_nk=Y[-index];R_nk=R[-index];delta_nk=delta[-index]
    Y_mtrain=Y_nk[which(R_nk==1)]
    #propensity score models
    pred_pi[2,index]=rep(mean(R[-index]),M) #MCAR
    lambda.range=sqrt(log(p)/((N-M)*mean(R[-index])))*10^seq(from=-3,to=0,length.out=100)
    fit_Log=cv.glmnet(x=X[-index,],y=R[-index],family="binomial",offset=rep(log(mean(R_nk)),length(R_nk)),lambda=lambda.range,nfolds=5)
    pred_pi[3,index]=predict(fit_Log,newx=X[index,],type="response",s="lambda.min",newoffset=rep(log(mean(R_nk)),length(index))) #offset
    fit_delta=cv.glmnet(x=X[-index,],y=delta_nk,family="binomial",nfolds=5)
    pred_pi_1=sum(delta_nk*R_nk)/sum(delta_nk);pred_pi_0=sum((1-delta_nk)*R_nk)/sum(1-delta_nk)
    pred_p_delta[index]=predict(fit_delta,newx=X[index,],type="response",s="lambda.min")
    pred_pi[4,index]=pred_pi_1*pred_p_delta[index]+pred_pi_0*(1-pred_p_delta[index]) #stratified
    #outcome models
    fit_LS=cv.glmnet(x=(X[-index,])[which(R_nk==1),],y=Y_mtrain,family="gaussian",alpha=1,nfolds=5)
    pred_m[2,index]=predict(fit_LS,newx=X[index,],s="lambda.min") #Lasso
    fit_poly2=cv.glmnet(x=cbind((X[-index,])[which(R_nk==1),],(X[-index,])[which(R_nk==1),]^2),y=Y_mtrain,family="gaussian",alpha=1,nfolds=5)
    pred_m[3,index]=predict(fit_poly2,newx=cbind(X[index,],X[index,]^2),s="lambda.min") #polynomial(degree 2, without interactions terms)
  }
  #check mse of the outcome model
  mse_m_t=colMeans((pred_m[1,]-t(pred_m[-1,]))^2)
  #check mse of the propensity score model
  mse_pi_t=colMeans((1-pred_pi[1,]/t(pred_pi[-1,]))^2)
  #mean estimator
  pred_theta_t=c();var_opt_t=c()
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    pred_theta_t=c(pred_theta_t,mean(pred_m[x[1],]+R*(Y-pred_m[x[1],])/pred_pi[x[2],]))
  }
  pred_theta_t=c(sum(R*Y)/sum(R),mean(Y),pred_theta_t)
  #asymptotic variance estimator
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    psi_mu=pred_m[x[1],]+R*(Y-pred_m[x[1],])/pred_pi[x[2],]-pred_theta_t[x[1]+num_mest*(x[2]-1)]
    var_opt_t=c(var_opt_t,mean(psi_mu^2))
  }
  var_opt_t=c(var(Y[which(R==1)])*N/sum(R),var(Y),var_opt_t)
  c(pred_theta=pred_theta_t,var_opt=var_opt_t,mse_m=mse_m_t,mse_pi=mse_pi_t)
}
pred_theta=results_par[,1:num_totest]
var_opt=results_par[,(num_totest+1):(2*num_totest)];
mse_m=results_par[,(2*num_totest+1):(2*num_totest+num_mest-1)]
mse_pi=results_par[,(2*num_totest+num_mest):(2*num_totest+num_mest+num_piest-2)]
mse_m_mean=colMeans(mse_m)
mse_m_sd=apply(mse_m,2,sd)
mse_pi_mean=colMeans(mse_pi)
mse_pi_sd=apply(mse_pi,2,sd)
bias=colMeans(pred_theta-theta)
RMSE=colMeans((pred_theta-theta)^2)^0.5
AL_opt=colMeans(2*sqrt(var_opt/N)*qnorm(1-0.05/2))
AC_opt=colMeans(abs(pred_theta-theta)<=sqrt(var_opt/N)*qnorm(1-0.05/2))
ESD=apply(pred_theta,2,sd)
ASD_opt=colMeans(var_opt^0.5)/sqrt(N)
table_result=function(x){
  paste(format(round(bias[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(RMSE[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(AL_opt[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(AC_opt[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(ESD[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(ASD_opt[x],digits=3),nsmall=3,scientific=FALSE),sep="")
}
result=apply(t(1:num_totest),2,table_result)
result=c(result,round(mse_m_mean,digits=3),round(mse_m_sd,digits=3),
         round(mse_pi_mean,digits=3),round(mse_pi_sd,digits=3))
result
proc.time()-ptm