#setting g1/g2/g3 high dim
#require(glmnet)
ptm=proc.time()
set.seed(123)
N=2000;p=100;pi_N=0.05;s=3;tmax=500;K=5;M=floor(N/K)
#g1
#c_signal=0.1 #ORgood_PSbad1 (->ORbad_PSbad)
#c_signal=0.3 #ORgood_PSbad2 (->ORbad_PSbad)
#c_signal=0.5 #ORgood_PSbad3 (->ORbad_PSbad)
#c_signal=1 #ORgood_PSbad4 (->ORbad_PSbad)
#beta=c(2,1*c(1,1,1),rep(0,p-3))
#gamma=0.3*c(0,1,1,1,rep(0,p-3)) 
#delta=0.3*c(0,1,1,1,rep(0,p-3))

#g3
#c_signal=1 #ORgood_PSgood1234
#beta=c(2,0*c(1,1,1),rep(0,p-3))
#gamma=1*c(0,1,1,1,rep(0,p-3))
#delta=0.05*c(0,1,1,1,rep(0,p-3)) #ORgood_PSgood1 (->ORbad_PSbad)
#delta=0.1*c(0,1,1,1,rep(0,p-3)) #ORgood_PSgood2 
#delta=0.2*c(0,1,1,1,rep(0,p-3)) #ORgood_PSgood3
#delta=0.3*c(0,1,1,1,rep(0,p-3)) #ORgood_PSgood4

#g2
beta=c(2,0.3*c(1,1,1),rep(0,p-3)) #ORbad_PSgood1234
gamma=0.5*c(0,1,1,1,rep(0,p-3)) 
#delta=0.05*c(0,1,1,1,rep(0,p-3)) #ORbad_PSgood1 (->ORbad_PSbad)
#delta=0.1*c(0,1,1,1,rep(0,p-3)) #ORbad_PSgood2 
#delta=0.2*c(0,1,1,1,rep(0,p-3)) #ORbad_PSgood3
delta=0.3*c(0,1,1,1,rep(0,p-3)) #ORbad_PSgood4
delta2=0.3*c(0,1,1,1,rep(0,p-3))

#theta=beta[1]+c_signal*sum(delta[-1]) #ORgood_PSbad and ORgood_PSgood
theta=beta[1]+sum(delta2[-1]) #ORbad_PSgood

func_link=function(x){exp(x)/(1+exp(x))}
#find the gamma(1) corresponds to gamma(-1) and pi_N
intercept_lower=-10;intercept_upper=10
N_sim=100000
Xs=matrix(rnorm(N_sim*s),nrow=N_sim)
Xgamma=Xs%*%gamma[2:(s+1)]+Xs^2%*%delta[2:(s+1)]+0*Xs[,1]*Xs[,2]
for (i in 1:4){
  intercept_range=seq(from=intercept_lower,to=intercept_upper,by=10^(-i))
  values=apply(t(intercept_range),2,function(x){mean(func_link(Xgamma+x))})
  intercept_upper=intercept_range[min(which(values-pi_N>0))]
  intercept_lower=intercept_upper-10^(-i)
  if (i==4){
    intercept=intercept_range[which.min(abs(values-pi_N))]
  }
}
gamma[1]=intercept
scale_para=sd(func_link(Xgamma+intercept))*5
est_mean=function(x){
  mean(pred_m[x[1],]+R*(Y-pred_m[x[1],])/pred_pi[x[2],])
}
est_var_opt=function(x){
  var(pred_m[x[1],]+R*(Y-pred_m[x[1],])/pred_pi[x[2],])
}
num_mest=2;num_piest=3;num_totest=(num_mest-1)*(num_piest-1)+2
pred_m=matrix(0,nrow=num_mest,ncol=N);pred_pi=matrix(0,nrow=num_piest,ncol=N)
pred_theta=matrix(0,nrow=tmax,ncol=num_totest+3)
var_opt=matrix(0,nrow=tmax,ncol=num_totest+3);var_mod=matrix(0,nrow=tmax,ncol=num_totest+3)
CI_lower_opt=matrix(0,nrow=tmax,ncol=num_totest+3);CI_upper_opt=matrix(0,nrow=tmax,ncol=num_totest+3)
mse_m=matrix(0,nrow=tmax,ncol=num_mest);mse_pi=matrix(0,nrow=tmax,ncol=num_piest)
for (t in 1:tmax){
  X=matrix(rnorm(N*p),nrow=N)
  Xv=cbind(rep(1,N),X)
  eps=rnorm(N)
  pi=func_link(Xv%*%gamma+X^2%*%delta[-1])
  #m=Xv%*%beta+c_signal*(X%*%gamma[-1]+X^2%*%delta[-1]) #ORgood_PSbad and ORgood_PSgood
  m=Xv%*%beta+X^2%*%delta2[-1] #ORbad_PSgood
  Y=m+eps
  R=rbinom(N,1,pi)
  pred_m[1,]=m #oracle
  pred_pi[1,]=pi #oracle
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    X_nk=X[-index,];Y_nk=Y[-index];R_nk=R[-index]
    fit_LS=cv.glmnet(x=X_nk[R_nk==1,],y=Y_nk[R_nk==1],parallel=TRUE)
    pred_m[2,index]=predict(fit_LS,newx=X[index,],s="lambda.min",nfolds=5) #Lasso
    #X_poly=cbind(X,X^2)
    #fit_poly=cv.glmnet(x=(X_poly[-index,])[R_nk==1,],y=Y_nk[R_nk==1],parallel=TRUE)
    #pred_m[3,index]=predict(fit_poly,newx=X_poly[index,],s="lambda.min",nfolds=5) #poly
    #propensity score models
    pred_pi[2,index]=rep(mean(R_nk),M) #MCAR
    fit_Log=cv.glmnet(x=X_nk,y=R_nk,family="binomial",parallel=TRUE)
    pred_pi[3,index]=predict(fit_Log,newx=X[index,],type="response",s="lambda.min",nfolds=5) #offset
    #fit_Log_poly=cv.glmnet(x=X_poly[-index,],y=R_nk,family="binomial",parallel=TRUE)
    #pred_pi[4,index]=predict(fit_Log_poly,newx=X_poly[index,],type="response",s="lambda.min",nfolds=5) #offset-poly
  }
  fit_Log=cv.glmnet(x=X,y=R,family="binomial",parallel=TRUE)
  pred_pi_IPW=predict(fit_Log,newx=X,type="response",s="lambda.min",nfolds=5) #offset
  pred_theta_IPW_POP=sum(R*Y/pred_pi_IPW)/sum(R/pred_pi_IPW) #IPW-POP
  pred_theta_IPW_NR=mean(R*Y)+mean(1-R)*sum(R*Y*(1-pred_pi_IPW)/pred_pi_IPW)/sum(R*(1-pred_pi_IPW)/pred_pi_IPW) #IPW-NR
  fit_LS=cv.glmnet(x=X[R==1,],y=Y[R==1],parallel=TRUE)
  pred_m_OLS=predict(fit_LS,newx=X,s="lambda.min",nfolds=5) #Lasso
  pred_theta_OLS=mean(pred_m_OLS)
  #check mse of the outcome model
  mse_m[t,]=c(mean((pred_m[1,]-pred_m[-1,])^2),mean((pred_m[1,]-pred_m_OLS)^2))
  #check mse of the propensity score model
  mse_pi[t,]=c(colMeans((1-pred_pi[1,]/t(pred_pi[-1,]))^2),mean((1-pred_pi[1,]/t(pred_pi_IPW))^2))
  #mean estimator
  grid=expand.grid(est_m=1:num_mest,est_pi=1:num_piest)
  w=which(grid$est_m==1|grid$est_pi==1)[-1]
  grid=grid[-w,]
  pred_theta[t,1]=sum(R*Y)/sum(R) #average over labeled samples
  pred_theta[t,2:num_totest]=apply(grid,1,est_mean)
  pred_theta[t,num_totest+1:3]=c(pred_theta_IPW_POP,pred_theta_IPW_NR,pred_theta_OLS)
  #asymptotic variance estimator
  var_opt[t,2:num_totest]=apply(grid,1,est_var_opt) #optimal
  var_opt[t,1]=var(Y[which(R==1)])*N/sum(R)
  var_opt[t,num_totest+1:3]=c(var(R*Y/pred_pi_IPW)/mean(R/pred_pi_IPW),
                              var(R*Y+mean(1-R)/mean(R*(1-pred_pi_IPW)/pred_pi_IPW)*R*Y*(1-pred_pi_IPW)/pred_pi_IPW),
                              var(pred_m_OLS))
  #CI
  CI_lower_opt[t,]=pred_theta[t,]-sqrt(var_opt[t,]/N)*qnorm(1-0.05/2)
  CI_upper_opt[t,]=pred_theta[t,]+sqrt(var_opt[t,]/N)*qnorm(1-0.05/2)
}
bias=colMeans(pred_theta-theta)
RMSE=colMeans((pred_theta-theta)^2)^0.5
AL_opt=colMeans(CI_upper_opt-CI_lower_opt)
AC_opt=colMeans((theta<=CI_upper_opt)&(theta>=CI_lower_opt))
ESD=apply(pred_theta,2,sd)
ASD_opt=colMeans(var_opt^0.5)/sqrt(N)
table_result=function(x){
  paste(format(round(bias[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(RMSE[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(AL_opt[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(AC_opt[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(ESD[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(ASD_opt[x],digits=3),nsmall=3,scientific=FALSE),sep="")
}
result=apply(t(1:(num_totest+3)),2,table_result)
N;pi_N;p
result
colMeans(mse_m)
colMeans(mse_pi)
proc.time()-ptm

