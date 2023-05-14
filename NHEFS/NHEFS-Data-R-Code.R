#NHEFS
ptm=proc.time()
require(readxl);require(glmnet);require(grf);require(KRLS)
NHEFS=read_excel("NHEFS.xls",sheet="NHEFS")
features=c("qsmk","alcoholpy","wt82_71","sex","race","age","education","smokeintensity","smokeyrs","exercise","active","wt71")
#categorical with more than 2 levels: education, exercise, active
data=NHEFS[,features]
ind_miss=which(is.na(data$wt82_71)|(data$alcoholpy==2))
data=data[-ind_miss,]
data$qsmk=factor(data$qsmk)
data$alcoholpy=factor(data$alcoholpy)
data$sex=factor(data$sex)
data$race=factor(data$race)
data$education=factor(data$education)
data$exercise=factor(data$exercise)
data$active=factor(data$active)
dim(data)
#table1(wt82_71~.,data)
set.seed(123)
#ind.shuffle=sample(nrow(data))
#data=data[ind.shuffle,]
X0=model.matrix(~-1+.,data=data[,c("sex","race","education","exercise","active","age","smokeintensity","smokeyrs","wt71")])[,-1]
X0=scale(X0)
Y0=data$wt82_71
qsmk=1;alcoholpy=0 #ATE2
qsmk=0;alcoholpy=0 #ATE1
T0=(data$qsmk==qsmk)*(data$alcoholpy==alcoholpy)
mean(T0)
K=5;N=nrow(X0);M=ceiling(N/K);p=ncol(X0);B=20
num_mest=5;num_piest=4
pred_m0=matrix(0,nrow=num_mest*B,ncol=N);pred_m1=matrix(0,nrow=num_mest*B,ncol=N);pred_pi=matrix(0,nrow=num_piest*B,ncol=N)
RF_opt_par=function(x,y,num_iter){
  best_param=list()
  best_loss=Inf
  ptm0=proc.time()
  for (i in 1:num_iter){
    param <- list(ntree=sample(500:2000,1),mtry=sample(1:p,1),nodesize=sample(3:20,1),
                  sample.fraction=runif(1,0.5,0.95),honesty.fraction=runif(1,0.5,0.9),
                  alpha=runif(1,0,0.2),seed.number=sample.int(10000,1)[[1]])
    fit_RF_iter=regression_forest(X=x,Y=y,tune.parameters="none",
                                  num.trees=param$ntree,mtry=param$mtry,min.node.size=param$nodesize,
                                  ci.group.size=1,alpha=param$alpha,seed=param$seed.number,sample.fraction=param$sample.fraction,
                                  honesty.fraction=param$honesty.fraction,honesty.prune.leaves=FALSE)
    loss=mean((fit_RF_iter$predictions-y)^2)
    if (!is.na(loss)&(loss<best_loss)) {
      best_loss=loss
      best_param=param
    }
  }
  fit_RF_best=regression_forest(X=x,Y=y,tune.parameters="none",
                                num.trees=best_param$ntree,mtry=best_param$mtry,min.node.size=best_param$nodesize,
                                ci.group.size=1,alpha=best_param$alpha,seed=best_param$seed.number,sample.fraction=best_param$sample.fraction,
                                honesty.fraction=best_param$honesty.fraction,honesty.prune.leaves=FALSE)
  proc.time()-ptm0
  return(fit_RF_best)
}
opt_param_RF0=matrix(0,nrow=K,ncol=6);opt_param_RF1=matrix(0,nrow=K,ncol=6)
pred_theta=matrix(0,nrow=B,ncol=num_mest*num_piest+1);var_opt=matrix(0,nrow=B,ncol=num_mest*num_piest+1)
grid=expand.grid(est_m=1:num_mest,est_pi=1:num_piest)
mse_m0=rep(0,num_mest);mse_m1=rep(0,num_mest);mse_pi=rep(0,num_piest)
for (b in 1:B){
  ind.shuffle=sample(nrow(data),replace=FALSE)
  X=X0[ind.shuffle,];Y=Y0[ind.shuffle];T=T0[ind.shuffle]
  for (k in 1:K){
    index=((k-1)*M+1):min(k*M,N);Mk=length(index)
    X_nk=X[-index,];Y_nk=Y[-index];T_nk=T[-index]
    #propensity score estimators
    ind_removepi=apply(X_nk,2,sd)==0
    X_nkpi=X_nk[,!ind_removepi];X_kpi=X[index,!ind_removepi]
    pred_pi[(b-1)*num_piest+1,index]=rep(mean(T_nk),Mk) #MCAR
    fit_Log=glm(T_nk~.-pi_N+offset(log(pi_N)),data=data.frame(x=X_nkpi,pi_N=rep(mean(T_nk),N-Mk)),family=binomial)
    pred_pi[(b-1)*num_piest+2,index]=predict(fit_Log,newdata=data.frame(x=X_kpi,pi_N=rep(mean(T_nk),Mk)),type="response") #offset
    lambda.range=sqrt(log(p)/((N-Mk)*mean(T_nk)))*10^seq(from=-3,to=0,length.out=100)
    fit_Loglasso=cv.glmnet(x=X_nkpi,y=T_nk,family="binomial",offset=rep(log(mean(T_nk)),N-Mk),lambda=lambda.range,nfolds=5)
    pred_pi[(b-1)*num_piest+3,index]=predict(fit_Loglasso,newx=X_kpi,type="response",s="lambda.min",newoffset=rep(log(mean(T_nk)),length(index))) #offset
    fit_RF=RF_opt_par(x=X_nkpi,y=T_nk,num_iter=100)
    pred_pi[(b-1)*num_piest+4,index]=predict(fit_RF,newdata=X_kpi)$predictions #RF
    #outcome models
    ind_remove0=apply(X_nk[T_nk==0,],2,sd)==0
    X_nk0=X_nk[T_nk==0,!ind_remove0];X_k0=X[index,!ind_remove0]
    ind_remove1=apply(X_nk[T_nk==1,],2,sd)==0
    ind_poly2=(11:14)[(11:14)%in%which(!ind_remove1)]
    X_nk1=X_nk[T_nk==1,!ind_remove1];X_k1=X[index,!ind_remove1]
    X_nk1_poly=cbind(X_nk1,X_nk[T_nk==1,ind_poly2]^2);X_k1_poly=cbind(X_k1,X[index,ind_poly2]^2)
    #LS
    fit_LS0=lm(y~.,data=data.frame(x=X_nk0,y=Y_nk[T_nk==0]))
    pred_m0[(b-1)*num_mest+1,index]=predict(fit_LS0,newdata=data.frame(x=X_k0))
    fit_LS1=lm(y~.,data=data.frame(x=X_nk1,y=Y_nk[T_nk==1]))
    pred_m1[(b-1)*num_mest+1,index]=predict(fit_LS1,newdata=data.frame(x=X_k1))
    #Lasso
    fit_Lasso0=cv.glmnet(x=X_nk0,y=Y_nk[T_nk==0],family="gaussian",alpha=1,nfolds=5)
    pred_m0[(b-1)*num_mest+2,index]=predict(fit_Lasso0,newx=X_k0,s="lambda.min")
    fit_Lasso1=cv.glmnet(x=X_nk1,y=Y_nk[T_nk==1],family="gaussian",alpha=1,nfolds=5)
    pred_m1[(b-1)*num_mest+2,index]=predict(fit_Lasso1,newx=X_k1,s="lambda.min")
    #poly-Lasso
    fit_polyLasso0=cv.glmnet(x=cbind(X_nk0,X_nk0[,11:14]^2),y=Y_nk[T_nk==0],family="gaussian",alpha=1,nfolds=5)
    pred_m0[(b-1)*num_mest+3,index]=predict(fit_polyLasso0,newx=cbind(X_k0,X_k0[,11:14]^2),s="lambda.min")
    fit_polyLasso1=cv.glmnet(x=X_nk1_poly,y=Y_nk[T_nk==1],family="gaussian",alpha=1,nfolds=5)
    pred_m1[(b-1)*num_mest+3,index]=predict(fit_polyLasso1,newx=X_k1_poly,s="lambda.min")
    #RF
    fit_RF0=RF_opt_par(x=X_nk0,y=Y_nk[T_nk==0],num_iter=100)
    opt_param_RF0[k,]=c(fit_RF0$`_num_trees`,fit_RF0$tunable.params$mtry,fit_RF0$tunable.params$min.node.size,
                        fit_RF0$tunable.params$sample.fraction,fit_RF0$tunable.params$honesty.fraction,fit_RF0$tunable.params$alpha)
    pred_m0[(b-1)*num_mest+4,index]=predict(fit_RF0,newdata=X_k0)$predictions
    fit_RF1=RF_opt_par(x=X_nk1,y=Y_nk[T_nk==1],num_iter=100)
    opt_param_RF1[k,]=c(fit_RF1$`_num_trees`,fit_RF1$tunable.params$mtry,fit_RF1$tunable.params$min.node.size,
                        fit_RF1$tunable.params$sample.fraction,fit_RF1$tunable.params$honesty.fraction,fit_RF1$tunable.params$alpha)
    pred_m1[(b-1)*num_mest+4,index]=predict(fit_RF1,newdata=X_k1)$predictions
    #krls
    fit_krls0=krls(X=X_nk0,y=Y_nk[T_nk==0],whichkernel="gaussian",derivative=FALSE,binary=FALSE,vcov=FALSE)
    pred_m0[(b-1)*num_mest+5,index]=predict(fit_krls0,newdata=X_k0)$fit
    fit_krls1=krls(X=X_nk1,y=Y_nk[T_nk==1],whichkernel="gaussian",derivative=FALSE,binary=FALSE,vcov=FALSE)
    pred_m1[(b-1)*num_mest+5,index]=predict(fit_krls1,newdata=X_k1)$fit
  }
  #check the mse
  #ATE estimator
  pred_theta_b=c();var_opt_b=c()
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    pred_theta0=mean(pred_m0[(b-1)*num_mest+x[1],]+(1-T)*(Y-pred_m0[(b-1)*num_mest+x[1],])/(1-pred_pi[(b-1)*num_piest+x[2],]))
    pred_theta1=mean(pred_m1[(b-1)*num_mest+x[1],]+T*(Y-pred_m1[(b-1)*num_mest+x[1],])/pred_pi[(b-1)*num_piest+x[2],])
    pred_theta_b=c(pred_theta_b,pred_theta1-pred_theta0)
  }
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    psi_mu=(pred_m1[(b-1)*num_mest+x[1],]+T*(Y-pred_m1[(b-1)*num_mest+x[1],])/pred_pi[(b-1)*num_piest+x[2],])-(pred_m0[(b-1)*num_mest+x[1],]+(1-T)*(Y-pred_m0[(b-1)*num_mest+x[1],])/(1-pred_pi[(b-1)*num_piest+x[2],]))-pred_theta[x[1]+num_mest*(x[2]-1)]
    var_opt_b=c(var_opt_b,mean(psi_mu^2))
  }
  pred_theta[b,]=c(mean(Y[T==1])-mean(Y[T==0]),pred_theta_b)
  var_opt[b,]=c(var(Y[T==1])/sum(T)*N+var(Y[T==0])/sum(1-T)*N,var_opt_b)
  #AL_opt[b,]=2*sqrt(var_opt[b,]/N)*qnorm(1-0.05/2)
  #CI_opt[b,]=cbind(pred_theta[b,]-AL_opt[b,]/2,pred_theta[b,]+AL_opt[b,]/2)
  mse_m0=mse_m0+colMeans((t(pred_m0[(b-1)*num_mest+1:num_mest,T==0])-Y[T==0])^2)/B
  mse_m1=mse_m1+colMeans((t(pred_m1[(b-1)*num_mest+1:num_mest,T==1])-Y[T==1])^2)/B
  mse_pi=mse_pi+colMeans((t(pred_pi[(b-1)*num_piest+1:num_piest,])-T)^2)/B
}
med_pred_theta=apply(pred_theta,2,median)
var_med=apply(var_opt+(pred_theta-matrix(rep(med_pred_theta,each=B),nrow=B))^2,2,median)
AL_med=2*sqrt(var_med/N)*qnorm(1-0.05/2)
CI_lower=med_pred_theta-AL_med/2
CI_upper=med_pred_theta+AL_med/2
result=cbind(med_pred_theta,CI_lower,CI_upper,AL_med)
colnames(result)=c("pred_theta","CI_lower","CI_upper","AL")
est_m=factor(grid$est_m,levels=1:num_mest,labels=c("LS","Lasso","poly-Lasso","RF","krls"))
est_pi=factor(grid$est_pi,levels=1:num_piest,labels=c("constant","logistic","log-Lasso","RF"))
rownames(result)=c("empirical",apply(t(1:nrow(grid)),2,function(x){paste(est_pi[x],est_m[x],sep="+")}))
#mse_m0=c(mean((mean(Y[T==0])-Y[T==0])^2),colMeans((t(pred_m0[,T==0])-Y[T==0])^2))
#mse_m1=c(mean((mean(Y[T==1])-Y[T==1])^2),colMeans((t(pred_m1[,T==1])-Y[T==1])^2))
mse_m0=c(mean((mean(Y0[T0==0])-Y0[T0==0])^2),mse_m0)
mse_m1=c(mean((mean(Y0[T0==1])-Y0[T0==1])^2),mse_m1)
#mse_pi=colMeans((t(pred_pi)-T)^2)
result=list(result=result,mse_m0=mse_m0,mse_m1=mse_m1,mse_pi=mse_pi,pihat_N=mean(T))
result
proc.time()-ptm
pdf(paste("qsmk",qsmk,"_alcoholpy",alcoholpy,".pdf",sep=""))
hist(pred_pi[0:(B-1)*num_piest+3,],30,col=rgb(1,0,0,0.5),border=rgb(1,0,0,0.5),xlim=c(0,max(pred_pi)*1.1),xlab="",main="Histogram of predicted propensity scores",freq=FALSE)
hist(pred_pi[0:(B-1)*num_piest+2,],30,add=TRUE,col=rgb(0,0,1,0.5),border=rgb(0,0,1,0.5),freq=FALSE)
hist(pred_pi[0:(B-1)*num_piest+4,],30,add=TRUE,col=rgb(0,1,0,0.5),border=rgb(0,1,0,0.5),freq=FALSE)
dev.off()
filename=paste("qsmk",qsmk,"_alcoholpy",alcoholpy,"_med.RData",sep="")
save.image(filename)

