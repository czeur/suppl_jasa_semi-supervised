library("LaplacesDemon")
library("fMultivar")
library("evd")
library("condmixt")
library("ismev")
library("mev")

rm(list=ls())
repl=500
gamma=matrix(NA,repl,2)
mse=matrix(NA,repl,3)

###################inputs################
n=500
m=1000
g_t=0
x=matrix(NA,n,repl)
y_tilde=matrix(NA,n,repl)
y=matrix(NA,n,repl)

for(le in 1:repl){
  ##############################data preperation##################################
  sz=100000
  hh1<-rcauchy2d(sz, rho =0.8)
  zz=matrix(NA,sz,2)
  for (q in 1:sz) {if (hh1[q,1]>0 && hh1[q,2]>0) {zz[q,]=hh1[q,]}}
  hh=na.omit(zz)
  ##################x_prep###############
  pre1_x=hh[1:n,1]
  
  # rho=0.8
  pre2_x=(1/(pi-atan(0.75)))*(atan(pre1_x)+atan(sqrt((25*(pre1_x^2))+9)/4)-atan(0.75)) 
  #if g_t=! 0 
  #x1=(1-((1-pre2_x)^(-g_t)))/-g_t
  #if g_t=0
  x1=-log(1-pre2_x)
  x[,le]=x1
  #################y_prep###############
  y_all=hh[1:(n+m),2]
  y_=hh[1:n,2]
  y[,le]=y_
  ###############gamma=0###############
  y3_rank=rank(y_all)
  F_nm<-rep(NA,(n+m))
  for(t in 1:(n+m)){F_nm[t]=y3_rank[t]/(n+m)}
  F_nm1=F_nm-(1/(2*(n+m)))
  y1=-log(1-F_nm1) 
  y_tilde1=y1[1:n]
  y_tilde[,le]=y_tilde1}

for(k in 1:499){
  for (re in 1:repl){
    
    #########################################improved estimator###########################################
    ###############tail index estimators###########
    
    x_sort=sort(x[,re])
    ml=gp.fit(x[,re],x_sort[n-k],method="zhang")
    gamma_1=ml$approx.mean[2] #gamma_1_estimate
    
    y_tilde_sort=sort(y_tilde[,re])
    ml1=gp.fit(y_tilde[,re],y_tilde_sort[n-k],method="zhang")
    gamma_2=ml1$approx.mean[2] #gamma_2_estimate
    ###########chen(0)###############
    y_sort=sort(y[,re])
    x_rank=rank(x[,re])
    y1_rank=rank(y[,re])
    indicator<-rep(NA,n)
    for (j in 1:n){if (x[j,re]>x_sort[n-k] && y[j,re]>y_sort[n-k]) {indicator[j]=1} else{indicator[j]=0}}
    R1_1<-sum(indicator)/k
    
    indicator1=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator1[i]=((n-x_rank[i]+1)/k)^gamma_1} else{indicator1[i]=0}}
    R1_intg=(1/gamma_1)*(R1_1-(sum(indicator1)/k)) 
    
    indicator2=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator2[i]=(log((n-y1_rank[i]+1)/k))} else{indicator2[i]=0}}
    R2_intg=(-(sum(indicator2)/k))
    g1=(gamma_1/(gamma_1+1))
    g2=((2*gamma_1)+1)
    R_1=(g1*((g2*R1_intg)-R2_intg))-R1_1
    ###########Adapted estimator################
    gamma_1_tilde=gamma_1+((gamma_1+1)*R_1*gamma_2)
    #########################################################################################
    gamma[re,1]=(gamma_1-g_t)^2
    gamma[re,2]=(gamma_1_tilde-g_t)^2
  }
  mse[k,1]=k
  mse[k,2]=sum(gamma[1:500,1])/500
  mse[k,3]=sum(gamma[1:500,2])/500
}
rmse=matrix(NA,499,2)
rmse[,1]=sqrt(mse[1:499,2])
rmse[,2]=sqrt(mse[1:499,3])


plot(mse[50:499,1],rmse[50:499,1],type='l', lty=2, lwd=1,xlab='k',ylim=c(0.03,0.16),ylab='RMSE',cex.lab=1.4)
lines(rmse[50:499,2],col="green", type="l", lty=3, lwd=2)
legend("topright", legend=c("MLE", "AMLE"),col=c("black","green"), lty=2:3, cex=0.8)





