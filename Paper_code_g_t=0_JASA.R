########libraries to install######
library("LaplacesDemon")
library("fMultivar")
library("evd")
library("condmixt")
library("ismev")
library("mev")
rm(list=ls())

###################inputs################
repl=10000
results_tab=matrix(NA,repl,9)
n=500
m=1000
k=125
p=(1/n) 
x=matrix(NA,n,repl)
y_tilde=matrix(NA,n,repl)
y=matrix(NA,n,repl)
z_tilde=matrix(NA,n,repl)
z=matrix(NA,n,repl)
################Generating simulation data###############

for(le in 1:repl){
  ##############################data preperation##################################
  ### cauchy 
  mu<-matrix(c(0,0,0))
  mu<-t(mu)
  sz<-20000
  s_=0
  r=0
  S=matrix(c(1,s_,s_,s_,1,r,s_,r,1),nrow=3,ncol=3) 
  hh1<-rmvc(n=sz, mu, S)
  zz=matrix(NA,sz,3)
  for (q in 1:sz) {if (hh1[q,1]>0 && hh1[q,2]>0 && hh1[q,3]>0) {zz[q,]=hh1[q,]}}
  hh=na.omit(zz)
  pre1_x=hh[1:n,1]
  #uncomment if cauchy_par=0
  pre2_x=(2/pi)*(atan(pre1_x))
  
  #uncomment if cauchy_par=0.8
  #pre2_x=(1/(pi-atan(0.75)))*(atan(pre1_x)+atan(sqrt((25*(pre1_x^2))+9)/4)-atan(0.75)) 
  
  x1=-log(1-pre2_x)
  
  x[,le]=x1
  y_all=hh[1:(n+m),2]
  y_=hh[1:n,2]
  y[,le]=y_
  
  z_all=hh[1:(n+m),3]
  z_=hh[1:n,3]
  z[,le]=z_
  ###############data_preperation_gamma=0###############
  y3_rank=rank(y_all)
  F_nm<-rep(NA,(n+m))
  for(t in 1:(n+m)){F_nm[t]=y3_rank[t]/(n+m)}
  F_nm1=F_nm-(1/(2*(n+m)))
  y1=-log(1-F_nm1) 
  y_tilde1=y1[1:n]
  y_tilde[,le]=y_tilde1
  z3_rank=rank(z_all)
  Q_nm<-rep(NA,(n+m))
  for(t in 1:(n+m)){Q_nm[t]=z3_rank[t]/(n+m)}
  Q_nm1=Q_nm-(1/(2*(n+m)))
  z1=-log(1-Q_nm1) 
  z_tilde1=z1[1:n]
  z_tilde[,le]=z_tilde1}

##########input for generating the improved esimatiors#########
a_vec=matrix(c(-0.25,-0.125,0,0.125,0.25),nrow=1,ncol=5) #different values of g
reduction_2d=matrix(NA,length(a_vec),3)
reduction_3d=matrix(NA,length(a_vec),3)

for (l in 1:length(a_vec)){
  if(a_vec[l]==0) {for (re in 1:repl){
    
    #########################################improved estimator###########################################
    ###############tail index estimators###########
    #############gamma_1 estimate#########
    x_sort=sort(x[,re])
    ml=gp.fit(x[,re],x_sort[n-k],method="zhang")
    gamma_1=ml$approx.mean[2]
    ############gamma_2 estimate########
    y_tilde_sort=sort(y_tilde[,re])
    ml1=gp.fit(y_tilde[,re],y_tilde_sort[n-k],method="zhang")
    gamma_2=ml1$approx.mean[2]
    ############gamma_3 estimate########
    z_tilde_sort=sort(z_tilde[,re])
    ml2=gp.fit(z_tilde[,re],z_tilde_sort[n-k],method="zhang")
    gamma_3=ml2$approx.mean[2]

    y_sort=sort(y[,re])
    z_sort=sort(z[,re])
    x_rank=rank(x[,re])
    y1_rank=rank(y[,re])
    z1_rank=rank(z[,re])
    
    ########calculation of tail dependence#########
    indicator<-rep(NA,n)
    for (j in 1:n){if (x[j,re]>x_sort[n-k] && y[j,re]>y_sort[n-k]) {indicator[j]=1} else{indicator[j]=0}}
    R1_y<-sum(indicator)/k #R_xy(1,1)
    
    indicator11<-rep(NA,n)
    for (j in 1:n){if (x[j,re]>x_sort[n-k] && z[j,re]>z_sort[n-k]) {indicator11[j]=1} else{indicator11[j]=0}}
    R1_z<-sum(indicator11)/k #R_xz(1,1)
    
    indicator21<-rep(NA,n)
    for (j in 1:n){if (y[j,re]>y_sort[n-k] && z[j,re]>z_sort[n-k]) {indicator21[j]=1} else{indicator21[j]=0}}
    R_y_z<-sum(indicator21)/k #R_yz(1,1)
    
    indicator1=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator1[i]=((n-x_rank[i]+1)/k)^gamma_1} else{indicator1[i]=0}}
    R_s_y_int=(1/gamma_1)*(R1_y-(sum(indicator1)/k)) #R_xy(s,1)/s^{1-\gamma}
    
    indicator2=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator2[i]=((n-x_rank[i]+1)/k)^gamma_1} else{indicator2[i]=0}}
    R_s_z_int=(1/gamma_1)*(R1_z-(sum(indicator2)/k)) #R_xz(s,1)/s^{1-\gamma}
    
    indicator3t=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator3t[i]=(log((n-x_rank[i]+1)/k))} else{indicator3t[i]=0}}
    R_s_xy=(-(sum(indicator3t)/k)) #R_xy(s,1)/s
    
    indicator3ttt=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator3ttt[i]=(log((n-x_rank[i]+1)/k))} else{indicator3ttt[i]=0}}
    R_s_xz=(-(sum(indicator3ttt)/k)) #R_xz(s,1)/s
    
    indicator4t=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator4t[i]=(log((n-y1_rank[i]+1)/k))} else{indicator4t[i]=0}}
    R_t_y=(-(sum(indicator4t)/k)) #R_xy(1,t)/t
    
    indicator4ttt=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator4ttt[i]=(log((n-z1_rank[i]+1)/k))} else{indicator4ttt[i]=0}}
    R_t_z=(-(sum(indicator4ttt)/k)) #R_xz(1,t)/t
    
    indicator3tt=matrix(NA,n,1)
    for(i in 1:n){ if(z[i,re]>z_sort[n-k] && y[i,re]>y_sort[n-k]){indicator3tt[i]=(log((n-y1_rank[i]+1)/k))} else{indicator3tt[i]=0}}
    R_s_yz=(-(sum(indicator3tt)/k)) #R_yz(s,1)/s
    
    indicator4tt=matrix(NA,n,1)
    for(i in 1:n){ if(z[i,re]>z_sort[n-k] && y[i,re]>y_sort[n-k]){indicator4tt[i]=(log((n-z1_rank[i]+1)/k))} else{indicator4tt[i]=0}}
    R_t_yz=(-(sum(indicator4tt)/k)) #R_yz(1,t)/t
    
    ###########improved estimator of the extreme value index################
    ###########2d###########
    R_1=R1_y
    gamma_imp_2d=gamma_1-(R_1*gamma_2)
    
    ###########3d###########
    R_a=R1_y
    R_b=R1_z
    R_ab=R_y_z
    c2=((R_a-(R_ab*R_b))/(1-(R_ab^2)))
    c3=((R_b-(R_ab*R_a))/(1-(R_ab^2)))
    gamma_imp_3d=gamma_1+(c2*(-gamma_2))+(c3*(-gamma_3))
    
    #############sigma_1 estimate#########
    sigma_1=ml$approx.mean[1]
    #############sigma_2 estimate#########
    sigma_2=ml1$approx.mean[1]
    #############sigma_3 estimate#########
    sigma_3=ml2$approx.mean[1]
    
    indicator3=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator3[i]=((log((n-y1_rank[i]+1)/k)^2)/2)} else{indicator3[i]=0}}
    R_lnt_y=(-(sum(indicator3)/k))
    
    indicator3l=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator3l[i]=((log((n-z1_rank[i]+1)/k)^2)/2)} else{indicator3l[i]=0}}
    R_lnt_z=(-(sum(indicator3l)/k))
    
    indicator1d=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator1d[i]=((n-x_rank[i]+1)/k)^gamma_imp_2d} else{indicator1d[i]=0}}
    R_s_y_int_2d=(1/gamma_imp_2d)*(R1_y-(sum(indicator1d)/k)) #R_xy(s,1)/s^{1-\gamma_2d}
    
    
    indicator1dd=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator1dd[i]=((n-x_rank[i]+1)/k)^gamma_imp_3d} else{indicator1dd[i]=0}}
    R_s_y_int_3d=(1/gamma_imp_3d)*(R1_y-(sum(indicator1dd)/k)) #R_xy(s,1)/s^{1-\gamma_3d}
    
    indicator2d=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator2d[i]=((n-x_rank[i]+1)/k)^gamma_imp_3d} else{indicator2d[i]=0}}
    R_s_z_int_3d=(1/gamma_imp_3d)*(R1_z-(sum(indicator2d)/k)) #R_xz(s,1)/s^{1-\gamma_3d}
    
    ###########improved estimator of the sigma################
    ###########2d###########
    R_s=(4*R1_y)-R_t_y-(R_s_xy)
    sigma_imp_2d=sigma_1*(1-((R_s/2)*((sigma_2)-1)))
    
    ###########3d###########
    z12=(4*R1_y)-R_t_y-(R_s_xy)
    z13=(4*R1_z)-R_t_z-(R_s_xz)
    z23=(4*R_y_z)-(R_s_yz+R_t_yz)
    con1=((2*z12)-(z13*z23))/(4-(z23^2))
    con2=((2*z13)-(z12*z23))/(4-(z23^2))
    sigma_imp_3d=sigma_1*(1-((con1*((sigma_2)-1))+(con2*((sigma_3)-1))))
   
    ##################################improved estimator extreme quantile####################################
    ###########standard quantile########## 
    Q11=(sigma_1*((((k/(n*p))^gamma_1)-1)/gamma_1))
    Q1=x_sort[n-k]+Q11 #2d_quant
    
    ###########improved quantile-2d##########
    Q12_2=(sigma_imp_2d*((((k/(n*p))^gamma_imp_2d)-1)/gamma_imp_2d))
    Q_imp_2d=x_sort[n-k]+Q12_2 #2d_with_gam_tilde_sigma_mult
    
    ###########improved quantile-3d##########
    Q12_33=(sigma_imp_3d*((((k/(n*p))^gamma_imp_3d)-1)/gamma_imp_3d)) #3d_with_gam_tilde_sigma_mult
    Q_imp_3d=x_sort[n-k]+Q12_33 #3d_with_gam_tilde_sigma_div
    #########################################################################################
    results_tab[re,1]=gamma_1
    results_tab[re,2]=gamma_imp_2d
    results_tab[re,3]=gamma_imp_3d
    results_tab[re,4]=sigma_1
    results_tab[re,5]=sigma_imp_2d
    results_tab[re,6]=sigma_imp_3d
    results_tab[re,7]=Q1
    results_tab[re,8]=Q_imp_2d
    results_tab[re,9]=Q_imp_3d
    }
    reduction<-matrix(NA,2,3)
    reduction[1,1]=(1-(var(results_tab[,2])/var(results_tab[,1])))*100 #reduction gamma 2d
    reduction[1,2]=(1-(var(results_tab[,5])/var(results_tab[,4])))*100 #reduction sigma 2d
    reduction[1,3]=(1-(var(results_tab[,8])/var(results_tab[,7])))*100 #reduction quantile 2d
    
    reduction[2,1]=(1-(var(results_tab[,3])/var(results_tab[,1])))*100 #reduction gamma 3d
    reduction[2,2]=(1-(var(results_tab[,6])/var(results_tab[,4])))*100 #reduction sigma 3d
    reduction[2,3]=(1-(var(results_tab[,9])/var(results_tab[,7])))*100 #reduction quantile 2d
  } 
  else{a=a_vec[l]
  y_tilde2=matrix(NA,n,repl)
  y_tilde3=matrix(NA,n,repl)
  z_tilde2=matrix(NA,n,repl)
  z_tilde3=matrix(NA,n,repl)
  for (lk in 1:repl){ 
    y_tilde2[,lk]=1-exp(-y_tilde[,lk])
    y_tilde3[,lk]=(1-((1-y_tilde2[,lk])^(-a)))/-a
    z_tilde2[,lk]=1-exp(-z_tilde[,lk])
    z_tilde3[,lk]=(1-((1-z_tilde2[,lk])^(-a)))/-a} 
  
  for (re in 1:repl){
    #########################################improved estimator###########################################
    ###############tail index estimators###########
    #############gamma_1 estimate#########
    x_sort=sort(x[,re])
    ml=gp.fit(x[,re],x_sort[n-k],method="zhang") 
    gamma_1=ml$approx.mean[2]
    ############gamma_2 estimate########
    y_tilde_sort=sort(y_tilde3[,re])
    ml1=gp.fit(y_tilde3[,re],y_tilde_sort[n-k],method="zhang") 
    gamma_2=ml1$approx.mean[2]
    ############gamma_3 estimate########
    z_tilde_sort=sort(z_tilde3[,re])
    ml2=gp.fit(z_tilde3[,re],z_tilde_sort[n-k],method="zhang")
    gamma_3=ml2$approx.mean[2]
    ###########dependence calculations###############
    y_sort=sort(y[,re])
    z_sort=sort(z[,re])
    x_rank=rank(x[,re])
    y1_rank=rank(y[,re])
    z1_rank=rank(z[,re])
    
    indicator<-rep(NA,n)
    for (i in 1:n){if (x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]) {indicator[i]=1} else{indicator[i]=0}}
    R1_y<-sum(indicator)/k #R_xy(1,1)
    
    indicator11<-rep(NA,n)
    for (i in 1:n){if (x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]) {indicator11[i]=1} else{indicator11[i]=0}}
    R1_z<-sum(indicator11)/k #R_xz(1,1)
    
    indicator21<-rep(NA,n)
    for (i in 1:n){if (y[i,re]>y_sort[n-k] && z[i,re]>z_sort[n-k]) {indicator21[i]=1} else{indicator21[i]=0}}
    R_y_z<-sum(indicator21)/k #R_yz(1,1)
    
    indicator1=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator1[i]=((n-x_rank[i]+1)/k)^gamma_1} else{indicator1[i]=0}}
    R_s_y_int=(1/gamma_1)*(R1_y-(sum(indicator1)/k)) #R_xy(s,1)/s^{1-\gamma}
    
    indicator2=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator2[i]=((n-x_rank[i]+1)/k)^gamma_1} else{indicator2[i]=0}}
    R_s_z_int=(1/gamma_1)*(R1_z-(sum(indicator2)/k)) #R_xz(s,1)/s^{1-\gamma}
    
    indicator3=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator3[i]=((n-y1_rank[i]+1)/k)^a} else{indicator3[i]=0}}
    R_t_y_int=(1/a)*(R1_y-(sum(indicator3)/k)) #R_xy(1,t)/t^{1-g}
    
    indicator4=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator4[i]=((n-z1_rank[i]+1)/k)^a} else{indicator4[i]=0}}
    R_t_z_int=(1/a)*(R1_z-(sum(indicator4)/k)) #R_xz(1,t)/t^{1-g}
    
    indicator5=matrix(NA,n,1)
    for(i in 1:n){ if(y[i,re]>y_sort[n-k] && z[i,re]>z_sort[n-k]){indicator5[i]=((n-y1_rank[i]+1)/k)^a} else{indicator5[i]=0}}
    R_s_yz_int=(1/a)*(R_y_z-(sum(indicator5)/k))  #R_yz(s,1)/s^{1-g}
    
    indicator6=matrix(NA,n,1)
    for(i in 1:n){ if(y[i,re]>y_sort[n-k] && z[i,re]>z_sort[n-k]){indicator6[i]=((n-z1_rank[i]+1)/k)^a} else{indicator6[i]=0}}
    R_t_zy_int=(1/a)*(R_y_z-(sum(indicator6)/k)) #R_yz(1,t)/t^{1-g}
    
    indicator3t=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator3t[i]=(log((n-x_rank[i]+1)/k))} else{indicator3t[i]=0}}
    R_s_xy=(-(sum(indicator3t)/k)) #R_xy(s,1)/s
    
    indicator3ttt=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator3ttt[i]=(log((n-x_rank[i]+1)/k))} else{indicator3ttt[i]=0}}
    R_s_xz=(-(sum(indicator3ttt)/k)) #R_xz(s,1)/s
    
    indicator4t=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator4t[i]=(log((n-y1_rank[i]+1)/k))} else{indicator4t[i]=0}}
    R_t_y=(-(sum(indicator4t)/k)) #R_xy(1,t)/t
    
    indicator4ttt=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator4ttt[i]=(log((n-z1_rank[i]+1)/k))} else{indicator4ttt[i]=0}}
    R_t_z=(-(sum(indicator4ttt)/k)) #R_xz(1,t)/t
    
    indicator3tt=matrix(NA,n,1)
    for(i in 1:n){ if(z[i,re]>z_sort[n-k] && y[i,re]>y_sort[n-k]){indicator3tt[i]=(log((n-y1_rank[i]+1)/k))} else{indicator3tt[i]=0}}
    R_s_yz=(-(sum(indicator3tt)/k)) #R_yz(s,1)/s
    
    indicator4tt=matrix(NA,n,1)
    for(i in 1:n){ if(z[i,re]>z_sort[n-k] && y[i,re]>y_sort[n-k]){indicator4tt[i]=(log((n-z1_rank[i]+1)/k))} else{indicator4tt[i]=0}}
    R_t_yz=(-(sum(indicator4tt)/k)) #R_yz(1,t)/t
    
    ###########improved extreme value index estimator################
    ############2d#############
    R_1=R1_y+(((a)/(a+1))*((R_s_xy-(((2*a)+1)*R_t_y_int))))
    gamma_imp_2d=gamma_1+(((1)/(a+1))*R_1*(a-gamma_2))
    ############3d#############
    R_a=(R1_y-((a/(a+1))*(R_s_xy-(((2*a)+1)*R_t_y_int))))
    R_b=(R1_z-((a/(a+1))*(R_s_xz-(((2*a)+1)*R_t_z_int))))
    R_ab=R_y_z
    g2=(1+a)
    c2=((1)/(1+a))*((R_a-(R_ab*R_b))/(1-(R_ab^2)))
    c3=((1)/(1+a))*((R_b-(R_ab*R_a))/(1-(R_ab^2)))  
    gamma_imp_3d=gamma_1+(c2*(a-gamma_2))+(c3*(a-gamma_3))
    
    #########################################improved estimator sigma###########################################
    ###############sigma estimators###########
    #############sigma_1 estimate#########
    sigma_1=ml$approx.mean[1]
    #############sigma_2 estimate#########
    sigma_2=ml1$approx.mean[1]
    #############sigma_3 estimate#########
    sigma_3=ml2$approx.mean[1]
    
    indicator3=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && y[i,re]>y_sort[n-k]){indicator3[i]=((log((n-x_rank[i]+1)/k)^2)/2)} else{indicator3[i]=0}}
    R_lnt_y=(-(sum(indicator3)/k))
    
    indicator3l=matrix(NA,n,1)
    for(i in 1:n){ if(x[i,re]>x_sort[n-k] && z[i,re]>z_sort[n-k]){indicator3l[i]=((log((n-x_rank[i]+1)/k)^2)/2)} else{indicator3l[i]=0}}
    R_lnt_z=(-(sum(indicator3l)/k))
    #########################################improved estimator sigma###########################################
    #############2d#############
    R_s=(((3*a)-1)*R_s_xy)+(a*R_lnt_y)-((((2*a)+1)^2)*R_t_y_int)+((2*(a+2))*R1_y)
    #c=(R_s/2)
    sigma_imp_2d=sigma_1*(1-((R_s/(1+((1+a)^2)))*((sigma_2/((n/k)^a))-1)))
    
    #############3d#############
    z12=(((3*a)-1)*R_s_xy)+(a*R_lnt_y)-((((2*a)+1)^2)*R_t_y_int)+((2*(a+2))*R1_y)
    z13=(((3*a)-1)*R_s_xz)+(a*R_lnt_z)-((((2*a)+1)^2)*R_t_z_int)+((2*(a+2))*R1_z)
    z23=(((a^2)+(4*a)+4)*R_y_z)-((a+1)*(R_s_yz+R_t_yz))
    g6= (1+((1+a)^2))
    con1=((g6*z12)-(z13*z23))/((g6^2)-(z23^2))
    con2=((g6*z13)-(z12*z23))/((g6^2)-(z23^2))
    sigma_imp_3d=sigma_1*(1-(con1*((sigma_2/(n/k)^a)-1)+(con2*((sigma_3/(n/k)^a)-1))))
    
    ##################################improved estimator extreme quantile####################################
    ###########standard quantile##########    
    Q11=(sigma_1*((((k/(n*p))^gamma_1)-1)/gamma_1))
    Q1=x_sort[n-k]+Q11 
    
    ###########improved quantile-2d##########
    Q12_2=(sigma_imp_2d*((((k/(n*p))^gamma_imp_2d)-1)/gamma_imp_2d))
    Q_imp_2d=x_sort[n-k]+Q12_2 
 
    ###########improved quantile-3d##########
    Q12_33=(sigma_imp_3d*((((k/(n*p))^gamma_imp_3d)-1)/gamma_imp_3d)) #3d_with_gam_tilde_sigma_mult
    Q_imp_3d=x_sort[n-k]+Q12_33
    
    #########################################################################################
    results_tab[re,1]=gamma_1
    results_tab[re,2]=gamma_imp_2d
    results_tab[re,3]=gamma_imp_3d
    results_tab[re,4]=sigma_1
    results_tab[re,5]=sigma_imp_2d
    results_tab[re,6]=sigma_imp_3d
    results_tab[re,7]=Q1
    results_tab[re,8]=Q_imp_2d
    results_tab[re,9]=Q_imp_3d
    }
  reduction<-matrix(NA,2,3)
  reduction[1,1]=(1-(var(results_tab[,2])/var(results_tab[,1])))*100 #reduction gamma 2d
  reduction[1,2]=(1-(var(results_tab[,5])/var(results_tab[,4])))*100 #reduction sigma 2d
  reduction[1,3]=(1-(var(results_tab[,8])/var(results_tab[,7])))*100 #reduction quantile 2d
  
  reduction[2,1]=(1-(var(results_tab[,3])/var(results_tab[,1])))*100 #reduction gamma 3d
  reduction[2,2]=(1-(var(results_tab[,6])/var(results_tab[,4])))*100 #reduction sigma 3d
  reduction[2,3]=(1-(var(results_tab[,9])/var(results_tab[,7])))*100 #reduction quantile 2d
  
  }
  reduction_2d[l,]=reduction[1,]
  reduction_3d[l,]=reduction[2,]}

simulation_res_2d=data.frame(cbind(t(a_vec),reduction_2d))
simulation_res_2d=setNames(simulation_res_2d,c("g","red_gamma","red_sigma","red_quantile"))
simulation_res_3d=data.frame(cbind(t(a_vec),reduction_3d))
simulation_res_3d=setNames(simulation_res_3d,c("g","red_gamma","red_sigma","red_quantile"))
