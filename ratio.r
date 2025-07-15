#source("~/CODES/vizu_bivariate_ratio.R")
rm(list = ls())
library(plot3D)
library(numDeriv)
library(xtable)
#
#Directory for output figures files
directory_="~/Documents/TEST/OUTPUT_RATIO/";

resultats = matrix(0,nrow=4,ncol=50);

## ----------------------------------------------------------------
dd_ = 0;
for (data_ in c("D1","D2","D3","D4")) {
  names_res = NULL;
  dd_ = dd_+1;
  
  cat(("###################################################"))
  
  #as = seq(-7,-0.05,length.out = 5000);
  #as = seq(-1.2,-0.05,length.out = 5000);
  #as = seq(0.05,+1.2,length.out = 5000);
  as_left = seq(-1.2,-0.05,length.out = 5000);
  as_right = seq(0.05,+1.2,length.out = 5000);
  #as_left = seq(-3.5,-0.05,length.out = 5000);
  #as_right = seq(0.05,+3.5,length.out = 5000);
  as = c(as_left, as_right);
  as00=as;
  
  #data_ = "D1"
  #data_ = "D2"
  #data_ = "D3"
  #data_ = "D4"
  #data_ = "DI"
  
  ##
  if (data_=="D1") {
    N=100
    n=29
    Y_bar   = 2.3640
    X1_bar  = 2.9250
    X2_bar  = 5.2390
    Cy_2    = 2.5582
    Cx1_2   = 0.0661
    Cx2_2   = 0.0461
    rhoyx1  = 0.1602
    rhoyx2  = 0.0829
    rhox1x2 = 0.0846
  } else if (data_=="D2") {
    N=97
    n=30
    Y_bar  = 3135.6186
    X1_bar = 3050.2784
    X2_bar = 2743.9587
    Cy_2 = 4.8674
    Cx1_2 = 5.4812
    Cx2_2 = 6.2422
    rhoyx1 = 0.8072
    rhoyx2 = 0.8501
    rhox1x2 = 0.6122
  } else if (data_=="D3") {
    N=332
    n=80
    Y_bar  = 1093.1
    X1_bar = 181.57
    X2_bar = 143.31
    Cy_2 = 0.7626
    Cx1_2 = 0.7684
    Cx2_2 = 0.7616
    rhoyx1 = 0.973
    rhoyx2 = 0.862
    rhox1x2 = 0.842
    Sy_2   = Cy_2*Y_bar^2
    Sx1_2  = Cx1_2 * X1_bar^2
    Sx2_2  = Cx2_2 * X2_bar^2
    Sx1y   = rhoyx1 * sqrt(Sx1_2) * sqrt(Sy_2)
    Sx2y   = rhoyx2 * sqrt(Sx2_2) * sqrt(Sy_2)
    Sx1x2  = rhox1x2 * sqrt(Sx1_2) * sqrt(Sx2_2)
  } else if (data_=="D4") {
    N=18;
    n=8;
    Y_bar  = 13.797
    X1_bar = 2.444
    X2_bar = 38.444
    Sy_2   = 35.486
    Sx1_2  = 74.679
    Sx2_2  = 174.967
    Sx1y   = 42.262
    Sx2y   = 46.512
    Sx1x2  = 11.546
    rhoyx1 = Sx1y   / sqrt(Sx1_2) / sqrt(Sy_2)
    rhoyx2 = Sx2y   / sqrt(Sx2_2) / sqrt(Sy_2)
    rhox1x2 = Sx2y  / sqrt(Sx2_2) / sqrt(Sx1_2)
    Cy_2 = Sy_2 / Y_bar^2
    Cx1_2 = Sx1_2 / X1_bar^2
    Cx2_2 = Sx2_2 / X2_bar^2
  } else if (data_=="DI") {
    N=200
    n=50
    Y_bar  = 500
    X1_bar = 25
    X2_bar = 30
    Cy_2 = 5
    Cx1_2 = 2
    Cx2_2 = 3
    rhoyx1 = 0.1
    rhoyx2 = 0.6
    rhox1x2 = 0.5
  }
  Sy_2   = Cy_2*Y_bar^2
  Sx1_2  = Cx1_2 * X1_bar^2
  Sx2_2  = Cx2_2 * X2_bar^2
  Sx1y   = rhoyx1 * sqrt(Sx1_2) * sqrt(Sy_2)
  Cyx1=Sx1y/Y_bar/X1_bar;
  Sx2y   = rhoyx2 * sqrt(Sx2_2) * sqrt(Sy_2)
  Cyx2=Sx2y/Y_bar/X2_bar;
  Sx1x2  = rhox1x2 * sqrt(Sx1_2) * sqrt(Sx2_2)
  Cx1x2=Sx1x2/X1_bar/X2_bar;
  f  = n/N;
  lb = (N-n)/n/N;
  cat("data=",data_,"\n");
  cat("N =",N,"|","n =",n,"|","f =",f,"|","Y_bar =",Y_bar,"|","X1_bar =",
      X1_bar,"|","X2_bar =",X2_bar,"Sy_2 =",Sy_2,"\n");
  cat("Sx1_2 =",Sx1_2,"|","Sx2_2 =",Sx2_2,"|","Sx1x2 =",Sx1x2,"|","Sx1y =",
      Sx1y,"|","Sx2y =",Sx2y,"\n");
  
  ## ----------------------------------------------------------------
  
  alpha=0;
  a=-1/2;
  b=3/8-alpha/4;
  
  Y_RPR03 <- function(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) {
    A = Y_bar^2 + Y_bar^2*lb*(Cy_2+(a^2+2*b)*Cx2_2+4*a*Cyx2);
    B = lb*Cx1_2;
    C = 2*Y_bar*lb*(a*Cx1x2+0.5*Cyx1);
    D0 = +Y_bar^2 + Y_bar^2*lb*(b*Cx2_2+a*Cyx2);
    D1 = +a*Y_bar*lb*Cx1x2;
    E = Y_bar^2;
    
    k0_ = (B*D0-C*D1)/(A*B-C^2);
    k1_ = (A*D1-C*D0)/(A*B-C^2);
    
    MSE_ = k0_^2*A+k1_^2*B+2*k0_*k1_*C-2*k0_*D0-2*k1_*D1+E;
    BIAS_ = k0_*Y_bar*( (k0_-1)/k0_ + b*lb*Cx2_2+a*lb*Cyx2 ) + k1_*a*lb*Cx1x2;
    return(list(MSE=MSE_,BIAS=BIAS_,k0=k0_, k1=k1_));
  }
  
  Y_RPR03_k01 <- function(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) {
    A = Y_bar^2 + Y_bar^2*lb*(Cy_2+(a^2+2*b)*Cx2_2+4*a*Cyx2);
    B = lb*Cx1_2;
    C = 2*Y_bar*lb*(a*Cx1x2+0.5*Cyx1);
    D0 = +Y_bar^2 + Y_bar^2*lb*(b*Cx2_2+a*Cyx2);
    D1 = +a*Y_bar*lb*Cx1x2;
    E = Y_bar^2;
    
    k0_ = 1;
    k1_ = (D1-C)/B;
    
    MSE_ = k0_^2*A+k1_^2*B+2*k0_*k1_*C-2*k0_*D0-2*k1_*D1+E;
    #cat("diffMSE=",MSE_-(A-2*D0+E-(D1-C)^2/B),"\n");
    BIAS_ = k0_*Y_bar*( (k0_-1)/k0_ + b*lb*Cx2_2+a*lb*Cyx2 ) + k1_*a*lb*Cx1x2;
    return(list(MSE=MSE_,BIAS=BIAS_,k0=k0_, k1=k1_));
  }
  
  Y_MSE_MIN_SUMf1f2 <- 
    function(a1,b1,a2,b2,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) { # k1, k2 no constrained
      A = a1^2*Cx1_2;
      B = a2^2*Cx2_2;
      C = a1*a2*Cx1x2;
      D1 = -a1*Cyx1;
      D2 = -a2*Cyx2;
      E = Cy_2;
      
      k1 = (B*D1-C*D2)/(A*B-C^2);
      k2 = (A*D2-C*D1)/(A*B-C^2);
      
      MSE_ = 
        Y_bar^2 * lb * (k1^2*A+k2^2*B+2*k1*k2*C-2*k1*D1-2*k2*D2+E );
      MSE2_ = 
        Y_bar^2 * lb * ( Cy_2 - (B*D1^2-2*C*D1*D2+A*D2^2)/(A*B-C^2) );
      MSE3_ = Y_bar^2 * lb * 
        ( Cy_2 - (Cx2_2*Cyx1^2-2*Cx1x2*Cyx1*Cyx2+Cx1_2*Cyx2^2)/(Cx1_2*Cx2_2-Cx1x2^2) );
      BIAS_ = Y_bar * lb * 
        ( k1* (a1*Cyx1 + b1*Cx1_2 ) + k2* (a2*Cyx2 + b2*Cx2_2 ) ) ;
      return(list(MSE=MSE_,MSE2_=MSE2_,MSE3_=MSE3_,BIAS=BIAS_,alpha1=k1,alpha2=k2));
    }
  
  Y_MSE_MIN_SUMf1f2_extendolkin <- 
    function(a1,b1,a2,b2,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) { # k1+k2=1
      A = a1^2*Y_bar^2*lb*Cx1_2;
      B = a2^2*Y_bar^2*lb*Cx2_2;
      C = a1*a2*Y_bar^2*lb*Cx1x2;
      D1 = -a1*Y_bar^2*lb*Cyx1;
      D2 = -a2*Y_bar^2*lb*Cyx2;
      E = Y_bar^2*lb*Cy_2;
      
      U = a2*Cyx2-a1*Cyx1+a2^2*Cx2_2-a1*a2*Cx1x2;
      V = a1^2*Cx1_2+a2^2*Cx2_2-2*a1*a2*Cx1x2;
      
      k1 = U/V;
      k2 = 1-k1;
      
      MSE_ = ( k1^2*A+k2^2*B+2*k1*k2*C-2*k1*D1-2*k2*D2+E );
      BIAS_ = Y_bar * lb * ( k1* (a1*Cyx1 + b1*Cx1_2 ) + 
                               k2* (a2*Cyx2 + b2*Cx2_2 ) ) ;
      return(list(MSE=MSE_,BIAS=BIAS_,alpha1=k1,alpha2=k2));
    }
  
  plotpdf <- function(xs,ys,xlab_,abval_,namefile_) {
    pdf(namefile_);
    par(lwd=3)
    plot(xs,ys,type='l',xlab=xlab_, ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
    abline(v=abval_,col="red");
    dev.off();
  }
  
  ## ----------------------------------------------------------------
  
  resu2_=Y_RPR03(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
  
  MSE_CRR_ab2 = resu2_$MSE;
  
  V_YO        = lb*Y_bar^2*Cy_2
  
  cat("V_YO/MSE_CRR_ab2=",V_YO/MSE_CRR_ab2,"\n");
  
  ## ----------------------------------------------------------------
  
  # f = 1 + a d_x + b d_x^2
  #as = seq(-0.7,-0.05,length.out = 5000);
  #as = seq(-3.0,-0.05,length.out = 5000);
  l=length(as);
  mses1=numeric(l);
  bias1=numeric(l);
  t=0;
  for (a in as) {
    t=t+1;
    resu_ = Y_RPR03(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
    mses1[t] = resu_$MSE;
    bias1[t] = resu_$BIAS;
  }
  as1=as;
  tmin = which.min(mses1);
  MSE_min1 = mses1[tmin]
  cat("V_YO/MSE_min1=",V_YO/MSE_min1,"\n");
  aopt=as[tmin];
  cat("aopt=",as[tmin],"  ( b=",b,")\n");
  
  ##section surface
  # f = 1 + a d_x + b d_x^2
  vect_bs = c(0.05,0.15,0.25,0.35,0.45,0.55);
  vect_mses1 = list();
  vect_bias1 = list();
  cc=0;
  for (bb in vect_bs) {
    cc=cc+1;
    l_bb=length(as);
    mses1_bb=numeric(l_bb);
    bias1_bb=numeric(l_bb);
    t=0;
    for (a in as) {
      t=t+1;
      resu__bb = Y_RPR03(a,bb,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
      mses1_bb[t] = resu__bb$MSE;
      bias1_bb[t] = resu__bb$BIAS;
    }
    vect_mses1[[as.character(cc)]] = mses1_bb;
    vect_bias1[[as.character(cc)]] = bias1_bb;
  }
  
  ## ----------------------------------------------------------------
  ##k0=1
  # f = 1 + a d_x + b d_x^2
  #as = seq(-0.7,-0.05,length.out = 5000);
  l=length(as);
  mses1_k01=numeric(l);
  bias1_k01=numeric(l);
  t=0;
  for (a in as) {
    t=t+1;
    resu_k01_ = 
      Y_RPR03_k01(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
    mses1_k01[t] = resu_k01_$MSE;
    bias1_k01[t] = resu_k01_$BIAS;
  }
  as1=as;
  tmin_k01 = which.min(mses1_k01);
  MSE_min1_k01 = mses1[tmin_k01]
  cat("V_YO/MSE_min1_k01=",V_YO/MSE_min1_k01,"\n");
  aopt_k01=as[tmin_k01];
  cat("aopt_k01=",as[tmin_k01],"  ( b=",b,")\n");
  cat("aopt_k01_formula=", (Cyx1*Cx1x2-Cx1_2*Cyx2)/(Cx1_2*Cx2_2-Cx1x2^2),
      " ( b=",b,")\n");
  aopt_k01_formula = (Cyx1*Cx1x2-Cx1_2*Cyx2)/(Cx1_2*Cx2_2-Cx1x2^2);
  
  ## ----------------------------------------------------------------
  
  # f = (X+tau)/(x+tau)
  #gs = seq(3/7*X2_bar,9.5*X2_bar,length.out = 5000);
  gs_left = seq( (-1-(min(as_left)))/(min(as_left))*X2_bar,(-1-(max(as_left)))/(max(as_left))*X2_bar,length.out = 5000);
  gs_right = seq( (-1-(min(as_right)))/(min(as_right))*X2_bar,(-1-(max(as_right)))/(max(as_right))*X2_bar,length.out = 5000);
  gs = c(gs_left, gs_right);
  l=length(gs);
  mses3=numeric(l);
  bias3=numeric(l);
  as=numeric(l);
  t=0;
  for (g in gs) {
    t=t+1;
    a = -X2_bar/(X2_bar+g);
    b = ( (X2_bar)/(X2_bar+g) )^2;
    resu_    = Y_RPR03(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
    mses3[t] = resu_$MSE;
    bias3[t] = resu_$BIAS;
    as[t]    = a;
  }
  as2=as;
  tmintau = which.min(mses3);
  MSE_min3 = mses3[tmintau]
  cat("V_YO/MSE_min3=",V_YO/MSE_min3,"\n");
  cat("gammaopt_tau=",gs[tmintau],"\n");
  aopt_tau=as2[tmintau];
  cat("aopt_tau=",as2[tmintau],"\n");
  MSE_min3_b0.375 = 
    Y_RPR03(aopt_tau,0.375,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)$MSE;
  cat("V_YO/MSE_min3_b0.375=",V_YO/MSE_min3_b0.375,"\n");
  
  as = as1;
  
  ## ----------------------------------------------------------------
  as=as00;
  # resu_extolk_ = 
  #   Y_MSE_MIN_SUMf1f2_extendolkin(-0.5,1,-0.5,1,
  #                                 Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
  # MSE_extolk   = resu_extolk_$MSE;
  # BIAS_extolk  = resu_extolk_$BIAS;
  # cat("V_YO/MSE_extolk=",V_YO/MSE_extolk,"\n");
  
  # alpha*f1+(1-alpha)*f2
  #f1=f2
  #as = seq(-0.7,-0.05,length.out = 5000);
  #as = seq(-3.0,-0.05,length.out = 5000);
  l=length(as);
  mses5=numeric(l);
  bias5=numeric(l);
  t=0;
  for (a in as) {
    t=t+1;
    resu_ = Y_MSE_MIN_SUMf1f2_extendolkin(a,1,a,1,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
    mses5[t] = resu_$MSE;
    bias5[t] = resu_$BIAS;
  }
  as5=as;
  tmin5 = which.min(mses5);
  MSE_min5 = mses5[tmin5]
  cat("V_YO/MSE_min5=",V_YO/MSE_min5,"\n");
  aopt5=as[tmin5];
  cat("aopt5=",as[tmin5],"  ( b1=b2=",1,")\n");
  
  ## ----------------------------------------------------------------
  
  divmse_a_b0.375 <- function(x) {# "b" is set to 0.375
    return(V_YO/Y_RPR03(x,0.375,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)$MSE);
  }
  divmse_b <- function(x) {# "a" is a global var
    return(V_YO/Y_RPR03(a,x,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)$MSE);
  }
  
  cat("numDeriv::grad(divmse_a_b0.375,x=-0.5)=",numDeriv::grad(divmse_a_b0.375,x=-0.5),"\n");
  cat("numDeriv::grad(divmse_a_b0.375,x=aopt)=",numDeriv::grad(divmse_a_b0.375,x=aopt),"   (aopt=", aopt,"\n");
  cat("numDeriv::grad(divmse_a_b0.375,x=aopt)=",numDeriv::grad(divmse_a_b0.375,x=aopt_k01_formula),"   (aopt=", aopt_k01_formula,"\n");
  a=-0.5;
  cat("numDeriv::grad(divmse_b,x=0.375)=",numDeriv::grad(divmse_b,x=0.375), "a=",a,"\n");
  a=aopt;
  cat("numDeriv::grad(divmse_b,x=0.375)=",numDeriv::grad(divmse_b,x=0.375), "a=aopt=",a,"\n");
  
  #
  resultats[dd_,1] = round(numDeriv::grad(divmse_a_b0.375,x=-0.5),8);
  names_res <- c(names_res, "deriv1");
  #
  resultats[dd_,2] = round(numDeriv::grad(divmse_a_b0.375,x=aopt_k01_formula),8);
  names_res <- c(names_res, "deriv2");
  #
  a=-0.5;
  resultats[dd_,3] = round(numDeriv::grad(divmse_b,x=0.375),8);
  names_res <- c(names_res, "deriv3_a=-0.5");
  #
  a=aopt_k01_formula;
  resultats[dd_,4] = round(numDeriv::grad(divmse_b,x=0.375),8);
  names_res <- c(names_res, "deriv4_a=aopt_k01_formula");
  # 
  
  ## ----------------------------------------------------------------
  
  MSE_yreg = lb*Sy_2*(1-rhoyx1^2-rhoyx2^2+2*rhoyx1*rhoyx2*rhox1x2);
  cat("V_YO/MSE_yreg=",V_YO/MSE_yreg,"\n");
  
  MSE_ydif = lb*Sy_2*(1-(rhoyx1^2+rhoyx2^2-2*rhoyx1*rhoyx2*rhox1x2)/(1-rhox1x2^2));
  cat("V_YO/MSE_ydif=",V_YO/MSE_ydif,"\n");   
  
  #MSE_ydif = lb*Sy_2*(1-(rhoyx1^2+rhoyx2^2-2*rhoyx1*rhoyx2*rhox1x2)/(1-rhox1x2^2));
  #cat("V_YO/MSE_ydif=",V_YO/MSE_ydif,"\n");
  L = 1-rhox1x2^2-(rhoyx1^2+rhoyx2^2-2*rhoyx1*rhoyx2*rhox1x2);
  A = 1-rhox1x2^2
  dt = (1-f)/n
  
  #Efficient Estimator of a Finite Population Mean Using Two Auxiliary Variables and 
  #Numerical Application in Agricultural, Biomedical, and Power Engineering
  #y_RE = w 1 y + w 2 (X 1 − x 1 ) + w 3 (X 2 − x 2 )
  MSE_yRE = dt*Sy_2*L/(A+dt*Cy_2*L);
  cat("V_YO/MSE_yRE=",V_YO/MSE_yRE,"\n");   
  
  MSE_ypr = 0.25*dt*Y_bar^2*(4*Cy_2*L-dt*A*(Cx1_2+Cx2_2)^2)/(A+dt*Cy_2*L+dt*A*(Cx1_2+Cx2_2));
  cat("V_YO/MSE_ypr=",V_YO/MSE_ypr,"\n");   
  
  ## ----------------------------------------------------------------
  
  resu_add12_ = 
    Y_MSE_MIN_SUMf1f2(-0.5,1,-0.5,1,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)
  MSE_add12 =resu_add12_$MSE;
  BIAS_add12 =resu_add12_$BIAS;
  cat("V_YO/MSE_add12=",V_YO/MSE_add12,"\n");
  
  resu_extolk_ = 
    Y_MSE_MIN_SUMf1f2_extendolkin(-0.5,1,-0.5,1,
                                  Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
  MSE_extolk   = resu_extolk_$MSE;
  BIAS_extolk  = resu_extolk_$BIAS;
  cat("V_YO/MSE_extolk=",V_YO/MSE_extolk,"\n");
  
  b=3/8;
  MSE_improved_Rcm = 
    Y_RPR03(aopt_k01_formula,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)$MSE;
  cat("V_YO/MSE_improved_Rcm=",V_YO/MSE_improved_Rcm,"\n");
  
  ## ----------------------------------------------------------------
  
  Cy = sqrt(Cy_2);
  Cx = sqrt(Cx1_2);
  Cz = sqrt(Cx1_2);
  Cxy = Cyx1;
  Cyz = Cyx2;
  Cxz = Cx1x2
  bopt = (Cxz^4*Cy^2-2*Cxy*Cxz^3*Cyz+Cx^2*Cxz^2*Cyz^2+
            (-2*Cx^4*Cyz^2+4*Cx^2*Cxy*Cxz*Cyz-2*Cx^2*Cxz^2*Cy^2)*Cz^2+
            (Cx^4*Cy^2-Cx^2*Cxy^2)*Cz^4)/(Cx^4*Cz^6-2*Cx^2*Cxz^2*Cz^4+Cxz^4*Cz^2)
  cat("Data = ", data_, "  bopt = ", bopt, "\n");
  
  resu_improved_a_b_Rcm = 
    Y_RPR03(aopt_k01_formula,bopt,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)
  
  MSE_improved_a_b_Rcm = 
    Y_RPR03(aopt_k01_formula,bopt,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)$MSE;
  cat("V_YO/MSE_improved_a_b_opt_Rcm=",V_YO/MSE_improved_a_b_Rcm,"\n");
  cat("k0=",resu_improved_a_b_Rcm$k0," ", "k1=,",resu_improved_a_b_Rcm$k1);
  cat("\n");
  
  ## ----------------------------------------------------------------
  library(latex2exp)
  #ymin_=-5;
  #ymax_=30;
  ymin_=V_YO/max(mses1,mses1_k01,MSE_add12,MSE_extolk,mses3,mses5); #
  ymax_=V_YO/min(mses1,mses1_k01,MSE_add12,MSE_extolk,mses3,mses5); #
  #namefile_=paste(directory_,"PRE_ALL_",data_,".pdf",sep="");
  #pdf(namefile_);
  namefile_=paste(directory_,"PRE_ALL_",data_,".png",sep="");
  png(namefile_);
  par(lwd=3)
  plot(as1,V_YO/mses1,type='l',xlab="a", ylab= " ",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=1,ylim = c(ymin_, ymax_),col="white")
  #vect_bs = c(0.05,0.15,0.25,0.35,0.45,0.55);
  cc=0;
  for (bb in vect_bs) {
    cc=cc+1;
    points(as1,V_YO/vect_mses1[[as.character(cc)]],type='l',xlab="a", ylab= " ",
           cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=1,col="grey")
  }
  points(as1,V_YO/mses1,type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=1,col="black")
  points(as2,V_YO/mses3,type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=4)
  points(as1,V_YO/mses1_k01,type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=3)
  points(as1,V_YO/mses5,type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=6,col="magenta")
  abline(v=-0.5,col="red");
  abline(h = V_YO/MSE_add12,col="blue",lty=8);
  #abline(h = V_YO/MSE_extolk,col="magenta",lty=6);
  if (data_=="D4")
    legend("topright",legend=c(TeX('$\\bar{y}_1$'),TeX('$\\bar{y}_2$'),
                               TeX('$\\bar{y}_5$'),TeX('$\\bar{y}_6$'),
                               TeX('$\\bar{y}_7$'),TeX('$\\bar{y}_9$')), 
           col=c("blue","magenta","black","black","black","grey"),lty=c(8,6,1,4,3,1), 
           ncol=1, cex=1.8);
  dev.off();
  
  #max(30,
  #min(0,)
  ymin_=min(mses1,mses3,mses1_k01,MSE_add12,MSE_extolk,mses5)/V_YO;
  ymax_=max(mses1,mses3,mses1_k01,MSE_add12,MSE_extolk,mses5)/V_YO;
  namefile_=paste(directory_,"PSE_ALL_",data_,".pdf",sep="");
  pdf(namefile_);
  par(lwd=3)
  plot(as1,mses1/V_YO,type='l',xlab="a", ylab= " ",cex.lab=2, cex.axis=2, cex.main=2,
       cex.sub=2,lty=1,ylim = c(ymin_, ymax_))
  points(as2,mses3/V_YO,type='l',xlab="a", ylab= " ",cex.lab=2, cex.axis=2,
         cex.main=2, cex.sub=2,lty=2)
  points(as1,mses1_k01/V_YO,type='l',xlab="a", ylab= " ",cex.lab=2, cex.axis=2,
         cex.main=2, cex.sub=2,lty=3)
  points(as1,mses5/V_YO,type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=6,col="magenta")
  abline(v=-0.5,col="red");
  abline(h = MSE_add12/V_YO,col="blue",lty=8);
  #abline(h = MSE_extolk/V_YO,col="magenta",lty=6);
  dev.off()
  
  ## ----------------------------------------------------------------
  
  #   
  #   # f = 1 + a d_x + b d_x^2
  #   as=as00;
  #   #as = seq(-1.2,-0.05,length.out = 5000);
  #   #as = seq(-0.7,-0.05,length.out = 300);
  #   bs = seq(0.05,0.55,length.out = 300);
  #   length_as=length(as);
  #   length_bs=length(bs);
  #   mses_2d=matrix(0,nrow=length_as,ncol=length_bs);
  #   bias_2d=matrix(0,nrow=length_as,ncol=length_bs);
  #   t=u=0;
  #   for (a in as) {
  #     t=t+1;
  #     u=0;
  #     for (b in bs) {
  #       u=u+1;
  #       resu_ = Y_RPR03(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
  #       mses_2d[t,u] = resu_$MSE;
  #       bias_2d[t,u] = resu_$BIAS;
  #     }
  #   }
  #   cdivmse2d = c(V_YO/mses_2d);
  #   breaks <- seq(from=min(cdivmse2d), to=max(cdivmse2d),length.out = 10)
  #   png(paste(directory_,"PRE3D_",data_,".png",sep=""))
  #   persp3D(V_YO/mses_2d, x = as, y = bs, col = jet.col(length(breaks)-1),
  #           main = " ", clab = " ", breaks = breaks, xlab="a", ylab="b", zlab = " ", 
  #           ticktype = "detailed", theta=25)
  #   dev.off()
  
}


MM1 <- resultats[,1:4];
rownames(MM1) <- paste("D",1:4,sep="")
xresultats1 = xtable(MM1,digits=c(8));
print(names_res[1:4]);
print(xresultats1, floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=TRUE, include.colnames=FALSE)




##############################################################










