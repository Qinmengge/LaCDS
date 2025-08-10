#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
#include <string>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;

// all the input must be sorted according to t, large->small.

// [[Rcpp::export]]
Eigen::MatrixXd pAUC(Eigen::MatrixXd Truth,Eigen::MatrixXd Pred, int n){
  Eigen::MatrixXd rst((n-1)*n/2,2);
  rst.setZero((n-1)*n/2,2);
  int k=0;
  for (int i = 0 ; i < n-1 ; i++){
    for (int j = i+1 ; j < n ; j++){
      if(Truth(i,0)==Truth(j,0)){
        rst(k,0)=1;
      }
      else{
        rst(k,0)=0;
      }
      rst(k,1)=Pred(i,0)*Pred(j,0)+(1-Pred(i,0))*(1-Pred(j,0));
      k+=1;
      Rcout << k;
      Rcout << "\n";
    }
  }
  return(rst);
}

// [[Rcpp::export]]
double Max(double a, double b){
  if(a>b){
    return(a);
  }
  else{
    return(b);
  }
}



// [[Rcpp::export]]
List Update(Eigen::MatrixXd wker,Eigen::MatrixXd gtau,Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd ker_t, Eigen::MatrixXd beta_v, int max_t, int c, int r,int nker){
  double loglik=0;
  double xbeta=0;
  int tmp=0;
  Eigen::MatrixXd score_t(nker,1);
  score_t.setZero(nker,1);
  Eigen::MatrixXd score_v(c,1);
  score_v.setZero(c,1);
  Eigen::MatrixXd info_t(nker,nker);
  info_t.setZero(nker,nker);
  Eigen::MatrixXd info_v(c,c);
  info_v.setZero(c,c);
  Eigen::MatrixXd info_tv(nker,c);
  info_tv.setZero(nker,c);

  Eigen::MatrixXd Lambda(max_t,r);
  Lambda.setZero(max_t,r);

  Eigen::MatrixXd LambdaP(max_t,r);
  LambdaP.setZero(max_t,r);

  Eigen::MatrixXd Ind(max_t,r);
  Ind.setZero(max_t,r);

  Eigen::MatrixXd beta_t(max_t,1);
  beta_t=wker*ker_t;


  for (int i = 0 ; i < r ; i++){
    xbeta=(X.row(i)*beta_v)(0,0);
    for (int s = 1 ; s <= max_t ; s++){
      if (t(i,0) >= s){
        Lambda(s-1,i)=1/(1+exp(-beta_t(s-1,0)-xbeta));
        LambdaP(s-1,i)=Lambda(s-1,i)*(1-Lambda(s-1,i));
        if (t(i,0) == s and ind(i,0) == 1){
          loglik=loglik+gtau(i,0)*log(Lambda(s-1,i));
        }
        else {
          loglik=loglik+gtau(i,0)*log(1-Lambda(s-1,i));
        }
      }
      else{
        break;
      }
    }
  }

  for (int i = 0 ; i < r ; i++){
    tmp=t(i,0)-1;
    Ind(tmp,i)=ind(i,0);
  }

  score_t=wker.transpose()*(Ind-Lambda)*gtau;
  score_v=(((Ind-Lambda)*(gtau.asDiagonal())*X).colwise().sum()).transpose();
  info_t=(wker.transpose())*((LambdaP*gtau).asDiagonal())*wker;
  info_tv=wker.transpose()*LambdaP*(gtau.asDiagonal())*X;
  info_v=(X.transpose())*(((LambdaP*(gtau.asDiagonal())).colwise().sum()).asDiagonal())*X;

  List result;
  result["loglik"]=loglik;
  result["score_t"]=score_t;
  result["score_v"]=score_v;
  result["info_t"]=info_t;
  result["info_tv"]=info_tv;
  result["info_v"]=info_v;
  return result;
}

// [[Rcpp::export]]
double Like(Eigen::MatrixXd wker,Eigen::MatrixXd gtau,Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd ker_t, Eigen::MatrixXd beta_v, int max_t, int c, int r){
  double loglik=0;
  double xbeta=0;

  Eigen::MatrixXd Lambda(max_t,r);
  Lambda.setZero(max_t,r);

  Eigen::MatrixXd beta_t(max_t,1);
  beta_t=wker*ker_t;

  for (int i = 0 ; i < r ; i++){
    xbeta=(X.row(i)*beta_v)(0,0);
    for (int s = 1 ; s <= max_t ; s++){
      if (t(i,0) >= s){
        Lambda(s-1,i)=1/(1+exp(-beta_t(s-1,0)-xbeta));
        if (t(i,0) == s and ind(i,0) == 1){
          loglik=loglik+gtau(i,0)*log(Lambda(s-1,i));
        }
        else {
          loglik=loglik+gtau(i,0)*log(1-Lambda(s-1,i));
        }
      }
      else{
        break;
      }
    }
  }

  return loglik;
}



// [[Rcpp::export]]
List NW(Eigen::MatrixXd wker,Eigen::MatrixXd gtau, Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd ker_t, Eigen::MatrixXd beta_v,double tol, int max_iter, double alpha){
  int r = X.rows();
  int c = X.cols();
  int max_t=t.maxCoeff();
  int nker=wker.cols();
  double step=1;
  Eigen::MatrixXd L(c,c);
  L.setOnes(c,c);
  Eigen::MatrixXd L1(c,c);
  L1.setOnes(c,c);
  Eigen::MatrixXd Lab(c,nker);
  Lab.setOnes(c,nker);
  Eigen::MatrixXd Ls(c,1);
  Ls.setOnes(c,1);
  Eigen::MatrixXd Ls1(nker,1);
  Ls1.setOnes(nker,1);
  Ls1=Ls1/10000000000;
  Eigen::MatrixXd Ls2(c,1);
  Ls2.setOnes(c,1);
  Eigen::MatrixXd Ms(c,c);
  Ms=Ls2.asDiagonal();
  Ms=Ms/10000000000;// don't add too larger.
  Eigen::MatrixXd ker2_t=ker_t;
  Eigen::MatrixXd beta2_v=beta_v;
  Eigen::MatrixXd dirc_t=ker_t;
  Eigen::MatrixXd dirc_v=beta_v;
  Eigen::MatrixXd A_inv(nker,nker);
  A_inv.setZero(nker,nker);
  Eigen::MatrixXd A(nker,nker);
  A.setZero(nker,nker);
  Eigen::MatrixXd B(nker,c);
  B.setZero(nker,c);
  Eigen::MatrixXd C(c,c);
  C.setZero(c,c);
  Eigen::MatrixXd schur(c,c);
  schur.setZero(c,c);
  Eigen::MatrixXd score_t(nker,1);
  score_t.setZero(nker,1);
  Eigen::MatrixXd score_v(c,1);
  score_v.setZero(c,1);
  double ll,ll2;
  double diff_v=0;
  List update=Update(wker,gtau,t,X,ind,ker_t,beta_v,max_t,c,r,nker);
  ll=update["loglik"];
  ll2=ll;
  for (int i = 0 ; i<=max_iter ; i++){
    step=1;
    A=update["info_t"];
    B=update["info_tv"];
    C=update["info_v"];
    score_t=update["score_t"];
    score_v=update["score_v"];
    A_inv=A.inverse();
    const Eigen::LLT<MatrixXd> llt(C-(B.transpose())*(A_inv)*(B));
    L=MatrixXd(llt.matrixL());
    L1=L.inverse();
    Lab=L1*(B.transpose())*A_inv;
    Ls=L1*score_v;
    dirc_t=((A_inv+(Lab.transpose())*Lab)*score_t-(Lab.transpose())*Ls);
    /* alpha backtracking, 1/2, */
    dirc_v=((L1.transpose())*Ls-(L1.transpose())*Lab*score_t);

    while (Like(wker,gtau,t,X,ind,ker_t+step*dirc_t,beta_v+step*dirc_v,max_t,c,r)<(Like(wker,gtau,t,X,ind,ker_t,beta_v,max_t,c,r)+0.5*step*((score_t.transpose()*dirc_t)(0,0)+(score_v.transpose()*dirc_v)(0,0)))){
      step*=alpha;
      Rcout<<"Step size control activated\n";
    }

    ker2_t=ker_t+step*dirc_t;
    beta2_v=beta_v+step*dirc_v;


    diff_v=(beta2_v-beta_v).norm();
    update=Update(wker,gtau,t,X,ind,ker2_t,beta2_v,max_t,c,r,nker);
    ll2=update["loglik"];
    if((ll2-ll)*(ll2-ll)<tol || diff_v<tol){
      List result;
      result["ker_t"]=ker2_t;
      result["beta_v"]=beta2_v;
      return result;
    }
    ll=ll2;
    ker_t=ker2_t;
    beta_v=beta2_v;
  }
  List result;
  result["ker_t"]=ker2_t;
  result["beta_v"]=beta2_v;
  return result;
}


// [[Rcpp::export]]
List lik(Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd Beta_t, Eigen::MatrixXd Beta_v, Eigen::MatrixXd pi,int max_t, int r, int c, int group){
  Eigen::MatrixXd tau(r,group);
  tau.setZero(r,group);
  double deno=0;
  double XBv=0;
  double pif=0;
  double lambda=0;
  double likeli2=0;
  for(int j=0 ; j<r ; j++){
    for(int g=0 ; g<group ; g++){
      XBv=((X.block(j,0,1,c))*(Beta_v.block(0,g,c,1)))(0,0);
      pif=pi(g,0);
      for(int s=0 ; s<t(j,0)-1 ; s++){
        lambda=1/(1+exp(-Beta_t(s,g)-XBv));
        pif=pif*(1-lambda);
      }
      int s=t(j,0)-1;
      lambda=1/(1+exp(-Beta_t(s,g)-XBv));
      if(ind(j,0)==1){
        pif=pif*(lambda);
      }
      else{
        pif=pif*(1-lambda);
      }
      tau(j,g)=pif;
      deno=deno+pif;
    }
    likeli2=likeli2+log(deno);
    tau.block(j,0,1,group)=tau.block(j,0,1,group)/deno;

    deno=0;
  }
  pi=tau.colwise().sum().transpose()/r;

  List result;
  result["likeli2"]=likeli2;
  result["tau"]=tau;
  result["pi"]=pi;

  return(result);
}

/* complete data likelihood appr*/

// [[Rcpp::export]]
List Inf2(Eigen::MatrixXd wker,Eigen::MatrixXd pi, Eigen::MatrixXd tau,Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd Ker_t, Eigen::MatrixXd Beta_v, int max_t, int c, int r, int group,int nker){
  double xbeta=0;
  double score_pi=0;
  int tmp=0;
  Eigen::MatrixXd Beta_t(max_t,group);
  Beta_t=wker*Ker_t;
  Eigen::MatrixXd Score(r,group*(1+c+nker)-1);
  Score.setZero(r,group*(1+c+nker)-1);
  Eigen::MatrixXd Info(group*(1+c+nker)-1,group*(1+c+nker)-1);
  Info.setZero(group*(1+c+nker)-1,group*(1+c+nker)-1);

  Eigen::MatrixXd gtau(r,1);
  gtau.setZero(r,1);

  Eigen::MatrixXd ker_t(nker,1);
  ker_t.setZero(nker,1);
  Eigen::MatrixXd beta_t(max_t,1);
  beta_t.setZero(max_t,1);
  Eigen::MatrixXd beta_v(c,1);
  beta_v.setZero(c,1);

  Eigen::MatrixXd Lambda(max_t,1);
  Lambda.setZero(max_t,1);

  Eigen::MatrixXd Ind(max_t,1);
  Ind.setZero(max_t,1);

  Eigen::MatrixXd ones(r,1);
  ones.setOnes(r,1);

  for (int i = 0 ; i < r ; i++){
    for (int g=0; g<group; g++){
      score_pi=0;
      gtau=tau.col(g);
      ker_t=Ker_t.col(g);
      beta_t=Beta_t.col(g);
      beta_v=Beta_v.col(g);
      Lambda.setZero(max_t,1);
      Ind.setZero(max_t,1);
      tmp=t(i,0)-1;
      Ind(tmp,0)=ind(i,0);
      if (g<group-1){
        Score(i,g)=gtau(i,0)/pi(g,0)-tau(i,group-1)/pi(group-1,0);
      }
      xbeta=(X.row(i)*beta_v)(0,0);
      for (int s = 1 ; s <= max_t ; s++){
        if (t(i,0) >= s){
          Lambda(s-1,0)=1/(1+exp(-beta_t(s-1,0)-xbeta));

        }
      }
      Score.block(i,group-1+g*(c+nker),1,nker)=((wker.transpose())*(Ind-Lambda)*gtau(i,0)).transpose();
      Score.block(i,group-1+g*(c+nker)+nker,1,c)=((Ind-Lambda).colwise().sum()(0,0)*gtau(i,0)*X.row(i));
    }
  }

  Info=Score.transpose()*Score;


  List result;
  result["Score"]=Score;
  result["Info"]=Info;


  return result;
}




// [[Rcpp::export]]
List EM(Eigen::MatrixXd wker,Eigen::MatrixXd Z,int group,Eigen::MatrixXd t, Eigen::MatrixXd X, Eigen::MatrixXd ind, Eigen::MatrixXd beta_t, Eigen::MatrixXd beta_v,double tol, double tol2, int max_iter, int min_iter, double alpha){

  int iter=0;
  bool NAN_ind=false;

  double DAEM=1;
  double DAEM_rate=1.3;
  int r = X.rows();
  int c = X.cols();
  int max_t=t.maxCoeff();
  int nker=wker.cols();

  Eigen::MatrixXd pi(group,1);
  pi.setOnes(group,1);
  pi=pi/group;

  /*Eigen::MatrixXd Z(r,group);
  Z.setZero(r,group);
  Z.block(0,0,r/2,1)=MatrixXd::Ones(r/2,1);
  Z.block(r/2,1,r/2,1)=MatrixXd::Ones(r/2,1);*/

  Eigen::MatrixXd Z2(r,group);
  Z2.setZero(r,group);

  Eigen::MatrixXd rtau(r,group);
  rtau.setZero(r,group);
  rtau.col(0)=MatrixXd::Ones(r,1);

  /*Eigen::MatrixXd tau(r,group);
  tau.setZero(r,group);
  tau.col(0)=MatrixXd::Ones(r,1);*/

  Eigen::MatrixXd tau(r,group);
  tau=Z;

  Eigen::MatrixXd eg(r,1);
  eg.setZero(r,1);

  Eigen::MatrixXd Ker_t(nker,group);
  Ker_t.setZero(nker,group);

  Eigen::MatrixXd Beta_t(max_t,group);
  Beta_t=wker*Ker_t;

  Eigen::MatrixXd Beta_v(c,group);
  Beta_v.setZero(c,group);

  Eigen::MatrixXd Ker_t2(nker,group);
  Ker_t2.setZero(nker,group);

  Eigen::MatrixXd Beta_v2(c,group);
  Beta_v2.setZero(c,group);

  Eigen::MatrixXd likelihood(1,max_iter);
  likelihood.setZero(1,max_iter);

  double XBv=0;
  double lambda=0;
  double pif=0;
  double deno=0;
  double likeli=0;
  double likeli2=0;

  for (int i = 0 ; i<=max_iter ; i++){
    iter+=1;
    for (int g = 0; g<group ; g++){
      Eigen::MatrixXd gtau=tau.col(g);
      //Eigen::MatrixXd gtau=Z.col(g);
      if(gtau==eg){
        continue;
      }
      List BETA;
      BETA=NW(wker,gtau, t, X, ind, Ker_t.col(g), Beta_v.col(g), tol, max_iter, alpha);
      Eigen::MatrixXd ker_t_temp=BETA["ker_t"];
      Eigen::MatrixXd beta_v_temp=BETA["beta_v"];
      Ker_t.col(g)=ker_t_temp;
      Beta_v.col(g)=beta_v_temp;
    }

    //Rcout << "beta_v*****\n";
    //Rcout << Beta_v;
    //Rcout << "*****\n";

    // update pi

    Beta_t=wker*Ker_t;

    List likl=lik(t, X, ind, Beta_t, Beta_v, pi, max_t, r, c, group);

    pi=likl["pi"];
    likeli2=likl["likeli2"];
    tau=likl["tau"];
    /*tau=tau.array().pow(DAEM);*/
    for(int j=0 ; j<r ; j++){
      tau.block(j,0,1,group)=tau.block(j,0,1,group)/tau.block(j,0,1,group).rowwise().sum()(0,0);
    }

    if (DAEM>=1){
      DAEM=1;
    }
    else{
      DAEM=DAEM*DAEM_rate;
    }

    double diff_v=(Beta_v-Beta_v2).norm();
    double diff_t=(Ker_t-Ker_t2).norm();
    double diff_lik=abs(likeli2-likeli)/abs(likelihood(0,0)-likeli2);

    //Rcout << "\nv\n";
    //Rcout << diff_v;
    //Rcout << "\nt\n";
    //Rcout << diff_t;
    Rcout << "\ndiff_lik\n";
    Rcout << diff_lik;
    Rcout << "\ndiff_v\n";
    Rcout << diff_v;

    NAN_ind=Rcpp::NumericVector::is_na(diff_lik);
    Rcout << "\nNAN_ind\n";
    Rcout << NAN_ind;
    if(NAN_ind==1){
      List result;
      result["NAN_ind"]=1;
      return(result);
    }

    //if(iter>50 && (diff_v<tol2 || (diff_lik)<tol2))
    if(iter>min_iter && (diff_v<tol2)){
      break;
     List result;
     result["likelihood"]=likelihood;
     result["Z"]=tau;
     result["Ker_t"]=Ker_t;
     Beta_t=wker*Ker_t;
     result["Beta_t"]=Beta_t;
     result["Beta_v"]=Beta_v;
     //Rcout<<"\ninf_start\n";
     List Inference=Inf2(wker,pi, tau,t, X, ind, Ker_t, Beta_v, max_t, c, r, group,nker);
     //Rcout<<"\ninf_done\n";
     Eigen::MatrixXd Info=Inference["Info"];
     Eigen::MatrixXd Score=Inference["Score"];
     result["Info"]=Info;
     result["Score"]=Score;
     result["NAN_ind"]=0;
     return(result);
    }


    likelihood(0,i)=likeli2;

    Ker_t2=Ker_t;
    Beta_v2=Beta_v;
    likeli=likeli2;



  }
  List result;
  result["likelihood"]=likelihood;
  result["Z"]=tau;
  result["Ker_t"]=Ker_t;
  Beta_t=wker*Ker_t;
  result["Beta_t"]=Beta_t;
  result["Beta_v"]=Beta_v;
  //Rcout<<"\ninf_start\n";
  List Inference=Inf2(wker,pi, tau,t, X, ind, Ker_t, Beta_v, max_t, c, r, group,nker);
  //Rcout<<"\ninf_done\n";
  Eigen::MatrixXd Info=Inference["Info"];
  Eigen::MatrixXd Score=Inference["Score"];
  result["Info"]=Info;
  result["Score"]=Score;
  result["NAN_ind"]=0;
  return(result);
}

