#ifndef WEIGHT_FUNCTIONS_H
#define WEIGHT_FUNCTIONS_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>


class Weight_Functions{
 
 private:
   
   double ialpha = 10.0;
   double Router = 1.2;
   double Rinner = 0.8;
   double Rrange = 0.4;
   
 public:
  
  Weight_Functions(){
    Rinner = 1.2;
    Router = 0.8;
    Rrange = Router-Rinner;
  };
  
  ~Weight_Functions(){}
  
  void setting(double Rin, double Rout){
    Rinner = Rin;
    Router = Rout;
    Rrange = Rout - Rin;
  }
  
  void pair(double r, double &W, double &dW, bool flag = false){
    if (r > Router){
        W = 1.0;
       dW = 0.0;
    }else if (r < Rinner){
       W = 0.0;
      dW = 0.0;
    }else {
      double u  = (r - Rinner)/Rrange;
      double u2 = u*u;
      double u3 = u2*u;
       W = u3*( 6.0*u2 - 15.0*u + 10.0);
      dW = u2*(     u2 -  2.0*u +  1.0)*30.0/Rrange;
    }
    
    if(flag){
      W = 1-W; 
     dW = -dW;
    }
    
  }
  
  void select_auto_min(int num, double *x, double *rsq, double *f, double &E, double Einf=0.0) {
    
    double *dRd_r = new double[num];
    double *dRn_r = new double[num];
    
    double sum_Rd = 0.0;
    double sum_Rn = 0.0;
    for(int j = 0; j < num; j++){
      
      double r = sqrt(rsq[j]);
      double Rn   = exp(-r*ialpha);
      double Rd   = r*Rn;
      sum_Rn += Rn;
      sum_Rd += Rd;
      
      double inv_r = 1.0/r;
      dRn_r[j] = -ialpha*Rn*inv_r;
      dRd_r[j] = (1.0 - r*ialpha)*Rn*inv_r;
      
    }
    
    double U, W,dW;
    U = sum_Rd/sum_Rn;
    pair(U,W,dW);
    
    double Ea  = E - Einf;
    
    E = W*Ea + Einf;
    
    double dE_Rd = dW*(1.0/sum_Rn)*Ea;
    double dE_Rn = dW*(sum_Rd)*(-1.0/(sum_Rn*sum_Rn))*Ea;
    
    for (int j = 0; j < num; j++){
      double *xj = &x[j*3];
      double *fj = &f[j*3];
      double fcut = dE_Rd*dRd_r[j] + dE_Rn*dRn_r[j];
      for (int k = 0; k < 3; k++)fj[k] = W*fj[k] - fcut*xj[k];
    }
    
    delete [] dRd_r;
    delete [] dRn_r;
  
  }

  void select_min(int num, double *x, double *rsq, double *f, double &E, double Einf=0.0) {
    
    double Rmin = 100.0;
    int jmin;
    
    for(int j = 0; j < num; j++){
      double r = sqrt(rsq[j]);
      if (Rmin > r){
        Rmin = r;
        jmin = j;
      }
    }
    
    double W,dW;
    pair(Rmin,W,dW);
    
    double Ea  = E - Einf;
    
    E = W*Ea + Einf;
    
    double *xj = &x[jmin*3];
    double *fj = &f[jmin*3];
    double fcut = -Ea*dW/Rmin;
    for (int k = 0; k < 3; k++)fj[k] = W*fj[k] + fcut*xj[k];
  
  }

  void select_max(int num, double *x, double *rsq, double *f, double &E, double Einf=0.0) {
    
    double Rmax = 0.0;
    int jmax;
    
    for(int j = 0; j < num; j++){
      double r = sqrt(rsq[j]);
      if (Rmax < r){
        Rmax = r;
        jmax = j;
      }
    }
    
    double W,dW;
    bool flag = true;
    pair(Rmax,W,dW,flag);
    
    double Ea  = E - Einf;
    
    E = W*Ea + Einf;
    
    double *xj = &x[jmax*3];
    double *fj = &f[jmax*3];
    double fcut = -Ea*dW;
    for (int k = 0; k < 3; k++)fj[k] = W*fj[k] + fcut*xj[k];
  
  }

};

#endif
