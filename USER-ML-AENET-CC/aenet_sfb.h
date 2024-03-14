#ifndef AENET_SFB_H
#define AENET_SFB_H

#include <cmath>
#include <vector>

class AENET_SFB{
  
  int aenet_ver = 0;
  int fctype_r=0;
  int fctype_a=0;
  double h_r = 3.0;
  double h_a = 1.25;
  
  double Rc_max = 0.0;
  double Rc_r;
  double Rc_a;
  double Rc_r_sq;
  double Rc_a_sq;
  
  int nGtot;
  int nGr;
  int nGa;
  
  bool multi=false;
  
  int nenv;
  int *env_map;
  int *spin;
  
  void calc_chebyshev_polynomial(double x, double x0, double x1, int n, double *T){
    
    int i; 	
    double XX = (2.0*x - x0 - x1)/(x1 - x0);
   
    T[0] = 1.0;
    if (n > 1) {
      T[1] = XX;
      for (i=2; i<n; i++ ) T[i] = 2.0*XX*T[i-1] - T[i-2];
    }
    
    return;
  
  }
  
  void diff_chebyshev_polynomial(double x, double x0, double x1, int n, double *dT){
    
    int i; 
    double U1,U2,U3;
    double XX = (2.0*x - x0 - x1)/(x1 - x0);
   
    dT[0] = 0.0;
    if (n > 1) {
      dT[1] = 1.0;
      U1 = 1.0;
      U2 = 2.0*XX;
      for (i=2; i<n; i++ ) {
        dT[i] = U2*i;
        U3 = 2.0*XX*U2 - U1;
        U1 = U2;
        U2 = U3;
      }
    }
    
    for (i=0; i<n; i++ ) {dT[i] = dT[i]*2.0/(x1 - x0);}
    
    return;
    
  }
  
  void compute_sfb_rad(double r, double Rc, int nG, double *G, double *dGdr){
    
    calc_chebyshev_polynomial(r, 0.0, Rc, nG,  G);
    diff_chebyshev_polynomial(r, 0.0, Rc, nG, dGdr);
    
  }
  
  void compute_sfb_ang(int version, double cos_ijk, int nG, double *G, double *dGdcos){
  
    static constexpr double PI = 3.14159265358979323846;
    double h_cos;
    
    if(version == -1){
      calc_chebyshev_polynomial(cos_ijk, 0.0, PI, nG,  G);
      diff_chebyshev_polynomial(cos_ijk, 0.0, PI, nG, dGdcos);
    }else if(version == 1){
      h_cos = 0.50*(cos_ijk+1.0);
      calc_chebyshev_polynomial(h_cos, -1.0, 1.0, nG,  G);
      diff_chebyshev_polynomial(h_cos, -1.0, 1.0, nG, dGdcos);
      for (int i =0; i < nG; i++){dGdcos[i] = 0.5*dGdcos[i];}
    }else{
      calc_chebyshev_polynomial(cos_ijk, -1.0, 1.0, nG,  G);
      diff_chebyshev_polynomial(cos_ijk, -1.0, 1.0, nG, dGdcos);
    }
    
  }
  
  void f_c(int fctype, double r, const double rc, double h, double &fc, double &dfc){
    
    static constexpr double PI = 3.14159265358979323846;
    double x,y,dy;
    
    if (r > rc){
       fc = 0.0;
      dfc = 0.0;
    }else{
      if (fctype == 1){
        x   = r/rc;
         fc = 1.0 + h*pow(x,h+1) - (h+1)*pow(x,h);
        dy  = h*(h+1)*pow(x,h) - (h+1)*h*pow(x,h-1);
        dfc = dy/rc;
      } else if (fctype == 2){
        x   = (r - rc)/h;
        // fc = (x*x)/(1 + x*x);
        //dy  = (2*x)/((1+x*x)*(1+x*x));
	double x2 = x*x;
         fc = x2/(1+x2);
        dy  = (2*x)/((1+x2)*(1+x2));
        dfc = dy/h;
      } else {
        x   =  PI*r/rc;
         fc =  0.5*(cos(x) + 1.0);
        dy  = -0.5*sin(x);
        dfc = dy*(PI/rc);
      }
    }
    
    
  }
  
  public:
  
  ~AENET_SFB(){
    delete [] env_map;
    delete [] spin;    
  }
  
  void set_version(int version){
    
    if(version == -1) {
       aenet_ver = -1;
       fctype_r  = 0;
       fctype_a  = 0;
    }else if (version == 0){
       aenet_ver = 0;
       fctype_r  = 0;
       fctype_a  = 0;
    }else if (version == 1){
       aenet_ver = 1;
       fctype_r  = 0;
       fctype_a  = 0;
    }else if (version == 2){
       aenet_ver = 1;
       fctype_r  = 1;
       fctype_a  = 1;
       h_r = 3.0;
       h_a = 5.0;
    }else if (version == 3){
       aenet_ver = 1;
       fctype_r  = 1;
       fctype_a  = 2;
       h_r = 3.0;
       h_a = 1.25;
    }else if (version == 4){
       aenet_ver = 1;
       fctype_r  = 1;
       fctype_a  = 2;
       h_r = 3.0;
       h_a = 1.25*(Rc_a/4.5);
    }else {
       aenet_ver = 0;
       fctype_r  = 0;
       fctype_a  = 0;
    }
   
  }
  
  void set_nG(int nt, int nr, int na){nGtot=nt;nGr=nr;nGa=na;}
  
  void set_multi(){multi = true;}
  
  void set_Rc(double Rmax, double Rr, double Ra){
    
    Rc_max = Rmax;
    Rc_r = Rr;
    Rc_a = Ra;
    Rc_r_sq = Rc_r*Rc_r;
    Rc_a_sq = Rc_a*Rc_a;
    
  }
  
  void set_envtype(int n, int nelements, std::vector<int> lmap){
    
    nenv = n;
    
    spin = new int [nenv]();
    env_map = new int [nelements]();
    
    if(multi){
      int s = -nenv/2;
      for (int i = 0; i < nenv; i++){
        if (s == 0 && nenv%2 == 0){s++;}
        spin[i] = s;
        s++;
      }
      
      for (int i = 0; i < nelements; i++) {env_map[i] = lmap[i];}
    }
  }
  
  int get_nG(){return nGtot;}
  
  double get_Rc_max(){return Rc_max;}
  
  void compute(int jnum, double *x, double *rsq, int *type, double *Gall, double *dGdx){
  
    double *Gr     = new double [nGr]{};
    double *Ga     = new double [nGa]{};
    double *dGdr   = new double [nGr]{};
    double *dGdcos = new double [nGa]{};
    double *dGrdxj = new double [nGr*3]{};
    double *dGadxj = new double [nGa*3]{};
    double *dGadxk = new double [nGa*3]{};
   

    //j-loop start
    for (int j = 0; j < jnum; j++){
      
      double xj[3], dr_dxj[3];

      double *x_j = &x[3*j];
      for(int i = 0; i < 3; i++)xj[i] = x_j[i];
      double rsqj = rsq[j];
      double rj = sqrt(rsqj);
      double rinvj = 1.0/rj;
      for(int i = 0; i < 3; i++) dr_dxj[i] = xj[i]*rinvj;
      
      if (rsqj <= Rc_r_sq){
        
        compute_sfb_rad(rj, Rc_r, nGr, Gr, dGdr);
        double fcj,dfcj;
        f_c(fctype_r, rj, Rc_r, h_r, fcj, dfcj);
        
        for (int n = 0; n < nGr; n++){
          double dGdfc = fcj*dGdr[n] + dfcj*Gr[n];
          double *dGrdx_j = &dGrdxj[n*3];
          for(int i = 0; i < 3; i++) dGrdx_j[i] = dGdfc*dr_dxj[i];
        }
        
        int nstart = 0;
        
        double *Gj = &Gall[nstart];
        for (int n = 0; n < nGr; n++) Gj[n] += fcj*Gr[n];
        
        double *dGdxj = &dGdx[(j*nGtot + nstart)*3];
        for(int ni = 0; ni < nGr*3; ni++)dGdxj[ni] += dGrdxj[ni];
        
        if (multi) if(env_map[type[j]] >= 0){
                  
          int jspin = spin[env_map[type[j]]];
          double sfcj = jspin*fcj;
          
          int nstart = nGr+nGa;
          
          double *Gj = &Gall[nstart];
          for (int n = 0; n < nGr; n++) Gj[n] += sfcj*Gr[n];
          
          double *dGdxj = &dGdx[(j*nGtot + nstart)*3];
          for(int ni = 0; ni < nGr*3; ni++)dGdxj[ni] += jspin*dGrdxj[ni];
          
        }
        
      }

      if (rsqj > Rc_a_sq)continue;
      
      //k-loop start 
      for (int k = j+1; k < jnum; k++){
      
        double rsqk = rsq[k];
        if (rsqk > Rc_a_sq)continue;
        
        double xk[3], dr_dxk[3], dcos_dxj[3], dcos_dxk[3];

        double *x_k = &x[3*k];
        for(int i = 0; i < 3; i++) xk[i] = x_k[i];
        
        double rk = sqrt(rsqk);
        double rinvjk = 1.0/(rj*rk);
        double cos_jk = (xj[0]*xk[0] + xj[1]*xk[1] + xj[2]*xk[2]) * rinvjk;
        
        compute_sfb_ang(aenet_ver,cos_jk, nGa, Ga, dGdcos);
        double fcj,dfcj,fck,dfck;
        f_c(fctype_a, rj, Rc_a, h_a, fcj, dfcj);
        f_c(fctype_a, rk, Rc_a, h_a, fck, dfck);
                
        double fcjk = fcj*fck;
        
        double rinvk = 1.0/rk;
        for(int i = 0; i < 3; i++) dr_dxk[i] = xk[i]*rinvk;
        
        double rinvsqj = 1.0/rsqj;
        double rinvsqk = 1.0/rsqk;
        for(int i = 0; i < 3; i++) dcos_dxj[i] = -cos_jk*xj[i]*rinvsqj + xk[i]*rinvjk;
        for(int i = 0; i < 3; i++) dcos_dxk[i] = -cos_jk*xk[i]*rinvsqk + xj[i]*rinvjk;
        
        for (int n = 0; n < nGa; n++){
          double dG_fc  =  fcj*fck*dGdcos[n];
          double G_dfcj = dfcj*fck*Ga[n];
          double G_dfck = fcj*dfck*Ga[n];
          double *dGadx_j = &dGadxj[n*3];
          double *dGadx_k = &dGadxk[n*3];
          for(int i = 0; i < 3; i++) dGadx_j[i] = dG_fc*dcos_dxj[i] + G_dfcj*dr_dxj[i];
          for(int i = 0; i < 3; i++) dGadx_k[i] = dG_fc*dcos_dxk[i] + G_dfck*dr_dxk[i];
        }
        
        int nstart = nGr;
        
        double *Gjk = &Gall[nstart];
        for (int n = 0; n < nGa; n++) Gjk[n] += fcjk*Ga[n];
        
        double *dGdxj = &dGdx[(j*nGtot + nstart)*3];
        for (int ni = 0; ni < nGa*3; ni++)dGdxj[ni] += dGadxj[ni];
        
        double *dGdxk = &dGdx[(k*nGtot + nstart)*3];
        for (int ni = 0; ni < nGa*3; ni++)dGdxk[ni] += dGadxk[ni];
        
        if (multi) if(env_map[type[j]] >= 0 && env_map[type[k]] >= 0){
                    
          int jkspin = spin[env_map[type[j]]]*spin[env_map[type[k]]];          
          double sfcjk = jkspin*fcj*fck;
          
          int nstart = nGr+nGa+nGr;
          
          double *Gjk = &Gall[nstart];
          for (int n = 0; n < nGa; n++) Gjk[n] += sfcjk*Ga[n];
          
          double *dGdxj = &dGdx[(j*nGtot + nstart)*3];
          for (int ni = 0; ni < nGa*3; ni++)dGdxj[ni] += jkspin*dGadxj[ni];
          
          double *dGdxk = &dGdx[(k*nGtot + nstart)*3];
          for (int ni = 0; ni < nGa*3; ni++)dGdxk[ni] += jkspin*dGadxk[ni];
          
        }
        
      }//k-loop end
            
    }//j-loop end
    
    
    delete [] Gr;
    delete [] Ga;
    delete [] dGdr;
    delete [] dGdcos;
    delete [] dGrdxj;
    delete [] dGadxj;
    delete [] dGadxk;
  
  }


};

#endif

