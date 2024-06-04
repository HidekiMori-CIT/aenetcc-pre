#ifndef AENET_POT_H
#define AENET_POT_H

#include <string> 
#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include "ann.h"
#include "aenet_sfb.h"


template <typename T>
std::vector<T> read_Fbinary(std::ifstream &fin){
    int size;
    fin.read( ( char * ) &size, sizeof(int));
    std::vector<T> data(size/sizeof(T));
    fin.read( ( char * ) &data[0], size ); 
    fin.read( ( char * ) &size, sizeof(int));
    return data;
}


class AENET_POT {
  
  std::string my_atomtype;
  
  double Escale;
  double Eshift;
  double Eatom;
  double Einf;
  
  ANN net;
  AENET_SFB sfb;
  
  bool skip_flag = false;

  public:
  
  ~AENET_POT(){}
  
  void skip_flag_on(){skip_flag = true;}
  bool skip(){return skip_flag;}

  void  set_pot(
    int aenet_version,
    std::vector<char> atomtype,
    std::vector<int> nnodes,std::vector<int> f_a, std::vector<double> W,
    int nelements, char** elements,
    int nenv, std::vector<char> envtypes,
    std::vector<double> sfparam,
    std::vector<double> sfval_avg, std::vector<double> sfval_cov,
    double scale, double shift, std::vector<double> Ea,
    int nTypes,std::vector<char> typeName){
    
    my_atomtype = std::string(atomtype.begin(),atomtype.end());
    
    std::vector<double> Wb(W.size());
    
    int eshift = 0;
    for (int l = 0; l < nnodes.size()-1; l++){
      for (int i = 0; i < nnodes[l+1]; i++){
        for (int j = 0; j < nnodes[l]; j++){
          Wb[i*(nnodes[l]+1) + 1 + j + eshift] = W[j*nnodes[l+1] + i + eshift];
        } 
        Wb[i*(nnodes[l]+1) + eshift] = W[nnodes[l]*nnodes[l+1] + i + eshift];
      }
      eshift += nnodes[l]*nnodes[l+1] + nnodes[l+1];
    }
    
    net.set_ann(nnodes, f_a, Wb);
    
    
    if(nenv > 1)sfb.set_multi();
    
    std::string etype;
    int elen = envtypes.size()/nenv;
    std::vector<int> lmap(nelements);
    for (int i = 0; i < nelements; i++) {
      lmap[i] = -1;
      for (int j = 0; j < nenv; j++){
        etype = std::string(envtypes.begin()+elen*j,envtypes.begin()+elen*(j+1));
        etype = etype.erase(etype.find_last_not_of(" ") + 1);
        etype = etype.erase(0,etype.find_first_not_of(" "));
        if (strcmp(etype.c_str(),elements[i]) == 0){
          lmap[i] = j;
          break;
        }
      }
    }
    
    sfb.set_envtype(nenv, nelements, lmap);
    
    
    int nGr   = sfparam[1]+1;
    int nGa   = sfparam[3]+1;
    int nGtot = nGr + nGa;
    if (nenv > 1) {nGtot = nGtot + (nGr + nGa);}
    
    double Rc_r  = sfparam[0];
    double Rc_a  = sfparam[2];
    double Rc_max = (Rc_r > Rc_a) ? Rc_r : Rc_a;
    
    sfb.set_Rc(Rc_max,Rc_r,Rc_a);
    sfb.set_nG(nGtot, nGr, nGa);
    
    double *Gscale = new double [nGtot];
    double *Gshift = new double [nGtot];
    
    for (int i = 0; i < nGtot; i++){
      Gshift[i] = sfval_avg[i];
      Gscale[i] = 1.0/sqrt(sfval_cov[i] - Gshift[i]*Gshift[i]);
    }
    
    net.set_scale(nGtot, Gscale, Gshift);
    
    delete [] Gscale;
    delete [] Gshift;


    Escale = 1.0/scale;
    Eshift = shift;
    
    std::string ithtype;
    int tlen = typeName.size()/nTypes;
    for (int i = 0; i < nTypes; i++){
      ithtype = std::string(typeName.begin()+tlen*i,typeName.begin()+tlen*(i+1));
      if (my_atomtype == ithtype){
        Eatom = Ea[i];
        break;
      }
    }
    
    
    double *Gall  = new double [nGtot]{};
    double *dEdG  = new double [nGtot]{};
    
    double E;
    
    net.compute(Gall, E, dEdG);
    
    Einf = Escale*E + Eshift + Eatom;
    
    delete [] Gall;
    delete [] dEdG;
    
    
    sfb.set_version(aenet_version);
    
  }
  
  double get_Rc_max(){return sfb.get_Rc_max();}
  
  double get_Einf(){return Einf;}
  
  void compute(int jnum, double *x, double *rsq, int *type, double &E, double *f){
    
    double E0;
    
    int nGtot = sfb.get_nG();
    
    double *Gall  = new double [nGtot]{};
    double *dEdG  = new double [nGtot]{};
    double *dGdx  = new double [jnum*nGtot*3]{};
    
    //compute descriptor
    sfb.compute(jnum, x, rsq, type, Gall, dGdx);
    
    //ANN feedfoward & feedback
    net.compute(Gall, E0, dEdG);
    
    //scaling
    E += Escale*E0 + Eshift + Eatom;
    for (int n = 0; n < nGtot; n++) dEdG[n] = Escale*dEdG[n];
    
    for (int j = 0; j < jnum; j++){
      
      double *fj = &f[j*3];
      int jshift = j*nGtot*3;
      
      for (int n = 0; n < nGtot; n++){
        
        double *dGdxj = &dGdx[jshift + n*3];
        for(int i = 0; i < 3; i++)fj[i] -= dEdG[n]*dGdxj[i];
        
      }
    }
    
    
    delete [] Gall;
    delete [] dEdG;
    delete [] dGdx;
    
  }

};

#endif
