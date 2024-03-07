#ifndef ANN_H
#define ANN_H

#include <vector>
#include <cmath>

class ANN {
  
  int nlayers;
  int *activation_ftype;
  int *nnodes;
  double *Wb;
  double *TW;
  int Wb_size = 0;
  int TW_size = 0;
  int nn_input;
  
  double *Gscale;
  double *Gshift;
  
  double **nodes, **dnodes, **bnodes;

  double activation(int f_a, double x, double &deriv)
  {
    
    double y,dy;
    double tanhbx,sigmoid_bx,softplus_x,omega,delta,ysq;
    
    const double a = 1.7159;
    const double b = 0.666666666666667;
    const double c = 0.1;
    const double beta = 1.0;
    const double eps = 0.25/4.0;
    
    if (f_a == 0){
      y  = x;
      dy = 1.0;
    }else if (f_a == 1){
      y  = tanh(x);
      dy = 1.00 - y*y;
    }else if (f_a == 2){
      y  = 1.00/(1.00 + exp(-x));
      dy = y*(1.00 - y);
    }else if (f_a == 3){
      tanhbx = tanh(b*x);
      y  = a*tanhbx;
      dy = a*(1.0 - tanhbx*tanhbx)*b;
    }else if (f_a == 4){
      tanhbx = tanh(b*x);
      y  = a*tanhbx + c*x;
      dy = a*(1.0 - tanhbx*tanhbx)*b + c;
    }else if (f_a == 5){
      sigmoid_bx  = (1.0/(1.0 + exp(-beta*x)));
      y  = x*sigmoid_bx;
      dy = beta*y + sigmoid_bx*(1.0 - beta*y);
    }else if (f_a == 6){
      softplus_x = log(1+exp(x));
      y  = x*tanh(softplus_x);
      omega = 4*(x+1) + 4*exp(2*x) + exp(3*x) + exp(x)*(4*x+6);
      delta = 2*exp(x) + exp(2*x) + 2;
      dy = exp(x)*omega/(delta*delta);
    }else if (f_a == 7){
      y  = 0.5*(x + sqrt(x*x + 4.0*eps));
      ysq = y*y;
      dy = ysq/(ysq + eps);
    }else {
      y  = 0.0;
      dy = 0.0;
    }
    
    deriv = dy;
    
    return y;
  }
  
  public:
  
  ~ANN(){
    
    delete [] activation_ftype;
    delete [] nnodes;
    delete [] Wb;
    delete [] TW;
    
    for (int n = 0; n < nlayers; n++) {
     delete [] nodes[n];
     delete [] dnodes[n];
     delete [] bnodes[n];
    }
    
    delete [] nodes;
    delete [] dnodes;
    delete [] bnodes;
    
    delete [] Gscale;
    delete [] Gshift;
    
  }
 
  void compute(double *G0, double &E, double *dEdG)
  {
   
    int nG = nn_input;
    
    double *G = new double [nG]{};
    
    //scaling
    for (int i = 0; i < nG; i++)G[i] = Gscale[i]*(G0[i] - Gshift[i]);
    
    
    // feedfwoard
    
    // input G (layer0) to hidden layer1
    for (int n = 0; n < nnodes[0]; n++) {
      nodes[0][n] = Wb[n*(nG+1)];
      for (int i = 0; i < nG; i++)nodes[0][n] += Wb[n*(nG+1) + 1 + i] * G[i];
      nodes[0][n] = activation(activation_ftype[0], nodes[0][n], dnodes[0][n]);
    }
    
    // hidden layers ~output layer
    int k = 0;
    if (nlayers > 1) {
      k += (nG+1)*nnodes[0];
        for (int l = 1; l < nlayers; l++) {
          for (int n = 0; n < nnodes[l]; n++) {
            nodes[l][n] = Wb[k + n*(nnodes[l-1]+1)];
            for (int j = 0; j < nnodes[l-1]; j++)nodes[l][n] += Wb[k + n*(nnodes[l-1]+1) + 1 + j]*nodes[l-1][j];
            nodes[l][n] = activation(activation_ftype[l], nodes[l][n], dnodes[l][n]);
          }
          k += (nnodes[l-1]+1)*nnodes[l];
        }
    }
    
    //set atomic energy
    E = nodes[nlayers - 1][0];
    
    
    // feedback
    
    for (int n = 0; n < nnodes[nlayers-1]; n++){bnodes[nlayers-1][n] = dnodes[nlayers-1][n];}
    
    if (nlayers > 1) {
      k = TW_size;
      for (int l = nlayers-1; l > 0; l--) {
        k -= (nnodes[l-1])*nnodes[l];
        for (int n = 0; n < nnodes[l-1]; n++) {
          bnodes[l-1][n] = 0;
          for (int j = 0; j < nnodes[l]; j++)bnodes[l-1][n] += TW[k + n*(nnodes[l]) + j] * bnodes[l][j];
          bnodes[l-1][n] *= dnodes[l-1][n];
        }
      }
    }
    
    //set dEdG
    for (int i = 0; i < nG; i++)for (int j = 0; j < nnodes[0]; j++)dEdG[i] += TW[i*(nnodes[0]) + j] * bnodes[0][j];
    
    
    //scaling
    for (int i = 0; i < nG; i++)dEdG[i] = Gscale[i]*dEdG[i];
    
    delete [] G;
    
  }
  
  void set_ann(std::vector<int> nn, std::vector<int> f_a, std::vector<double> Wb_0){
    
    nlayers = nn.size()-1;
    
    nnodes  = new int[nlayers];
    activation_ftype = new int[nlayers];
    nn_input = nn[0];
    for (int l = 0; l < nlayers; l++){nnodes[l] = nn[l+1];}
    for (int l = 0; l < nlayers; l++){activation_ftype[l] = f_a[l];}
    
    nodes  = new double *[nlayers];
    dnodes = new double *[nlayers];
    bnodes = new double *[nlayers];
    for (int l = 0; l < nlayers; ++l) {
      nodes[l]  = new double [nnodes[l]];
      dnodes[l] = new double [nnodes[l]];
      bnodes[l] = new double [nnodes[l]];
    }
    
    Wb_size = 0;
    TW_size = 0;
    for (int l = 0; l < nlayers; l++){
      Wb_size += (nn[l]+1)*nnodes[l];
      TW_size += nn[l]*nnodes[l];
    }
    
    if(Wb_size != Wb_0.size()){
      std::cout<<"# of elements of Wb is different from input Wb"<<std::endl;
      std::exit(1);
    }
    
    Wb = new double[Wb_size]();
    TW = new double[TW_size]();
    
    for (int i = 0; i < Wb_size; i++) Wb[i] = Wb_0[i];
    
    int shift = 0;
    int TW_shift = 0;
    for (int l = 0; l < nlayers; l++){
      for (int i = 0; i < nn[l]; i++){
        for (int j = 0; j < nnodes[l]; j++){
          TW[i*(nnodes[l]) + j + TW_shift] = Wb[j*(1+nn[l]) + 1 + i + shift];
        } 
      }
      shift += nn[l]*nnodes[l] + nnodes[l];
      TW_shift += nn[l]*nnodes[l];
    }
    
  }
  
  void set_scale(int nG, double *scale, double *shift){
    
    //# of array elements check
    if(nG != nn_input){
      std::cout<<"# of elements of G is different from # of elements of input layer"<<std::endl;
      std::exit(1);
    }
    
    
    Gscale = new double [nG];
    Gshift = new double [nG];
    
    for (int n=0; n < nG; n++){
      Gscale[n] = scale[n];
      Gshift[n] = shift[n];
    }
    
  
  }
 
};

  

#endif
