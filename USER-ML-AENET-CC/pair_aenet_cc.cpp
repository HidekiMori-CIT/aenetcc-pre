
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include "pair_aenet_cc.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "potential_file_reader.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairAENET_CC::PairAENET_CC(LAMMPS *lmp) : Pair(lmp)
{

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;

  aenet_pot = nullptr;

  Eshifted_on = false;
  Alchemic_on = false;



}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairAENET_CC::~PairAENET_CC()
{
  if(copymode) return;
  int i, stat;
  char Error_message[256];

  delete [] aenet_pot;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairAENET_CC::compute(int eflag, int vflag)
{

  ev_init(eflag,vflag);
  //if (eflag || vflag) ev_setup(eflag,vflag);
  //else evflag = vflag_fdotr = vflag_atom = 0;

  double **x = atom->x;
  double **f = atom->f;
  int  *type = atom->type;
  
  int nlocal      = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  int   inum = list->inum;
  int *ilist = list->ilist;
  
  int    *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  
  
  for (int ii = 0; ii < inum; ii++) {
    
    int i      = ilist[ii];
    int itype  = map[type[i]];
    
    if(aenet_pot[itype].skip())continue;

    int   jnum = numneigh[i];
    int *jlist = firstneigh[i]; 
    
    int    *indx_j = new int[jnum]{};
    int    *type_j = new int[jnum]{};
    double    *x_j = new double[jnum*3]{};
    double  *rsq_j = new double[jnum]{};
    
    int num_j=0;
    
    for (int jj = 0; jj < jnum; jj++) {
      
      int j = jlist[jj];
      j &= NEIGHMASK;
      
      double xx[3];
      double rr = 0.0;

      for (int k = 0; k < 3; k++)xx[k] = x[j][k] - x[i][k];
      for (int k = 0; k < 3; k++)rr   += xx[k]*xx[k];

      if(rr > Rc_max_sq)continue;
      
      indx_j[num_j] = j;
      type_j[num_j] = map[type[j]];
      
      double *xj = &x_j[num_j*3];
      for (int k = 0; k < 3; k++)xj[k] = xx[k];

      rsq_j[num_j]  = rr;
      
      num_j = num_j + 1;
      
    }
    
    double E_i = 0.0;   
    double *f_j = new double[num_j*3]{};
    
    aenet_pot[itype].compute(num_j, x_j, rsq_j, type_j, E_i, f_j);
    double Einf = aenet_pot[itype].get_Einf();    
    if(Econnect_on)EW_Function.select_auto_min(num_j, x_j, rsq_j, f_j, E_i, Einf);
    if(Eshifted_on)E_i -= Einf; 
    
    if (eflag) {
      if (eflag_global) eng_vdwl += E_i;
      if (eflag_atom  ) eatom[i] += E_i;
    }
    
    for (int jj = 0; jj < num_j; jj++){
      
      int j  = indx_j[jj];
      
      double *xj = &x_j[jj*3];
      double *fj = &f_j[jj*3];
      
      for (int kk = 0; kk < 3; kk++)f[i][kk] -= fj[kk];
      for (int kk = 0; kk < 3; kk++)f[j][kk] += fj[kk];
      
      if(evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,0.0,
                              fj[0],fj[1],fj[2],
                              xj[0],xj[1],xj[2]);
      
    }
    
    delete [] indx_j;
    delete [] type_j;
    delete [] x_j;
    delete [] f_j;
    delete [] rsq_j;
    
  }
  
  if (vflag_fdotr) virial_fdotr_compute();
  
  
}

/* ---------------------------------------------------------------------- */

void PairAENET_CC::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[n+1];

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAENET_CC::settings(int narg, char **arg)
{
  if (narg == 0){}
  else if (narg < 6 ) {    
    int iarg = 0;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"Eshifted") == 0){
        Eshifted_on  = true;
        iarg += 1;
      } else if (strcmp(arg[iarg],"Alchemic") == 0){
        Alchemic_on  = true;
        iarg += 1;
      } else if (strcmp(arg[iarg],"Econnect") == 0) {
        Econnect_on = true;
        double R0 = utils::numeric(FLERR, arg[iarg+1], false, lmp);
        double R1 = utils::numeric(FLERR, arg[iarg+2], false, lmp);
        EW_Function.setting(R0, R1);
        iarg += 3;
      } else error->all(FLERR,"Incorrect args for pair coefficients");
    }
  } else error->all(FLERR,"Illegal pair_style command");
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAENET_CC::coeff(int narg, char **arg)
{

  if (narg < 5) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  int aenet_version = atoi(std::string(arg[2]).substr(1,2).c_str());

  if (nelements) {
    for (int i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }

  nelements = narg - 4 - atom->ntypes;
  if (nelements < 0) error->all(FLERR,"Incorrect args for pair coefficients");

  std::string file_mask = std::string(arg[3+nelements]);

  if (nelements == 0) {

    map_element2type(narg-4,arg+4);

  } else {

    const int ntypes  = atom->ntypes;
    elements = new char*[ntypes];
    for (int i = 0; i < ntypes; i++) elements[i] = nullptr;
  
    for (int i = 0; i < nelements; i++) {
      std::string entry = arg[i+3];
      elements[i] = utils::strdup(entry);
    }

    // read args that map atom types to ANN elements
    // map[i] = which element the Ith atom type is, -1 if not mapped  
    for (int i = 4 + nelements; i < narg; i++) {
      int m = i - (4+nelements) + 1;
      int j;
      for (j = 0; j < nelements; j++)
        if (strcmp(arg[i],elements[j]) == 0) break;
      if (j < nelements) map[m] = j;
      else if (strcmp(arg[i],"NULL") == 0) map[m] = -1;
      else error->all(FLERR,"Incorrect args for pair coefficients");
    }

    // clear setflag since coeff() called once with I,J = * *
    // set setflag i,j for type pairs where both are mapped to elements
    int count = 0;
    for (int i = 1; i <= ntypes; i++) {
      for (int j = i; j <= ntypes; j++) {
        setflag[i][j] = 0;
        if ((map[i] >= 0) && (map[j] >= 0)) {
          setflag[i][j] = 1;
          count++;
        }
      }
    }

    if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  }

  aenet_pot = new AENET_POT [nelements];
  
  std::vector<std::string> element_set;
  for (int i = 0; i < nelements; i++)element_set.push_back(std::string(elements[i]));
  
  int pot_num=0;
  std::string wc = "**";
  int len = wc.length();
  for (int i = 0; i < nelements; i++){
    std::string ith_file_name = file_mask;
    if (ith_file_name.find(wc) != std::string::npos) {
      ith_file_name.replace(ith_file_name.find(wc),len,element_set[i]);
    } else ith_file_name = element_set[i]+"."+ith_file_name;
    pot_num += load_pot_file(aenet_version, i, ith_file_name, element_set);
  }
  if (pot_num == 0) error->all(FLERR,"aenet/cc: parameter file does not exist");

}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAENET_CC::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style AENET requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style AENET requires newton pair on");

  // need a full neighbor list
  // request a full neighbor list

  //int irequest = neighbor->request(this,instance_me);
  //neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;

  neighbor->add_request(this, NeighConst::REQ_FULL);

}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAENET_CC::init_one(int /*i*/, int /*j*/)
{
  double Rc_max;

  Rc_max = 0.0;
  for (int i = 0; i < nelements; i++){
    if (Rc_max < aenet_pot[i].get_Rc_max()){
      Rc_max = aenet_pot[i].get_Rc_max();
    }
  }

  Rc_max_sq = Rc_max*Rc_max;

  return Rc_max;
}

/* ---------------------------------------------------------------------- */

int PairAENET_CC::load_pot_file(int aenet_version, int itype, std::string file_name, std::vector<std::string> elements){

    std::vector<int>    nlayers     ;
    std::vector<int>    nnodes_max  ;
    std::vector<int>    Wsize       ;
    std::vector<int>    nvalues     ;
    std::vector<int>    nnodes      ;
    std::vector<int>    f_a         ;
    std::vector<int>    iw          ;
    std::vector<int>    iv          ;
    std::vector<double> W           ;

    std::vector<char>   description ;
    std::vector<char>   atomtype    ;
    std::vector<int>    nenv        ;
    std::vector<char>   envtypes    ;
    std::vector<double> Rc_min      ;
    std::vector<double> Rc_max      ;
    std::vector<char>   sftype      ;
    std::vector<int>    nsf         ;
    std::vector<int>    nsfparam    ;
    std::vector<int>    sf          ;
    std::vector<double> sfparam     ;

    std::vector<int>    sfenv       ;
    std::vector<int>    neval       ;
    std::vector<double> sfval_min   ;
    std::vector<double> sfval_max   ;
    std::vector<double> sfval_avg   ;
    std::vector<double> sfval_cov   ;

    std::vector<char>   file        ;
    std::vector<int>    normalized  ;
    std::vector<double> E_scale     ;
    std::vector<double> E_shift     ;
    std::vector<int>    nTypes      ;
    std::vector<char>   typeName    ;
    std::vector<double> E_atom      ;
    std::vector<int>    nAtomsTot   ;
    std::vector<int>    nStrucs     ;
    std::vector<double> E_info      ;
    
    int skip_flag = 0;

    if (comm->me == 0) {

      std::ifstream fin(file_name, std::ios::in | std::ios::binary );

      if (!fin) {

        if (Alchemic_on) {
          skip_flag = 1;
          if(fin.is_open())fin.close();
        } else error->all(FLERR,"aenet/cc: parameter file does not exist");

      } else { 

        nlayers     = read_Fbinary<int>(fin);
        nnodes_max  = read_Fbinary<int>(fin);
        Wsize       = read_Fbinary<int>(fin);
        nvalues     = read_Fbinary<int>(fin);
        nnodes      = read_Fbinary<int>(fin);
        f_a         = read_Fbinary<int>(fin);
        iw          = read_Fbinary<int>(fin);
        iv          = read_Fbinary<int>(fin);
        W           = read_Fbinary<double>(fin);

        description = read_Fbinary<char>(fin);
        atomtype    = read_Fbinary<char>(fin);
        nenv        = read_Fbinary<int>(fin);
        envtypes    = read_Fbinary<char>(fin);
        Rc_min      = read_Fbinary<double>(fin);
        Rc_max      = read_Fbinary<double>(fin);
        sftype      = read_Fbinary<char>(fin);
        nsf         = read_Fbinary<int>(fin);
        nsfparam    = read_Fbinary<int>(fin);
        sf          = read_Fbinary<int>(fin);
        sfparam     = read_Fbinary<double>(fin);

        sfenv       = read_Fbinary<int>(fin);
        neval       = read_Fbinary<int>(fin);
        sfval_min   = read_Fbinary<double>(fin);
        sfval_max   = read_Fbinary<double>(fin);
        sfval_avg   = read_Fbinary<double>(fin);
        sfval_cov   = read_Fbinary<double>(fin);

        file        = read_Fbinary<char>(fin);
        normalized  = read_Fbinary<int>(fin);
        E_scale     = read_Fbinary<double>(fin);
        E_shift     = read_Fbinary<double>(fin);
        nTypes      = read_Fbinary<int>(fin);
        typeName    = read_Fbinary<char>(fin);
        E_atom      = read_Fbinary<double>(fin);
        nAtomsTot   = read_Fbinary<int>(fin);
        nStrucs     = read_Fbinary<int>(fin);
        E_info      = read_Fbinary<double>(fin);

        fin.close();

      }

   }

    MPI_Bcast(&skip_flag, 1, MPI_INT, 0, world);
    if(skip_flag == 1){
      aenet_pot[itype].skip_flag_on();
      return 0;
    }

    int v_size;

    v_size = sftype.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    sftype.resize(v_size);
    MPI_Bcast(&sftype[0], sftype.size(), MPI_CHAR, 0, world);
    std::string output_sftype(sftype.begin(),sftype.end());
    output_sftype = output_sftype.erase(output_sftype.find_last_not_of(" ") + 1);
    output_sftype = output_sftype.erase(0, output_sftype.find_first_not_of(" "));

    v_size = nlayers.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    nlayers.resize(v_size);
    MPI_Bcast(&nlayers[0], nlayers.size(), MPI_INT, 0, world);

    v_size = nnodes.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    nnodes.resize(v_size);
    MPI_Bcast(&nnodes[0], nnodes.size(), MPI_INT, 0, world);

    v_size = f_a.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    f_a.resize(v_size);
    MPI_Bcast(&f_a[0], f_a.size(), MPI_INT, 0, world);

    v_size = Wsize.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    Wsize.resize(v_size);
    MPI_Bcast(&Wsize[0], Wsize.size(), MPI_INT, 0, world);

    v_size = W.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    W.resize(v_size);
    MPI_Bcast(&W[0], W.size(), MPI_DOUBLE, 0, world);

    v_size = Rc_max.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    Rc_max.resize(v_size);
    MPI_Bcast(&Rc_max[0], Rc_max.size(), MPI_DOUBLE, 0, world);

    v_size = nenv.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    nenv.resize(v_size);
    MPI_Bcast(&nenv[0], nenv.size(), MPI_INT, 0, world);

    v_size = envtypes.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    envtypes.resize(v_size);
    MPI_Bcast(&envtypes[0], envtypes.size(), MPI_INT, 0, world);

    v_size = nsf.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    nsf.resize(v_size);
    MPI_Bcast(&nsf[0], nsf.size(), MPI_INT, 0, world);

    v_size = sfparam.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    sfparam.resize(v_size);
    MPI_Bcast(&sfparam[0], sfparam.size(), MPI_DOUBLE, 0, world);

    v_size = sfval_avg.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    sfval_avg.resize(v_size);
    MPI_Bcast(&sfval_avg[0], sfval_avg.size(), MPI_DOUBLE, 0, world);

    v_size = sfval_cov.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    sfval_cov.resize(v_size);
    MPI_Bcast(&sfval_cov[0], sfval_cov.size(), MPI_DOUBLE, 0, world);

    v_size = E_scale.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    E_scale.resize(v_size);
    MPI_Bcast(&E_scale[0], E_scale.size(), MPI_DOUBLE, 0, world);

    v_size = E_shift.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    E_shift.resize(v_size);
    MPI_Bcast(&E_shift[0], E_shift.size(), MPI_DOUBLE, 0, world);

    v_size = atomtype.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    atomtype.resize(v_size);
    MPI_Bcast(&atomtype[0], atomtype.size(), MPI_CHAR, 0, world);

    v_size = E_atom.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    E_atom.resize(v_size);
    MPI_Bcast(&E_atom[0], E_atom.size(), MPI_DOUBLE, 0, world);

    v_size = nTypes.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    nTypes.resize(v_size);
    MPI_Bcast(&nTypes[0], nTypes.size(), MPI_INT, 0, world);

    v_size = typeName.size();
    MPI_Bcast(&v_size, 1, MPI_INT, 0, world);
    typeName.resize(v_size);
    MPI_Bcast(&typeName[0], typeName.size(), MPI_CHAR, 0, world);

    aenet_pot[itype].set_pot(
      aenet_version,
      atomtype,
      nnodes,f_a,W,
      elements,
      nenv[0],envtypes,sfparam,sfval_avg,sfval_cov,
      E_scale[0],E_shift[0],E_atom,nTypes[0],typeName
    );

    return 1;
}



