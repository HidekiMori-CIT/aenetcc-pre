
#ifdef PAIR_CLASS

PairStyle(aenet/cc,PairAENETCC)

#else

#ifndef MLP_PAIR_AENETCC_H
#define MLP_PAIR_AENETCC_H

#include "aenet_pot.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairAENETCC : public Pair {
 private:
   AENET_POT *aenet_pot;
   void load_pot_file(int, int, std::string, int, char**);
   
 public:
  PairAENETCC(class LAMMPS *);
  ~PairAENETCC();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  virtual void allocate();
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements

  int *map;                     // mapping from atom types to elements

};

}

#endif
#endif

/* ERROR/WARNING messages:


*/
