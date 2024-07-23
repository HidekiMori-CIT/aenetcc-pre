
#ifdef PAIR_CLASS
// clang-format off
PairStyle(aenet/cc,PairAENET_CC);
// clang-format on
#else

#ifndef MLP_PAIR_AENETCC_H
#define MLP_PAIR_AENETCC_H

#include "aenet_pot.h"
#include "pair.h"
#include "weight_functions.h"

namespace LAMMPS_NS {

class PairAENET_CC : public Pair {
 protected:
   virtual void allocate();
   AENET_POT *aenet_pot = nullptr;
   int load_pot_file(int, int, std::string, std::vector<std::string>);
   double Rc_max_sq;
   
   bool Eshifted_on = false;
   bool Alchemic_on = false;

   bool Econnect_on;
   Weight_Functions EW_Function;
   
 public:
  PairAENET_CC(class LAMMPS *);
  ~PairAENET_CC() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;


};

}

#endif
#endif

/* ERROR/WARNING messages:


*/
