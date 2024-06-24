
#ifdef PAIR_CLASS
// clang-format off
PairStyle(aenet/cc,PairAENETCC);
// clang-format on
#else

#ifndef MLP_PAIR_AENETCC_H
#define MLP_PAIR_AENETCC_H

#include "aenet_pot.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairAENETCC : public Pair {
 private:
   AENET_POT *aenet_pot = nullptr;
   void load_pot_file(int, int, std::string, int, char**);
   double Rc_max_sq;
   
   bool Eshifted_on = false;
   bool Alchemic_on = false;
   
 public:
  PairAENETCC(class LAMMPS *);
  ~PairAENETCC() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  virtual void allocate();

};

}

#endif
#endif

/* ERROR/WARNING messages:


*/
