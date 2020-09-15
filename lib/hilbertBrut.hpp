#ifndef HILBERT_BRUT_H
#define HILBERT_BRUT_H

#include "hilbertBase.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ---------------------------------------- *** //
  // *** The most basic full Hilbert space class  *** //

  // *** Hilbert space is computed with brute force (which explains the name)
  // *** in addition to base, has create methods as well as
  // *** confint_vec which is a map index -> quantum state
  // *** and confint_map which is a map quantum state -> index
  class HilbertBrut: public HilbertBase{
  public:
	HilbertBrut(const urowvec& localhs_sizes_int);
	inline // get quantum state at the index i
	uword getConfInt(uword i) const {return this->confint_vec[i];}
	inline // get index corresponding to confint; needed for generalOp
	uword getConfPair(const uword& confint) const {return this->confint_map.at(confint);}
	void create();
	void createWithFixedQn(uword qn_tot=0);
	void print() const;
  protected:
	uvec confint_vec;
	HilMap confint_map;
  };

  // *** ---------------------------------------- *** //
  
} //* namespace quantum_solver_ed


#endif //* HILBERT_BRUT_H
