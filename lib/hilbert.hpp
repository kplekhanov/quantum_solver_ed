#ifndef HILBERT_H
#define HILBERT_H

#include "hilbertBrut.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------------------------------- *** //
  // *** advanced complete Hilbert space class  *** //

  // *** as compared to brute, create uses two confint maps
  // *** the Hilbert space is split into two to accelerate confint search
  class Hilbert: public HilbertBase{
  public:
	Hilbert(const urowvec& localhs_sizes_int);
	inline
	uword getConfInt(uword i) const {return this->confint_vec[i];}
	uword getConfPair(const uword& confint) const;
	void create();
	void createWithFixedQn(uword qn_tot=0);
	void print() const;
  protected:
	uword half_shift, mask_left, mask_rght;
	uvec confint_vec;
	HilMap confint_map_left, confint_map_rght;
  };

  // *** -------------------------------------- *** //
  
} //* namespace quantum_solver_ed


#endif //* HILBERT_H
