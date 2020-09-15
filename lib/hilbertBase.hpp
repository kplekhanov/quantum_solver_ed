#ifndef HILBERT_BASE_H
#define HILBERT_BASE_H

#include "hilbertBones.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------------------------ *** //
  // *** 2nd base space for the Hilbert space class *** //

  // *** as compared to bones, it already has a well defined mask structure
  // *** it has localhs attricubtes which tell what is the local structure of the Hilbert space
  // *** based on these attributes, one can get universal methods
  // *** to read / change quantum numers at a given site
  // *** however, creation methods are not specified yet
  // *** the idea is that it is a rather generic Hilbert space
  // *** if one knows exactly the Hilbert psace structure, like HC bosons, for example
  // *** one could hard-code most of the methods in a more efficient way
  // *** this spans the triee of HilbertBrut, Hilbert, HilbertSym
  // *** which basically only differ by their creation methods
  class HilbertBase: public HilbertBones{
  public:
	HilbertBase(const urowvec& localhs_sizes_int);
	inline
	uword getNumStates() const {return this->num_states;}
	inline
	uword getNumSites() const {return this->num_sites;}
	inline
	uword getLocalhsSize(uword i) const {return this->localhs_sizes_int[i];}
	inline // polymorphism for HilbertSym
	virtual bool checkSymmetric() const {return false;}
	urowvec confInt2confVec(const uword& confint) const;
	uword confVec2confInt(const urowvec& confvec) const;
	uword readQnAtSite(uword i, const uword& confint) const;
	uword changeQnAtSite(uword i, uword qn, const uword& confint) const;
	void configureMasks();
	virtual void create() = 0;
	virtual void createWithFixedQn(uword qn_tot) = 0;
	virtual void print() const = 0;
  protected:
	uword num_states, num_sites;
	urowvec localhs_sizes_int;
	urowvec localhs_sizes_bit;
	urowvec localhs_shifts;
	urowvec localhs_masks, localhs_nonmasks;
  };

  // *** ------------------------------------------ *** //
  
} //* namespace quantum_solver_ed

  
#endif //* HILBERT_BASE_H
