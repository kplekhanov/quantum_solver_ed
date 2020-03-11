#ifndef HILBERT_BONES_H
#define HILBERT_BONES_H

#include "headers.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  class HilbertBones{
  public:
	HilbertBones();
	virtual const HilbertBones* getHilPtr() const = 0;
	// size of the Hilbert space (in the symmetry sector for HilbertSym)
	virtual uword getNumStates() const = 0;
  
	// size of the lattice
	virtual uword getNumSites() const = 0;
  
	// size of the local Hilbert space at site i (for generic HS)
	virtual uword getLocalhsSize(uword i) const = 0;
  
	// 0 if don't use symmetries, 1 otherwise (for genericOp)
	virtual bool checkSymmetric() const = 0;
  
	// to be used for printing functions (operators should work with confints)
	virtual urowvec confInt2confVec(const uword& confint) const = 0;
	virtual uword confVec2confInt(const urowvec& confvec) const = 0;
  
	// to be used by operators on a generic HS (like Bose Hubbard but not hcb)
	virtual uword readQnAtSite(uword i, const uword& confint) const = 0;
	virtual uword changeQnAtSite(uword i, uword qn, const uword& confint) const = 0;
  
	// some functions to get index -> confint and vice versa
	virtual uword getConfInt(uword i) const = 0;
	virtual uword getConfPair(const uword& confint) const = 0;
  protected:
	uword num_states, num_sites;
  };

  HilbertBones::HilbertBones(){}

} //* namespace quantum_solver_ed

  
#endif //* HILBERT_BONES_H
