#ifndef HILBERT_SYM_H
#define HILBERT_SYM_H

#include "hilbertBrut.hpp"
#include "symmetry.hpp"
#include "word.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------------------- *** //
  // *** the most sdvanced of Hilbert classes  *** //

  // *** as compared to others, uses symmetries
  // *** confint_vec maps index to a word corresponding to a symmetric confint
  // *** however, when generalOp acts on a confint, in general it does not produce
  // *** symmetric word as in confint_vec
  // *** hence, one needs to keep track of all the states in the original Hilbert space
  template <typename T=cx_double>
  class HilbertSym: public HilbertBase{
	typedef vector<Symmetry<T>*> VectorSymPtr;
  public:
	HilbertSym(const urowvec& localhs_sizes_int);
	inline
	bool checkSymmetric() const {return true;}
	inline
	void addSym(Symmetry<T>* symptr) {this->symptr_vec.push_back(symptr);}
	void readSymFromFile(string file_name, char delimiter=';');
	inline
	uword getConfInt(uword i) const {return this->confint_vec[i];}
	uword getConfPair(const uword& confint) const;
	PairUwordT<T> getSymPair(const uword& confint) const;
	PairUword getRep(const uword& confint) const;
	int testRep(const uword& confint) const;
	void fill(uword index, const uword& confint, HilMap* confint_map_tmp_ptr);
	void create();
	void createWithFixedQn(uword qn_tot=0);
	void print() const;
	void printSymmetries() const;
  private:
	uword num_states_nonsym;
	uword tot_shift, half_shift, mask_left, mask_rght;
	VectorSymPtr symptr_vec;
	uvec confint_vec, rep_lookup_vec;
	HilMap confint_map_left, confint_map_rght;
  };

  // *** ------------------------------------- *** //
  
} //* namespace quantum_solver_ed


#endif //* HILBERT_SYM_H
