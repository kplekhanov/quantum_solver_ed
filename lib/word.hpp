#ifndef WORD_H
#define WORD_H

#include "headers.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***       Word class         *** //

  // *** Word is basically a binary representation of the symmetric quantum state
  // *** in a form "rep state | degeneracy | sym index"
  // *** sym index tells whuch sym has to be applied to the rep state in order to get given state
  class Word{
  public:
	Word(uword wrd);
	Word(uword conf, uword deg, uword sym);
	static uword getWrd(uword rep, uword deg, uword sym);
	static uword getRep(uword wrd);
	inline
	uword getWrd() const {return this->w;}
	inline
	uword getRep() const {return this->w >> (2 * this->shift);}
	inline
	uword getDeg() const {return (this->w & this->mask_deg) >> this->shift;}
	inline
	uword getSym() const {return this->w & this->mask_sym;}
  private:
	uword w;
	static const uword shift = 8;
	static const uword mask_sym = (uword(1) << shift) - 1;
	static const uword mask_deg = mask_sym << shift;
  };

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed


#endif //* WORD_H
