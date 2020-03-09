#ifndef WORD_H
#define WORD_H

#include "headers.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***       Word class         *** //

  class Word{
  public:
	Word(uword wrd);
	Word(uword conf, uword deg, uword sym);
	static uword getWrd(uword rep, uword deg, uword sym);
	static uword getRep(uword wrd);
	uword getWrd() const;
	uword getRep() const;
	uword getDeg() const;
	uword getSym() const;
  private:
	uword w;
	static const uword shift = 8;
	static const uword mask_sym = (uword(1) << shift) - 1;
	static const uword mask_deg = mask_sym << shift;
  };

  Word::Word(uword wrd){
	this->w = wrd;
  }

  Word::Word(uword rep, uword deg, uword sym){
	this->w = rep << (2 * this->shift);
	this->w = this->w | deg << this->shift;
	this->w = this->w | sym;
  }

  uword Word::getWrd(uword rep, uword deg, uword sym){
	uword w = rep << (2 * shift);
	w = w | deg << shift;
	w = w | sym;
	return w;
  }

  inline
  uword Word::getRep(uword w){
	return w >> (2 * shift);
  }

  inline
  uword Word::getWrd() const{
	return this->w;
  }

  inline
  uword Word::getRep() const{
	return this->w >> (2 * this->shift);
  }

  inline
  uword Word::getDeg() const{
	return (this->w & this->mask_deg) >> this->shift;
  }

  inline
  uword Word::getSym() const{
	return this->w & this->mask_sym;
  }

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed


#endif //* WORD_H
