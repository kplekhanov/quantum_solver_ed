#include "../lib/word.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***       Word class         *** //

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

  uword Word::getRep(uword w){
	return w >> (2 * shift);
  }

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
