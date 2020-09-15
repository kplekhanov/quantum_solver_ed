#ifndef ELEMENTARY_OP_H
#define ELEMENTARY_OP_H

#include "hilbertBones.hpp"
#include "headers.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ---------------------------- *** //
  // *** class of elemetary operators *** //

  // *** base class for polymorphism
  // *** ElementaryOp acts on a confint, i.e. a binary representation of the state
  // *** example: cdag_i c_j -- creates a particle at i and annihilates at j
  class ElementaryOp{
  public:
	ElementaryOp(){}
	inline
	const HilbertBones* getHilPtr() const {return this->hil_ptr;}
	void print(ostream& os) const;
	virtual void fill(const HilbertBones* hil_ptr, vector<uword> indices_int) = 0;
  protected:
	const HilbertBones* hil_ptr;
	string name;
  };

  // *** diagonal operators (return coef when applied)
  class DiagOp: public ElementaryOp{
  public:
	DiagOp(){}
	virtual double apply(const uword& conf_in) const = 0;
  };

  // *** off-diagonal operators (return (conf, coef) when applied)
  template<typename T>
  class OffDiagOp: public ElementaryOp{
  public:
	OffDiagOp(){}
	virtual PairUwordT<T> apply(const uword& conf_in) const = 0;
  };

  // *** functions
  ostream& operator<<(ostream& os, const ElementaryOp& eop);

  // *** ---------------- *** //

} //* namespace quantum_solver_ed


#endif //* ELEMENTARY_OP_H
