#ifndef ELEMENTARY_OP_H
#define ELEMENTARY_OP_H

#include "headers.hpp"
#include "hilbertBones.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ---------------- *** //
  // *** ops act on confs *** //

  // *** base class for polymorphism
  class ElementaryOp{
  public:
	const HilbertBones* getHilPtr() const;
	void print(ostream& os) const;
  protected:
	const HilbertBones* hil_ptr;
	string name;
  };

  inline
  const HilbertBones* ElementaryOp::getHilPtr() const{
	return this->hil_ptr;
  }

  void ElementaryOp::print(ostream& os) const{
	os << this->name;
  }

  inline
  ostream& operator<<(ostream& os, const ElementaryOp& eop){
	eop.print(os);
	return os;
  }


  // *** diagonal operators (return coef when applied)
  class DiagOp: public ElementaryOp{
  public:
	virtual double apply(const uword& conf_in) const = 0;
  };


  // *** off-diagonal operators (return (conf, coef) when applied)
  template<typename T>
  class OffDiagOp: public ElementaryOp{
  public:
	virtual PairUwordT<T> apply(const uword& conf_in) const = 0;
  };

  // *** ---------------- *** //

} //* namespace quantum_solver_ed


#endif //* ELEMENTARY_OP_H
