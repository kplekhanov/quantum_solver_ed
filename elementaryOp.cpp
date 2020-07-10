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
	ElementaryOp(){};
	const HilbertBones* getHilPtr() const;
	void print(ostream& os) const;
	virtual void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
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

  void ElementaryOp::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->hil_ptr = hil_ptr;
	cout << "Warning: filling is not implemented for this operator" << endl;
	cout << "Reading indices: ";
	for (auto i=indices_int.begin(); i<indices_int.end(); i++){
	  cout << *i << " ";
	}
	cout << endl;
  }

  inline
  ostream& operator<<(ostream& os, const ElementaryOp& eop){
	eop.print(os);
	return os;
  }


  // *** diagonal operators (return coef when applied)
  class DiagOp: public ElementaryOp{
  public:
	DiagOp(){};
	virtual double apply(const uword& conf_in) const = 0;
  };


  // *** off-diagonal operators (return (conf, coef) when applied)
  template<typename T>
  class OffDiagOp: public ElementaryOp{
  public:
	OffDiagOp(){};
	virtual PairUwordT<T> apply(const uword& conf_in) const = 0;
  };

  // *** ---------------- *** //

} //* namespace quantum_solver_ed


#endif //* ELEMENTARY_OP_H
