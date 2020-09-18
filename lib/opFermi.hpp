#ifndef OP_FERMI_H
#define OP_FERMI_H

#include "elementaryOp.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** Fermioninc operators
  // *** the only difference wrt bosonic ones
  // *** is the presence of a fermionic sign -- the JW string
  // *** which ensures correct fermionic commutation relations
  
  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  class Ni_f: public DiagOp{
  public:
	Ni_f(){}
	Ni_f(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	inline
	double apply(const uword& conf_in) const {return double((conf_in >> i) & 1);}
  private:
	uword i;
  };

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** hopping operators
  template<typename T=double>
  class BdagiBj_f: public OffDiagOp<T>{
  public:
	BdagiBj_f(){}
	BdagiBj_f(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j, start, end;
  };

  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_FERMI_H
