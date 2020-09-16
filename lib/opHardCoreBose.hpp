#ifndef OP_HARD_CORE_BOSE_H
#define OP_HARD_CORE_BOSE_H

#include "elementaryOp.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  class Ni_2: public DiagOp{
  public:
	Ni_2(){}
	Ni_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
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
  class BdagiBj_2: public OffDiagOp<T>{
  public:
	BdagiBj_2(){}
	BdagiBj_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_HARD_CORE_BOSE_H
