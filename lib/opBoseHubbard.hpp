#ifndef OP_BOSE_HUBBARD_H
#define OP_BOSE_HUBBARD_H

#include "elementaryOp.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  class Ni: public DiagOp{
  public:
	Ni(){}
	Ni(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
  };

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** hopping operators
  template<typename T=double>
  class BdagiBj: public OffDiagOp<T>{
  public:
	BdagiBj(){}
	BdagiBj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_BOSE_HUBBARD_H
