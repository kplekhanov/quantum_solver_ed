#ifndef OP_HARD_CORE_BOSE_H
#define OP_HARD_CORE_BOSE_H

#include "elementaryOp.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  class Ni: public DiagOp{
  public:
	Ni(const HilbertBones* hil_ptr, uword i);
	double apply(const uword& conf_in) const;
  private:
	uword i;
  };

  Ni::Ni(const HilbertBones* hil_ptr, uword i){
	this->hil_ptr = hil_ptr;
	this->i = i;
	this->name = "N_" + to_string(i);
  }

  inline
  double Ni::apply(const uword& conf_in) const{
	return double((conf_in >> i) & 1);
  }

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** hopping operators
  template<typename T=double>
  class BdagiBj: public OffDiagOp<T>{
  public:
	BdagiBj(const HilbertBones* hil_ptr, uword i, uword j);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  BdagiBj<T>::BdagiBj(const HilbertBones* hil_ptr, uword i, uword j){
	this->hil_ptr = hil_ptr;
	this->i = i;
	this->j = j;
	this->name = "Bdag_" + to_string(i) + ".B_" + to_string(j);
  }

  template<typename T>
  PairUwordT<T> BdagiBj<T>::apply(const uword& conf_in) const{
	if ( (((conf_in >> i) & 1) == 0) & (((conf_in >> j) & 1) == 1) )
	  return PairUwordT<T>((conf_in & ~(1 << j)) | (1 << i), double(1));
	else
	  return PairUwordT<T>(0 , 0);
  }

  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_HARD_CORE_BOSE_H
