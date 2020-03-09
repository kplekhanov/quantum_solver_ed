#ifndef OP_BOSE_HUBBARD_H
#define OP_BOSE_HUBBARD_H

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

  double Ni::apply(const uword& conf_in) const{
	return double(this->hil_ptr->readQnAtSite(this->i, conf_in));
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
	uword qni = this->hil_ptr->readQnAtSite(this->i, conf_in);
	uword qnj = this->hil_ptr->readQnAtSite(this->j, conf_in);
	if ((qni == this->hil_ptr->getLocalhsSize(this->i)-1) | (qnj == 0))
	  return PairUwordT<T>(0, 0);
	uword conf_new = this->hil_ptr->changeQnAtSite(this->i, qni+1, conf_in);
	conf_new = this->hil_ptr->changeQnAtSite(this->j, qnj-1, conf_new);
	return PairUwordT<T>(conf_new, sqrt(double(qnj*(qni+1))));
  }

  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_BOSE_HUBBARD_H
