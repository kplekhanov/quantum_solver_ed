#ifndef OP_SPIN_ONE_HALF_H
#define OP_SPIN_ONE_HALF_H

#include "elementaryOp.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** Zeeman along the z-axis
  class Szi: public DiagOp{
  public:
	Szi(const HilbertBones* hil_ptr, uword i);
	double apply(const uword& conf_in) const;
  private:
	uword i;
  };

  Szi::Szi(const HilbertBones* hil_ptr, uword i){
	this->hil_ptr = hil_ptr;
	this->i = i;
	this->name = "Sz_" + to_string(i);
  }

  inline
  double Szi::apply(const uword& conf_in) const{
	return double((conf_in >> i) & 1) - 0.5;
  }

  // *** Exchange along the z-axis
  class SziSzj: public DiagOp{
  public:
	SziSzj(const HilbertBones* hil_ptr, uword i, uword j);
	double apply(const uword& conf_in) const;
  private:
	uword i;
	uword j;
  };

  SziSzj::SziSzj(const HilbertBones* hil_ptr, uword i, uword j){
	this->hil_ptr = hil_ptr;
	this->i = i;
	this->j = j;
	this->name = "Sz_" + to_string(i) + "Sz_" + to_string(j);
  }

  inline
  double SziSzj::apply(const uword& conf_in) const{
	return (double((conf_in >> i) & 1) - 0.5) * (double((conf_in >> j) & 1) - 0.5);
  }   

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** sweep gates
  template<typename T=double>
  class SpiSmj: public OffDiagOp<T>{
  public:
	SpiSmj(const HilbertBones* hil_ptr, uword i, uword j);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SpiSmj<T>::SpiSmj(const HilbertBones* hil_ptr, uword i, uword j){
	this->hil_ptr = hil_ptr;
	this->i = i;
	this->j = j;
	this->name = "Sp_" + to_string(i) + ".Sm_" + to_string(j);
  }

  template<typename T>
  PairUwordT<T> SpiSmj<T>::apply(const uword& conf_in) const{
	if ( (((conf_in >> i) & 1) == 0) & (((conf_in >> j) & 1) == 1) )
	  return PairUwordT<T>((conf_in & ~(1 << j)) | (1 << i), double(1));
	else
	  return PairUwordT<T>(0 , 0);
  }

  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_SPIN_ONE_HALF_H
