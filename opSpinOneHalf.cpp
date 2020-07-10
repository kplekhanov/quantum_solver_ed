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
	Szi(){};
	Szi(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
  };

  Szi::Szi(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void Szi::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 1){
	  throw logic_error( "In Szi::fill. Vector of indices should have size 1." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->name = "Sz_" + to_string(this->i);
  }

  inline
  double Szi::apply(const uword& conf_in) const{
	return double((conf_in >> i) & 1) - 0.5;
  }

  // *** Exchange along the z-axis
  class SziSzj: public DiagOp{
  public:
	SziSzj(){};
	SziSzj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
	uword j;
  };

  SziSzj::SziSzj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void SziSzj::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSzj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sz_" + to_string(this->j);
  }

  inline
  double SziSzj::apply(const uword& conf_in) const{
	return (double((conf_in >> i) & 1) - 0.5) * (double((conf_in >> j) & 1) - 0.5);
  }   

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** swap gates
  template<typename T=double>
  class SpiSmj: public OffDiagOp<T>{
  public:
	SpiSmj(){};
	SpiSmj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SpiSmj<T>::SpiSmj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SpiSmj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SpiSmj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sp_" + to_string(this->i) + ".Sm_" + to_string(this->j);
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
