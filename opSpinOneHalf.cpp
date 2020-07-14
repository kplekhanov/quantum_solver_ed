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


  // *** in-plane DM interactions -- Sz Sp
  template<typename T=double>
  class SziSpj: public OffDiagOp<T>{
  public:
	SziSpj(){};
	SziSpj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SziSpj<T>::SziSpj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSpj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSpj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sp_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSpj<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 0 )
	  return PairUwordT<T>(conf_in | (1 << j), double((conf_in >> i) & 1) - 0.5);
	else
	  return PairUwordT<T>(0 , 0);
  }

  
  // *** in-plane DM interactions -- Sz Sm
  template<typename T=double>
  class SziSmj: public OffDiagOp<T>{
  public:
	SziSmj(){};
	SziSmj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SziSmj<T>::SziSmj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSmj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSmj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sm_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSmj<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 1 )
	  return PairUwordT<T>(conf_in & ~(1 << j), double((conf_in >> i) & 1) - 0.5);
	else
	  return PairUwordT<T>(0 , 0);
  }

  
  // *** in-plane DM interactions -- Sz Sx
  template<typename T=double>
  class SziSxj: public OffDiagOp<T>{
  public:
	SziSxj(){};
	SziSxj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SziSxj<T>::SziSxj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSxj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSxj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sx_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSxj<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 0 )
	  return PairUwordT<T>(conf_in | (1 << j), (double((conf_in >> i) & 1) - 0.5) * 0.5 );
	else
	  return PairUwordT<T>(conf_in & ~(1 << j), (double((conf_in >> i) & 1) - 0.5) * 0.5 );
  }
  

  // *** in-plane DM interactions -- Sz Sy
  template<typename T>
  class SziSyj: public OffDiagOp<T>{
  public:
	SziSyj(){};
	SziSyj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SziSyj<T>::SziSyj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (is_same<T, double>::value){
	  throw logic_error( "In SziSyj::SziSyj. This operator can't be real" );
	}
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSyj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSyj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sy_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSyj<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 0 ){
	  return PairUwordT<T>(conf_in | (1 << j),
						   cx_double(0,-(double((conf_in >> i) & 1) - 0.5) * 0.5));
	}
	else
	  return PairUwordT<T>(conf_in & ~(1 << j),
						   cx_double(0, (double((conf_in >> i) & 1) - 0.5) * 0.5));
  }
  
  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_SPIN_ONE_HALF_H
