#include "../lib/opSpinOneHalf.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** Zeeman along the z-axis
  Szi_2::Szi_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void Szi_2::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 1){
	  throw logic_error( "In Szi::fill. Vector of indices should have size 1." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->name = "Sz_" + to_string(this->i);
  }

  // *** Exchange along the z-axis
  SziSzj_2::SziSzj_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void SziSzj_2::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSzj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sz_" + to_string(this->j);
  }

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** swap gates
  template<typename T>
  SpiSmj_2<T>::SpiSmj_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SpiSmj_2<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SpiSmj_2::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sp_" + to_string(this->i) + ".Sm_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SpiSmj_2<T>::apply(const uword& conf_in) const{
	if ( (((conf_in >> i) & 1) == 0) & (((conf_in >> j) & 1) == 1) )
	  return PairUwordT<T>((conf_in & ~(1 << j)) | (1 << i), double(1));
	else
	  return PairUwordT<T>(0 , 0);
  }


  // *** in-plane DM interactions -- Sz Sp
  template<typename T>
  SziSpj_2<T>::SziSpj_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSpj_2<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSpj_2::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sp_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSpj_2<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 0 )
	  return PairUwordT<T>(conf_in | (1 << j), double((conf_in >> i) & 1) - 0.5);
	else
	  return PairUwordT<T>(0 , 0);
  }

  
  // *** in-plane DM interactions -- Sz Sm
  template<typename T>
  SziSmj_2<T>::SziSmj_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSmj_2<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSmj_2::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sm_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSmj_2<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 1 )
	  return PairUwordT<T>(conf_in & ~(1 << j), double((conf_in >> i) & 1) - 0.5);
	else
	  return PairUwordT<T>(0 , 0);
  }

  
  // *** in-plane DM interactions -- Sz Sx
  template<typename T>
  SziSxj_2<T>::SziSxj_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSxj_2<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSxj_2::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sx_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSxj_2<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 0 )
	  return PairUwordT<T>(conf_in | (1 << j), (double((conf_in >> i) & 1) - 0.5) * 0.5 );
	else
	  return PairUwordT<T>(conf_in & ~(1 << j), (double((conf_in >> i) & 1) - 0.5) * 0.5 );
  }
  

  // *** in-plane DM interactions -- Sz Sy
  template<typename T>
  SziSyj_2<T>::SziSyj_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (is_same<T, double>::value){
	  throw logic_error( "In SziSyj_2::SziSyj_2. This operator can't be real" );
	}
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSyj_2<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSyj_2::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sy_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSyj_2<T>::apply(const uword& conf_in) const{
	if ( ((conf_in >> j) & 1) == 0 ){
	  return PairUwordT<T>(conf_in | (1 << j),
						   cx_double(0,-(double((conf_in >> i) & 1) - 0.5) * 0.5));
	}
	else
	  return PairUwordT<T>(conf_in & ~(1 << j),
						   cx_double(0, (double((conf_in >> i) & 1) - 0.5) * 0.5));
  }
  
  // *** ------------------ *** //

  // *** ------------------------ *** //
  // *** instatation of template  *** //
  
  template class SpiSmj_2<double>;
  template class SpiSmj_2<cx_double>;
  template class SziSpj_2<double>;
  template class SziSpj_2<cx_double>;
  template class SziSmj_2<double>;
  template class SziSmj_2<cx_double>;
  template class SziSxj_2<double>;
  template class SziSxj_2<cx_double>;
  template class SziSyj_2<cx_double>;

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
