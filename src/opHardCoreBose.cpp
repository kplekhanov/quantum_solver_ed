#include "../lib/opHardCoreBose.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  Ni_2::Ni_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void Ni_2::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 1){
	  throw logic_error( "In Ni_2::fill. Vector of indices should have size 1." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->name = "N_" + to_string(this->i);
  }

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** hopping operators
  template<typename T>
  BdagiBj_2<T>::BdagiBj_2(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void BdagiBj_2<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In BdagiBj_2::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Bdag_" + to_string(this->i) + ".B_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> BdagiBj_2<T>::apply(const uword& conf_in) const{
	if ( (((conf_in >> i) & 1) == 0) & (((conf_in >> j) & 1) == 1) )
	  return PairUwordT<T>((conf_in & ~(1 << j)) | (1 << i), double(1));
	else
	  return PairUwordT<T>(0, 0);
  }

  // *** ------------------ *** //

  // *** ------------------------ *** //
  // *** instatation of template  *** //
  
  template class BdagiBj_2<double>;
  template class BdagiBj_2<cx_double>;

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
