#include "../lib/opFermi.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  Ni_f::Ni_f(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void Ni_f::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 1){
	  throw logic_error( "In Ni_f::fill. Vector of indices should have size 1." );
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
  BdagiBj_f<T>::BdagiBj_f(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void BdagiBj_f<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In BdagiBj_f::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->start = min(this->i ,this->j);
	this->end = max(this->i ,this->j);
	this->name = "Bdag_" + to_string(this->i) + ".B_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> BdagiBj_f<T>::apply(const uword& conf_in) const{
	if ( (((conf_in >> i) & 1) == 1) || (((conf_in >> j) & 1) == 0) )
	  return PairUwordT<T>(0, 0);
	double jw_string = 1;
	for (uword k=this->start; k<this->end; ++k)
	  jw_string *= pow(-1, (conf_in >> k) & 1);
	return PairUwordT<T>((conf_in & ~(1 << j)) | (1 << i), jw_string);
  }

  // *** ------------------ *** //

  // *** ------------------------ *** //
  // *** instatation of template  *** //
  
  template class BdagiBj_f<double>;
  template class BdagiBj_f<cx_double>;

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
