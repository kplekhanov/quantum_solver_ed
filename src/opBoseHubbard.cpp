#include "../lib/opBoseHubbard.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  Ni::Ni(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void Ni::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 1){
	  throw logic_error( "In Ni::fill. Vector of indices should have size 1." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->name = "N_" + to_string(this->i);
  }
  
  double Ni::apply(const uword& conf_in) const{
	return double(this->hil_ptr->readQnAtSite(this->i, conf_in));
  }

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** hopping operators
  template<typename T>
  BdagiBj<T>::BdagiBj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void BdagiBj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In BdagiBj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Bdag_" + to_string(this->i) + ".B_" + to_string(this->j);
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

  // *** ------------------------ *** //
  // *** instatation of template  *** //
  
  template class BdagiBj<double>;
  template class BdagiBj<cx_double>;

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
