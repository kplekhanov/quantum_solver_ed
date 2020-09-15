#include "../lib/symmetry.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  Symmetry class          *** //

  template<typename T>
  Symmetry<T>::Symmetry(){}
  
  template<typename T>
  Symmetry<T>::Symmetry(uword n_sites){
	this->n_sites = n_sites;
	this->perm = urowvec(n_sites);
	for (uword i=0; i<n_sites; i++)
	  this->perm[i] = i;
	this->chi = 1;
	this->name = "id" + to_string(this->n_sites);
  }

  template<typename T>
  void Symmetry<T>::setName(string new_name){
	this->name = new_name;
  }

  template<typename T>
  void Symmetry<T>::setPerm(const urowvec& other_perm){
	this->perm = other_perm;
  }

  template<typename T>
  void Symmetry<T>::setChi(T new_chi){
	this->chi = new_chi;
  }

  template<typename T>
  urowvec Symmetry<T>::apply(const urowvec& conf_vec_in) const{
	uword n_sites = conf_vec_in.n_elem;
	if (n_sites != this->n_sites){
	  throw logic_error("In Symmetry::apply. Size problem.");
	}
	urowvec conf_vec_out(n_sites);
	for (uword i=0; i<n_sites; i++)
	  conf_vec_out[this->perm[i]] = conf_vec_in[i];
	return conf_vec_out;
  }

  template<typename T>
  void Symmetry<T>::print(ostream& os) const{
	os << this->getName() << " |";
	for (uword i=0; i<this->getNumSites(); i++)
	  os << " " << (*this)[i];
	os << " | " << this->getChi() <<  " |";
  }

  template<typename T>
  void Symmetry<T>::fill(vector<uword> indices_int, T amplitude_val){
	if (n_sites != indices_int.size()){
	  throw logic_error("In Symmetry::fill. Number of index entries shoud be equal to n_sites.");
	}
	for (uword i=0; i<indices_int.size(); ++i){
	  	this->perm[i] = indices_int[i];
	}
	this->chi = amplitude_val;
  }

  // *** ------------------------ *** //

  // *** ------------------------ *** //
  // *** instatation of template  *** //
  
  template class Symmetry<double>;
  template class Symmetry<cx_double>;

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
