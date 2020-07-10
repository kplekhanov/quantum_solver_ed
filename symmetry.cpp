#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "symmetryBones.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  Symmetry class          *** //

  template<typename T=cx_double>
  class Symmetry: public SymmetryBones<T>{
  public:
	Symmetry(uword n_sites);
	uword getNumSites() const;
	string getName() const;
	void setName(string new_name);
	urowvec getPerm() const;
	void setPerm(const urowvec& other_perm);
	T getChi() const;
	void setChi(T new_chi);
	const uword& operator[](uword pos) const;
	uword& operator[](uword pos);
	urowvec apply(const urowvec& conf_vec_in) const;
	void print(ostream& os) const;
	void fill(vector<uword> indices_int, T amplitude_val);
  private:
	uword n_sites;
	urowvec perm;
	T chi;
	string name;
  };

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
  inline
  uword Symmetry<T>::getNumSites() const{
	return this->n_sites;
  }

  template<typename T>
  inline
  string Symmetry<T>::getName() const{
	return this->name;
  }

  template<typename T>
  inline
  void Symmetry<T>::setName(string new_name){
	this->name = new_name;
  }

  template<typename T>
  inline
  urowvec Symmetry<T>::getPerm() const{
	return this->perm;
  }

  template<typename T>
  inline
  void Symmetry<T>::setPerm(const urowvec& other_perm){
	this->perm = other_perm;
  }

  template<typename T>
  inline
  T Symmetry<T>::getChi() const{
	return this->chi;
  }

  template<typename T>
  void Symmetry<T>::setChi(T new_chi){
	this->chi = new_chi;
  }

  template<typename T>
  inline
  const uword& Symmetry<T>::operator[](uword pos) const{
	return this->perm[pos];
  }

  template<typename T>
  inline
  uword& Symmetry<T>::operator[](uword pos){
	return this->perm[pos];
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
  // ***         Functions        *** //

  template<typename T>
  Symmetry<T> operator*(const Symmetry<T>& sym1, const Symmetry<T>& sym2){
	uword n_sites = sym1.getNumSites();
	if (n_sites != sym2.getNumSites()){
	  throw logic_error("In Symmetry::operator+. Size problem.");
	}
	Symmetry<T> sym_out(n_sites);
	sym_out.setName(sym1.getName() + "*" + sym2.getName());
	sym_out.setChi(sym1.getChi() * sym2.getChi());
	for (uword i=0; i<n_sites; i++)
	  sym_out[i] = sym1[sym2[i]];
	return sym_out;
  }

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed


#endif //* SYMMETRY_H
