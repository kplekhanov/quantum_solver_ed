#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "headers.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  Symmetry class          *** //

  template<typename T=cx_double>
  class Symmetry{
  public:
	Symmetry();
	Symmetry(uword n_sites);
	inline
	uword getNumSites() const {return this->n_sites;}
	inline
	string getName() const {return this->name;}
	void setName(string new_name);
	inline
	urowvec getPerm() const {return this->perm;}
	void setPerm(const urowvec& other_perm);
	inline
	T getChi() const {return this->chi;}
	void setChi(T new_chi);
	inline
	const uword& operator[](uword pos) const {return this->perm[pos];}
	inline
	uword& operator[](uword pos) {return this->perm[pos];}
	urowvec apply(const urowvec& conf_vec_in) const;
	void print(ostream& os) const;
	void fill(vector<uword> indices_int, T amplitude_val);
  private:
	uword n_sites;
	urowvec perm;
	T chi;
	string name;
  };
  // *** ------------------------ *** //

  // *** ------------------------ *** //
  // ***         Functions        *** //

  template<typename T>
  ostream& operator<<(ostream& os, const Symmetry<T>& sym){
	sym.print(os);
	return os;
  }
  
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
