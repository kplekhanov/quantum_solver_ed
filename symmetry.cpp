#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "symmetryBones.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  Symmetry class          *** //

  class Symmetry: public SymmetryBones{
  public:
	Symmetry(uword n_sites);
	uword getNumSites() const;
	string getName() const;
	void setName(string new_name);
	urowvec getPerm() const;
	void setPerm(const urowvec& other_perm);
	cx_double getChi() const;
	void setChi(cx_double new_chi);
	const uword& operator[](uword pos) const;
	uword& operator[](uword pos);
	urowvec apply(const urowvec& conf_vec_in) const;
	void print(ostream& os) const;
  private:
	uword n_sites;
	urowvec perm;
	cx_double chi;
	string name;
  };

  Symmetry::Symmetry(uword n_sites){
	this->n_sites = n_sites;
	perm = urowvec(n_sites);
	for (uword i=0; i<n_sites; i++)
	  this->perm[i] = i;
	this->chi = 1;
	this->name = "id" + to_string(this->n_sites);
  }

  inline
  uword Symmetry::getNumSites() const{
	return this->n_sites;
  }

  inline
  string Symmetry::getName() const{
	return this->name;
  }

  inline
  void Symmetry::setName(string new_name){
	this->name = new_name;
  }

  inline
  urowvec Symmetry::getPerm() const{
	return this->perm;
  }

  inline
  void Symmetry::setPerm(const urowvec& other_perm){
	this->perm = other_perm;
  }

  inline
  cx_double Symmetry::getChi() const{
	return this->chi;
  }

  void Symmetry::setChi(cx_double new_chi){
	this->chi = new_chi;
  }

  inline
  const uword& Symmetry::operator[](uword pos) const{
	return this->perm[pos];
  }

  inline
  uword& Symmetry::operator[](uword pos){
	return this->perm[pos];
  }

  urowvec Symmetry::apply(const urowvec& conf_vec_in) const{
	uword n_sites = conf_vec_in.n_elem;
	if (n_sites != this->n_sites){
	  throw logic_error("In Symmetry::apply. Size problem.");
	}
	urowvec conf_vec_out(n_sites);
	for (uword i=0; i<n_sites; i++)
	  conf_vec_out[this->perm[i]] = conf_vec_in[i];
	return conf_vec_out;
  }

  void Symmetry::print(ostream& os) const{
	os << this->getName() << "|";
	for (uword i=0; i<this->getNumSites(); i++)
	  os << " " << (*this)[i];
	os << " | " << this->getChi() <<  " |";
  }

  // *** ------------------------ *** //


  // *** ------------------------ *** //
  // ***         Functions        *** //

  Symmetry operator*(const Symmetry& sym1, const Symmetry& sym2){
	uword n_sites = sym1.getNumSites();
	if (n_sites != sym2.getNumSites()){
	  throw logic_error("In Symmetry::operator+. Size problem.");
	}
	Symmetry sym_out(n_sites);
	sym_out.setName(sym1.getName() + "*" + sym2.getName());
	sym_out.setChi(sym1.getChi() * sym2.getChi());
	for (uword i=0; i<n_sites; i++)
	  sym_out[i] = sym1[sym2[i]];
	return sym_out;
  }

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed


#endif //* SYMMETRY_H
