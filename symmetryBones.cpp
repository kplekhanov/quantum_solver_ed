#ifndef SYMMETRY_BONES_H
#define SYMMETRY_BONES_H

#include "headers.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  SymmetryBones class     *** //

  template<typename T=cx_double>
  class SymmetryBones{
  public:
	SymmetryBones();
	virtual string getName() const = 0;
	virtual void setName(string new_name) = 0;
	virtual T getChi() const = 0;
	virtual void setChi(T new_chi) = 0;
	virtual urowvec apply(const urowvec& conf_vec_in) const = 0;
	virtual void print(ostream& os) const = 0;
  private:
	T chi;
	string name;
  };

  template<typename T>
  SymmetryBones<T>::SymmetryBones(){}


  // *** ------------------------ *** //


  // *** ------------------------ *** //
  // ***         Functions        *** //

  template<typename T>
  ostream& operator<<(ostream& os, const SymmetryBones<T>& sym){
	sym.print(os);
	return os;
  }

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed


#endif //* SYMMETRY_BONES_H
