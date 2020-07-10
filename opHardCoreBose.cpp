#ifndef OP_HARD_CORE_BOSE_H
#define OP_HARD_CORE_BOSE_H

#include "elementaryOp.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** chemical potential
  class Ni: public DiagOp{
  public:
	Ni(){};
	Ni(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
  };

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

  inline
  double Ni::apply(const uword& conf_in) const{
	return double((conf_in >> i) & 1);
  }

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** hopping operators
  template<typename T=double>
  class BdagiBj: public OffDiagOp<T>{
  public:
	BdagiBj(){};
	BdagiBj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

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
	if ( (((conf_in >> i) & 1) == 0) & (((conf_in >> j) & 1) == 1) )
	  return PairUwordT<T>((conf_in & ~(1 << j)) | (1 << i), double(1));
	else
	  return PairUwordT<T>(0 , 0);
  }

  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_HARD_CORE_BOSE_H
