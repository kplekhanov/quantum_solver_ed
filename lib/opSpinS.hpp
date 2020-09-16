#ifndef OP_SPIN_H
#define OP_SPIN_H

#include "elementaryOp.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** Zeeman along the z-axis
  class Szi: public DiagOp{
  public:
	Szi(){}
	Szi(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
  };

  // *** Exchange along the z-axis
  class SziSzj: public DiagOp{
  public:
	SziSzj(){}
	SziSzj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
	uword j;
  };

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** swap gates
  template<typename T=double>
  class SpiSmj: public OffDiagOp<T>{
  public:
	SpiSmj(){}
	SpiSmj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  // *** in-plane DM interactions -- Sz Sp
  template<typename T=double>
  class SziSpj: public OffDiagOp<T>{
  public:
	SziSpj(){}
	SziSpj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };
  
  // *** in-plane DM interactions -- Sz Sm
  template<typename T=double>
  class SziSmj: public OffDiagOp<T>{
  public:
	SziSmj(){}
	SziSmj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };
  
  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_SPIN_H
