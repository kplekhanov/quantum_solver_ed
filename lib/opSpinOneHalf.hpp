#ifndef OP_SPIN_ONE_HALF_H
#define OP_SPIN_ONE_HALF_H

#include "elementaryOp.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** Zeeman along the z-axis
  class Szi_2: public DiagOp{
  public:
	Szi_2(){}
	Szi_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	inline
	double apply(const uword& conf_in) const{
	  return double((conf_in >> i) & 1) - 0.5;
	}
  private:
	uword i;
  };

  // *** Exchange along the z-axis
  class SziSzj_2: public DiagOp{
  public:
	SziSzj_2(){}
	SziSzj_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	inline
	double apply(const uword& conf_in) const{
	  return (double((conf_in >> i) & 1) - 0.5) * (double((conf_in >> j) & 1) - 0.5);
	}
  private:
	uword i;
	uword j;
  };

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** swap gates
  template<typename T=double>
  class SpiSmj_2: public OffDiagOp<T>{
  public:
	SpiSmj_2(){}
	SpiSmj_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  // *** in-plane DM interactions -- Sz Sp
  template<typename T=double>
  class SziSpj_2: public OffDiagOp<T>{
  public:
	SziSpj_2(){}
	SziSpj_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  
  // *** in-plane DM interactions -- Sz Sm
  template<typename T=double>
  class SziSmj_2: public OffDiagOp<T>{
  public:
	SziSmj_2(){}
	SziSmj_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  
  // *** in-plane DM interactions -- Sz Sx
  template<typename T=double>
  class SziSxj_2: public OffDiagOp<T>{
  public:
	SziSxj_2(){}
	SziSxj_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  // *** in-plane DM interactions -- Sz Sy
  template<typename T>
  class SziSyj_2: public OffDiagOp<T>{
  public:
	SziSyj_2(){}
	SziSyj_2(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

 
  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_SPIN_ONE_HALF_H
