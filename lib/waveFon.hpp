#ifndef WAVEFON_H
#define WAVEFON_H

#include "headers.hpp"
#include "hilbertBones.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  Wave-function class     *** //

  // *** basically a wrapper around armadillo vectors
  // *** additional stuff -- checking of the pointers to the Hilbert space,
  // *** as well as a proper scalar product
  template<typename T=double>
  class WaveFon{
  public:
	template<typename L> friend inline
	WaveFon<L> operator+(const WaveFon<L>& phi_in1, const WaveFon<L>& phi_in2);
	template<typename L> friend inline
	WaveFon<L> operator-(const WaveFon<L>& phi_in1, const WaveFon<L>& phi_in2);
	template<typename L> friend inline
	L operator*(const WaveFon<L>& phi_in1, const WaveFon<L>& phi_in2);
	template<typename L> friend inline
	ostream& operator<<(ostream& os, const WaveFon<L>& phi);
	template<typename L, typename M> friend inline
	WaveFon<L> operator*(const WaveFon<L>& phi_in, const M coef);
	template<typename L, typename M> friend inline
	WaveFon<L> operator*(const M coef, const WaveFon<L>& phi_in);
	template<typename L, typename M> friend inline
	WaveFon<L> operator/(const WaveFon<L>& phi_in, const M coef);
  
	WaveFon(const HilbertBones* hil_ptr);
	inline
	const HilbertBones* getHilPtr() const {return this->hil_ptr;}
	void randomise();
	void normalise();
	double getNorm() const;
	const MyVec<T>& getVec() const;
	void setVec(const MyVec<T>& new_vec);
	inline
	T& operator[](uword pos) {return this->vec[pos];}
	inline
	const T& operator[](uword pos) const {return this->vec[pos];};
	WaveFon<T>& operator=(const WaveFon<T>& phi);
  private:
	const HilbertBones* hil_ptr;
	MyVec<T> vec;
  };

  // *** ------------------------ *** //

  // *** ------------------------ *** //
  // ***         Functions        *** //

  // *** unfortunately, because functions are template, they have to be inside hpp
  // *** otherwise one needs to do a lot of ugly instantiation in the cpp
  // *** moreover, becase it's a header, all the functions have to be inline
  // *** if you know how to do it better, PLEASE LET ME KNOW!!!!!
  template<typename T> inline
  WaveFon<T> operator+(const WaveFon<T>& phi_in1, const WaveFon<T>& phi_in2){
	if (phi_in1.hil_ptr != phi_in2.hil_ptr){
	  throw logic_error("In WaveFon::operator+. Pointer problem.");
	}
	WaveFon<T> phi_out(phi_in1.hil_ptr);
	phi_out.vec = phi_in1.vec + phi_in2.vec;
	return phi_out;
  }

  template<typename T> inline
  WaveFon<T> operator-(const WaveFon<T>& phi_in1, const WaveFon<T>& phi_in2){
	if (phi_in1.hil_ptr != phi_in2.hil_ptr){
	  throw logic_error("In WaveFon::operator-. Pointer problem.");
	}
	WaveFon<T> phi_out(phi_in1.hil_ptr);
	phi_out.vec = phi_in1.vec - phi_in2.vec;
	return phi_out;
  }

  template<typename T> inline
  T operator*(const WaveFon<T>& phi_in1, const WaveFon<T>& phi_in2){
	if (phi_in1.hil_ptr != phi_in2.hil_ptr){
	  throw logic_error("In WaveFon::operator*. Pointer problem.");
	}
	return cdot(phi_in1.vec, phi_in2.vec);
  }

  template<typename T> inline
  ostream& operator<<(ostream& os, const WaveFon<T>& phi){
	os << phi.vec;
	return os;
  }

  template<typename L, typename M> inline
  WaveFon<L> operator*(const WaveFon<L>& phi_in, const M coef){
	if (is_same<L, double>::value and is_same<M, cx_double>::value){
	  throw logic_error("In WaveFon::operator*."
						"Can't multiply real WaveFon and a complex number.");
	}
	WaveFon<L> phi_out(phi_in.hil_ptr);
	phi_out.vec = phi_in.vec * coef;
	return phi_out;
  }

  template<typename L, typename M> inline
  WaveFon<L> operator*(const M coef, const WaveFon<L>& phi_in){
	if (is_same<L, double>::value and is_same<M, cx_double>::value){
	  throw logic_error("In WaveFon::operator*."
						"Can't multiply real WaveFon and a complex number.");
	}
	WaveFon<L> phi_out(phi_in.hil_ptr);
	phi_out.vec = phi_in.vec * coef;
	return phi_out;
  }

  template<typename L, typename M> inline
  WaveFon<L> operator/(const WaveFon<L>& phi_in, const M coef){
	if (is_same<L, double>::value and is_same<M, cx_double>::value){
	  throw logic_error("In WaveFon::operator/."
						"Can't multiply real WaveFon and a complex number.");
	}
	WaveFon<L> phi_out(phi_in.hil_ptr);
	phi_out.vec = phi_in.vec / coef;
	return phi_out;
  }

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed


#endif //* WAVEFON_H
