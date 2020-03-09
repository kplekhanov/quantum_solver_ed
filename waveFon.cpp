#ifndef WAVEFON_H
#define WAVEFON_H

#include "headers.hpp"
#include "hilbertBones.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  Wave-function class     *** //

  template<typename T=double>
  class WaveFon{
  public:
	template<typename L> friend
	WaveFon<L> operator+(const WaveFon<L>& phi_in1, const WaveFon<L>& phi_in2);
	template<typename L> friend
	WaveFon<L> operator-(const WaveFon<L>& phi_in1, const WaveFon<L>& phi_in2);
	template<typename L> friend
	L operator*(const WaveFon<L>& phi_in1, const WaveFon<L>& phi_in2);
	template<typename L> friend
	ostream& operator<<(ostream& os, const WaveFon<L>& phi);
	template<typename L, typename M> friend
	WaveFon<L> operator*(const WaveFon<L>& phi_in, const M coef);
	template<typename L, typename M> friend
	WaveFon<L> operator*(const M coef, const WaveFon<L>& phi_in);
	template<typename L, typename M> friend
	WaveFon<L> operator/(const WaveFon<L>& phi_in, const M coef);
  
	WaveFon(const HilbertBones* hil_ptr);
	const HilbertBones* getHilPtr() const;
	void randomise();
	void normalise();
	double getNorm() const;
	const MyVec<T>& getVec() const;
	void setVec(const MyVec<T>& new_vec);
	T& operator[](uword pos);
	const T& operator[](uword pos) const;
	WaveFon<T>& operator=(const WaveFon<T>& phi);
  private:
	const HilbertBones* hil_ptr;
	MyVec<T> vec;
  };

  template<typename T>
  WaveFon<T>::WaveFon(const HilbertBones* hil_ptr){
	this->hil_ptr = hil_ptr;
	this->vec = MyVec<T>(hil_ptr->getNumStates(), fill::zeros);
  }

  template<typename T> inline
  const HilbertBones* WaveFon<T>::getHilPtr() const{
	return this->hil_ptr;
  }

  template<typename T>
  void WaveFon<T>::randomise(){
	this->vec.randu();
	this->normalise();
  }

  template<typename T> inline
  void WaveFon<T>::normalise(){
	this->vec = arma::normalise(this->vec);
  }

  template<typename T> inline
  double WaveFon<T>::getNorm() const{
	return norm(this->vec);
  }

  template<typename T>
  const MyVec<T>& WaveFon<T>::getVec() const{
	return this->vec;
  }

  template<typename T>
  void WaveFon<T>::setVec(const MyVec<T>& new_vec){
	this->vec = new_vec;
  }

  template<typename T> inline
  T& WaveFon<T>::operator[](uword pos){
	return this->vec[pos];
  }

  template<typename T> inline
  const T& WaveFon<T>::operator[](uword pos) const{
	return this->vec[pos];
  }

  template<typename T>
  WaveFon<T>& WaveFon<T>::operator=(const WaveFon<T>& phi){
	if (phi.hil_ptr != this->hil_ptr){
	  throw logic_error("In WaveFon::operator=. Pointer problem.");
	}
	this->vec = phi.vec;
	return *this;
  }

  // *** ------------------------ *** //


  // *** ------------------------ *** //
  // ***         Functions        *** //

  template<typename T>
  WaveFon<T> operator+(const WaveFon<T>& phi_in1, const WaveFon<T>& phi_in2){
	if (phi_in1.hil_ptr != phi_in2.hil_ptr){
	  throw logic_error("In WaveFon::operator+. Pointer problem.");
	}
	WaveFon<T> phi_out(phi_in1.hil_ptr);
	phi_out.vec = phi_in1.vec + phi_in2.vec;
	return phi_out;
  }

  template<typename T>
  WaveFon<T> operator-(const WaveFon<T>& phi_in1, const WaveFon<T>& phi_in2){
	if (phi_in1.hil_ptr != phi_in2.hil_ptr){
	  throw logic_error("In WaveFon::operator-. Pointer problem.");
	}
	WaveFon<T> phi_out(phi_in1.hil_ptr);
	phi_out.vec = phi_in1.vec - phi_in2.vec;
	return phi_out;
  }

  template<typename T>
  T operator*(const WaveFon<T>& phi_in1, const WaveFon<T>& phi_in2){
	if (phi_in1.hil_ptr != phi_in2.hil_ptr){
	  throw logic_error("In WaveFon::operator*. Pointer problem.");
	}
	return cdot(phi_in1.vec, phi_in2.vec);
  }

  template<typename T>
  ostream& operator<<(ostream& os, const WaveFon<T>& phi){
	os << phi.vec;
	return os;
  }

  template<typename L, typename M>
  WaveFon<L> operator*(const WaveFon<L>& phi_in, const M coef){
	if (is_same<L, double>::value and is_same<M, cx_double>::value){
	  throw logic_error("In WaveFon::operator*."
							 "Can't multiply real WaveFon and a complex number.");
	}
	WaveFon<L> phi_out(phi_in.hil_ptr);
	phi_out.vec = phi_in.vec * coef;
	return phi_out;
  }

  template<typename L, typename M>
  WaveFon<L> operator*(const M coef, const WaveFon<L>& phi_in){
	if (is_same<L, double>::value and is_same<M, cx_double>::value){
	  throw logic_error("In WaveFon::operator*."
							 "Can't multiply real WaveFon and a complex number.");
	}
	WaveFon<L> phi_out(phi_in.hil_ptr);
	phi_out.vec = phi_in.vec * coef;
	return phi_out;
  }

  template<typename L, typename M>
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
