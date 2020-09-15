#include "../lib/waveFon.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // ***  Wave-function class     *** //

  template<typename T>
  WaveFon<T>::WaveFon(const HilbertBones* hil_ptr){
	this->hil_ptr = hil_ptr;
	this->vec = MyVec<T>(hil_ptr->getNumStates(), fill::zeros);
  }

  template<typename T>
  void WaveFon<T>::randomise(){
	this->vec.randu();
	this->normalise();
  }

  template<typename T>
  void WaveFon<T>::normalise(){
	this->vec = arma::normalise(this->vec);
  }

  template<typename T>
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
  // *** instatation of template  *** //
  
  template class WaveFon<double>;
  template class WaveFon<cx_double>;

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
