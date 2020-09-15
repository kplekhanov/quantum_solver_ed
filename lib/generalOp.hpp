#ifndef GENERAL_OP_H
#define GENERAL_OP_H

#include "hilbertSym.hpp"
#include "elementaryOp.hpp"
#include "waveFon.hpp"
#include "word.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // *** GeneralOp acts on states *** //

  // *** This is a class for general operator acting on wavefunctions
  // *** Such operators can be decomposed into elementary operators acting on confints
  // *** typical example is the Hamiltonian
  // *** has methods to get a sparse/dense matrix representation or doing Lanczos on fly
  // *** can be read from file; in this case, the template S is elemetaryOp
  // *** each file has to correspond to a given elementaryOp
  template<typename T=double>
  class GeneralOp{
	typedef vector<DiagOp*> VectorDiagOpPtr;
	typedef vector<OffDiagOp<T>*> VectorOffDiagOpPtr;
  public:
	GeneralOp(const HilbertBones* hil_ptr);
	void append(DiagOp* dop_ptr, double coef);
	void append(DiagOp* dop_ptr, cx_double coef);
	void append(OffDiagOp<T>* oop_ptr, T coef);
	WaveFon<T> apply(const WaveFon<T>& phi_in) const;
	Col<double> doLanczosOnFly(uword N, double err=1e-10);
	template<typename Q>
	void createMatrix(Q* matrix_ptr);
	void print() const;
	template <typename S>
	void readFromFile(string file_name, char delimiter=';');
  private:
	const HilbertBones* hil_ptr;
	VectorDiagOpPtr dop_ptr_vec;
	VectorOffDiagOpPtr oop_ptr_vec;
	vector<double> dop_coef_vec;
	vector<T> oop_coef_vec;
	WaveFon<T> applyNonSym(const WaveFon<T>& phi_in) const;
	WaveFon<T> applySym(const WaveFon<T>& phi_in) const;
	template<typename Q>
	void createMatrixNonSym(Q* matrix_ptr);
	template<typename Q>
	void createMatrixSym(Q* matrix_ptr);
  };

  // *** ---------------- *** //

} //* namespace quantum_solver_ed


#endif //* GENERAL_OP_H
