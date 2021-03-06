#include "../lib/generalOp.hpp"

namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ------------------------ *** //
  // *** GeneralOp acts on states *** //

  template<typename T>
  GeneralOp<T>::GeneralOp(const HilbertBones* hil_ptr){
	this->hil_ptr = hil_ptr;
  }

  template<typename T>
  void GeneralOp<T>::append(DiagOp* dop_ptr, double coef){
	if (dop_ptr->getHilPtr() != this->hil_ptr){
	  throw logic_error( "In GeneralOp::append. Pointer problem." );
	}
	this->dop_ptr_vec.push_back(dop_ptr);
	this->dop_coef_vec.push_back(coef);
  }

  template<typename T>
  void GeneralOp<T>::append(DiagOp* dop_ptr, cx_double coef){
	// if coef is complex, take the real part only
	if (dop_ptr->getHilPtr() != this->hil_ptr){
	  throw logic_error( "In GeneralOp::append. Pointer problem." );
	}
	this->dop_ptr_vec.push_back(dop_ptr);
	this->dop_coef_vec.push_back(coef.real());
  }

  template<typename T>
  void GeneralOp<T>::append(OffDiagOp<T>* oop_ptr, T coef){
	if (oop_ptr->getHilPtr() != this->hil_ptr){
	  throw logic_error( "In GeneralOp::append. Pointer problem." );
	}
	this->oop_ptr_vec.push_back(oop_ptr);
	this->oop_coef_vec.push_back(coef);
  }

  template<typename T>
  WaveFon<T> GeneralOp<T>::apply(const WaveFon<T>& phi_in) const{
	if (this->hil_ptr->checkSymmetric() == true)
	  return GeneralOp<T>::applySym(phi_in);
	else
	  return GeneralOp<T>::applyNonSym(phi_in);
  }

  template<typename T>
  WaveFon<T> GeneralOp<T>::applyNonSym(const WaveFon<T>& phi_in) const{
	if (phi_in.getHilPtr() != this->hil_ptr){
	  throw logic_error( "In GeneralOp::apply. Pointer problem." );
	}
	const HilbertBones* hil_ptr = this->hil_ptr;
	uword num_iters_per_thread = hil_ptr->getNumStates() / OMP_NUM_THREADS;
	WaveFon<T> phi_out(hil_ptr);
	uword dop_len = this->dop_ptr_vec.size();
#pragma omp parallel for schedule(static, num_iters_per_thread)
	for (uword i=0; i<hil_ptr->getNumStates(); ++i){
	  T coef_old = phi_in[i];
	  if (abs(coef_old) != 0){
		uword conf_old = hil_ptr->getConfInt(i);
		T coef_to_add = 0;
		for (uword j=0; j<dop_len; ++j)
		  coef_to_add += this->dop_coef_vec[j]
			* this->dop_ptr_vec[j]->apply(conf_old);
		phi_out[i] += coef_to_add * coef_old;
	  }
	}
	uword oop_len = this->oop_ptr_vec.size();
#pragma omp parallel for schedule(static, num_iters_per_thread)
	for (uword index_new=0; index_new<hil_ptr->getNumStates(); ++index_new){
	  uword conf_new = hil_ptr->getConfInt(index_new);
	  T coef_to_add = 0;
	  for (uword j=0; j<oop_len; ++j){
		PairUwordT<T> conf_op_pair = this->oop_ptr_vec[j]->apply(conf_new);
		if (abs(conf_op_pair.second) != 0){
		  uword index_old = hil_ptr->getConfPair(conf_op_pair.first);
		  coef_to_add += phi_in[index_old] * conf_op_pair.second
			* this->oop_coef_vec[j];
		}
	  }
	  phi_out[index_new] += coef_to_add;
	}
	return phi_out;
  }
  
  template<typename T>
  WaveFon<T> GeneralOp<T>::applySym(const WaveFon<T>& phi_in) const{
	if (phi_in.getHilPtr() != this->hil_ptr)
	  throw logic_error( "In GeneralOp::apply. Pointer problem." );
	const HilbertSym<T>* hil_ptr = dynamic_cast<const HilbertSym<T>*>(this->hil_ptr);
	if (! hil_ptr)
	  throw logic_error( "In GeneralOp::apply. HilbertSym problem." );
	uword num_iters_per_thread = hil_ptr->getNumStates() / OMP_NUM_THREADS;
	WaveFon<T> phi_out(hil_ptr);
	uword dop_len = this->dop_ptr_vec.size();
#pragma omp parallel for schedule(static, num_iters_per_thread)
	for (uword i=0; i<hil_ptr->getNumStates(); ++i){
	  T coef_old = phi_in[i];
	  if (abs(coef_old) != 0){
		uword conf_old = Word::getRep(hil_ptr->getConfInt(i));
		T coef_to_add = 0;
		for (uword j=0; j<dop_len; ++j)
		  coef_to_add += this->dop_coef_vec[j]
			* this->dop_ptr_vec[j]->apply(conf_old);
		phi_out[i] += coef_to_add * coef_old;
	  }
	}
	uword oop_len = this->oop_ptr_vec.size();
#pragma omp parallel for schedule(static, num_iters_per_thread)
	for (uword index_new=0; index_new<hil_ptr->getNumStates(); ++index_new){
	  Word w(hil_ptr->getConfInt(index_new));
	  uword w_rep = w.getRep();
	  double w_deg_sqrt = sqrt(double(w.getDeg()));
	  T coef_to_add = 0;
	  for (uword j=0; j<oop_len; ++j){
		pair<uword,T> conf_op_pair = this->oop_ptr_vec[j]->apply(w_rep);
		if (abs(conf_op_pair.second) != 0){
		  // getting first = rep index for a state (from above), and second = coef
		  PairUwordT<T> sympair = hil_ptr->getSymPair(conf_op_pair.first);
		  coef_to_add += this->oop_coef_vec[j] * phi_in[sympair.first]
			* conf_op_pair.second * sympair.second / w_deg_sqrt;
		}
	  }
	  phi_out[index_new] += coef_to_add;
	}
	return phi_out;
  }
  
  template<typename T>
  template<typename Q>
  void GeneralOp<T>::createMatrix(Q* matrix_ptr){
	if (is_same<Q, Mat<T>>::value){
	  Mat<T> mat(this->hil_ptr->getNumStates(), this->hil_ptr->getNumStates(), fill::zeros);
	  *matrix_ptr = mat;
	}
	else if (is_same<Q, SpMat<T>>::value){
	  SpMat<T> mat(this->hil_ptr->getNumStates(), this->hil_ptr->getNumStates());
	  *matrix_ptr = mat;
	}
	else
	  throw logic_error( "In GeneralOp::createMatrx. Matrix type is unknown." );
  	if (this->hil_ptr->checkSymmetric() == true)
	  this->createMatrixSym<Q>(matrix_ptr);
	else
	  this->createMatrixNonSym<Q>(matrix_ptr);
  }

  template<typename T>
  template<typename Q>
  void GeneralOp<T>::createMatrixNonSym(Q* matrix_ptr){
	const HilbertBones* hil_ptr = this->hil_ptr;
	uword dop_len = this->dop_ptr_vec.size();
	uword oop_len = this->oop_ptr_vec.size();
	for (uword j_conf=0; j_conf<hil_ptr->getNumStates(); ++j_conf){
	  uword conf_old = hil_ptr->getConfInt(j_conf);
	  // diagonal elements
	  for (uword j_dop=0; j_dop<dop_len; ++j_dop)
		(*matrix_ptr)(j_conf,j_conf) += this->dop_coef_vec[j_dop]
		  * this->dop_ptr_vec[j_dop]->apply(conf_old);
	  // off-diagonal elements
	  for (uword j_oop=0; j_oop<oop_len; ++j_oop){
		PairUwordT<T> conf_op_pair = this->oop_ptr_vec[j_oop]->apply(conf_old);
		if (abs(conf_op_pair.second) != 0){
		  uword k_conf = hil_ptr->getConfPair(conf_op_pair.first);
		  (*matrix_ptr)(j_conf,k_conf) += this->oop_coef_vec[j_oop] * conf_op_pair.second;
		}
	  }
	}
	if (!(*matrix_ptr).is_hermitian(1e-10))
	  throw logic_error( "In GeneralOp::createMatrixNonSym. Matrix non-hermitian." );
  }
  
  template<typename T>
  template<typename Q>
  void GeneralOp<T>::createMatrixSym(Q* matrix_ptr){
	const HilbertSym<T>* hil_ptr = dynamic_cast<const HilbertSym<T>*>(this->hil_ptr);
	if (!hil_ptr)
	  throw logic_error( "In GeneralOp::apply. HilbertSym problem." );
	uword dop_len = this->dop_ptr_vec.size();
	uword oop_len = this->oop_ptr_vec.size();
	for (uword j_conf=0; j_conf<hil_ptr->getNumStates(); ++j_conf){
	  Word w(hil_ptr->getConfInt(j_conf));
	  uword w_rep = w.getRep();
	  double w_deg_sqrt = sqrt(double(w.getDeg()));
	  // diagonal elements
	  for (uword j_dop=0; j_dop<dop_len; ++j_dop)
		(*matrix_ptr)(j_conf,j_conf) += this->dop_coef_vec[j_dop]
		  * this->dop_ptr_vec[j_dop]->apply(w_rep);
	  // off-diagonal elements
	  for (uword j_oop=0; j_oop<oop_len; ++j_oop){
		pair<uword,T> conf_op_pair = this->oop_ptr_vec[j_oop]->apply(w_rep);
		if (abs(conf_op_pair.second) != 0){
		  // getting first = rep index for a state (from above), and second = coef
		  PairUwordT<T> sympair = hil_ptr->getSymPair(conf_op_pair.first);
		  (*matrix_ptr)(j_conf,sympair.first) += this->oop_coef_vec[j_oop]
			* conf_op_pair.second * sympair.second / w_deg_sqrt;
		}
	  }
	}
	if (!(*matrix_ptr).is_hermitian(1e-10))
	  throw logic_error( "In GeneralOp::createMatrixSym. Matrix non-hermitian." );
  }

  template<typename T>
  Col<double> GeneralOp<T>::doLanczosOnFly(uword N, double err){
	Mat<double> tridiag = Mat<double>(N, N, fill::zeros);
	Mat<double> tridiag_temp;
	Col<double> eigvals;
	WaveFon<T> v_0(this->hil_ptr); // initialized to zero
	WaveFon<T> v_1(this->hil_ptr);
	WaveFon<T> w_1(this->hil_ptr);

	wall_clock timer; //
  
	v_1.randomise();
	double b_1 = 0;
	double a_1 = 0;
	double e_min = 1e16;
	uword i = 0;
	while (i < N){
	  timer.tic(); //
	  w_1 = this->apply(v_1) - b_1 * v_0;
	  a_1 = real(w_1 * v_1);
	  w_1 = w_1 - a_1 * v_1;
	  b_1 = w_1.getNorm();
	  v_0 = v_1;
	  v_1 = w_1 / b_1;
    
	  tridiag(i,i) = a_1;
	  if (i < N - 1){
		tridiag(i,i+1) = b_1;
		tridiag(i+1,i) = b_1;
	  }
	  tridiag_temp = tridiag.submat(0, 0, i, i);
	  eigvals = eig_sym(tridiag_temp);

	  if (abs(e_min - eigvals[0]) < err)
		break;
	  ++i;      
	  e_min = eigvals[0];
	}
	return eigvals;
  }

  template<typename T>
  void GeneralOp<T>::print() const{
	cout << "Printing GeneralOp" << endl;
	uword dop_len = this->dop_ptr_vec.size();
	cout << "Diagonal operators: ";
	for (uword i=0; i<dop_len; i++){
	  cout << "(" << this->dop_coef_vec[i]  << ") * " << *(this->dop_ptr_vec[i]);
	  if (i != dop_len-1) cout << " + ";
	}
	cout << endl;
	uword oop_len = this->oop_ptr_vec.size();
	cout << "Off-diagonal operators: ";
	for (uword i=0; i<oop_len; i++){
	  cout << "(" << this->oop_coef_vec[i]  << ") * " << *(this->oop_ptr_vec[i]);
	  if (i != oop_len-1) cout << " + ";
	}
	cout << endl;
  }

  // *** ------------------------ *** //

  // *** ------------------------ *** //
  // *** instatation of template  *** //
  
  template class GeneralOp<double>;
  template class GeneralOp<cx_double>;

  template void GeneralOp<double>::createMatrix(Mat<double>* matrix_ptr);
  template void GeneralOp<cx_double>::createMatrix(Mat<cx_double>* matrix_ptr);
  template void GeneralOp<double>::createMatrix(SpMat<double>* matrix_ptr);
  template void GeneralOp<cx_double>::createMatrix(SpMat<cx_double>* matrix_ptr);

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
