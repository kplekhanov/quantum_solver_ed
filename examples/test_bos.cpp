#include "../lib/hilbert.hpp"
#include "../lib/generalOp.hpp"
#include "../lib/opBoseHubbard.hpp"


using namespace quantum_solver_ed;
using namespace arma;
using namespace std;

int main(int argc, char *argv[]){
  // variables
  const double pi = std::acos(-1);
  uword N = 10;
  uword Qn = 5;
  wall_clock timer;
  urowvec localSizes = urowvec(N, fill::ones)*2;
  double k = 0; // wave-vector number

  // creating the HS
  cout << "Testing Hil" << endl;
  timer.tic();
  Hilbert hil(localSizes);
  hil.createWithFixedQn(Qn);
  const HilbertBones* hil_ptr = hil.getHilPtr();
  cout << "Creating Hil took " << timer.toc() << endl;
  //hil.print();
  cout << endl;

  // creating the Hamiltonian
  GeneralOp<double> ham(hil_ptr);
  for (uword i=0; i<N; i++){
	ham.append(new Ni(hil_ptr, {i}), 0.123);
	ham.append(new BdagiBj<double>(hil_ptr, {i, (i+1)%N}), -1.0);
	ham.append(new BdagiBj<double>(hil_ptr, {(i+1)%N, i}), -1.0);
  }
  ham.print();
  cout << endl;
  
  // testing Lanczos
  cout << "Testing Lanczos" << endl;
  cout.precision(8);
  timer.tic();
  ham.doLanczosOnFly(100, 1e-14);
  cout << "Testing Lanczos took " << timer.toc() << endl << endl;

  // testing dense matrix
  cout << "Creating dense matrix" << endl;
  cout.precision(8);
  timer.tic();
  Mat<double> ham_mat_2;
  ham.createMatrix(&ham_mat_2);
  cout << "Creating dense matrix took " << timer.toc() << endl << endl;
  //cout << ham_mat_2;

  // testing sparse matrix
  cout << "Creating sparse matrix" << endl;
  cout.precision(8);
  timer.tic();
  SpMat<double> ham_sp_mat_2;
  ham.createMatrix(&ham_sp_mat_2);
  cout << "Creating sparse matrix took " << timer.toc() << endl << endl;

  // getting eigenvalues
  vec eigval = eigs_sym(ham_sp_mat_2, 10);
  cout << eigval;

  return 0;
}
