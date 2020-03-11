#include "../hilbertSym.cpp"
#include "../waveFon.cpp"
#include "../symmetry.cpp"
#include "../generalOp.cpp"
#include "../opBoseHubbard.cpp"


using namespace quantum_solver_ed;
using namespace arma;
using namespace std;


int main(int argc, char *argv[]){
  // variables
  const double pi = std::acos(-1);
  const cx_double cx_i(0, 1);
  uword N = 24;
  uword Qn = 12;
  wall_clock timer;
  urowvec localSizes = urowvec(N, fill::ones)*2;
  double k = 0; // wave-vector number

  // creating the HS
  cout << "Testing HilSym" << endl;
  timer.tic();
  HilbertSym<cx_double> hil(localSizes);
  for (uword i=0; i<N; i++){
	Symmetry<cx_double>* sym_ptr = new Symmetry<cx_double>(N);
	sym_ptr->setName("t_"+std::to_string(i)+".k_"+std::to_string(k));
	sym_ptr->setChi(std::exp(2 * pi * i * k * cx_i / double(N)));
	for (uword j=0; j<N; j++)
	  (*sym_ptr)[j] = (j + i) % N;
	hil.addSym(sym_ptr);
  }
  hil.createWithFixedQn(Qn);
  const HilbertSym<cx_double>* hil_ptr = hil.getHilPtr();
  cout << "Creating HilSym took " << timer.toc() << endl << endl;

  // creating the Hamiltonian
  GeneralOp<cx_double> ham(hil_ptr);
  for (uword i=0; i<N; i++){
	ham.append(new Ni(hil_ptr, i), 0.123);
	ham.append(new BdagiBj<cx_double>(hil_ptr, i, (i+1)%N), -1.0);
	ham.append(new BdagiBj<cx_double>(hil_ptr, (i+1)%N, i), -1.0);
  }

  ham.print();

  // testing Lanczos
  cout << "Testing Lanczos" << endl;
  cout.precision(8);
  timer.tic();
  ham.doLanczosOnFly(100, 1e-14);
  cout << "Testing Lanczos took " << timer.toc() << endl << endl;

  // testing dense matrix
  /*
  cout << "Creating dense matrix" << endl;
  cout.precision(8);
  timer.tic();
  Mat<cx_double> ham_mat_2;
  ham.createMatrix(&ham_mat_2);
  cout << "Creating dense matrix took " << timer.toc() << endl << endl;
  */

  // testing sparse matrix
  cout << "Creating sparse matrix" << endl;
  cout.precision(8);
  timer.tic();
  SpMat<cx_double> ham_sp_mat_2;
  ham.createMatrix(&ham_sp_mat_2);
  cout << "Creating sparse matrix took " << timer.toc() << endl << endl;

  //
  cx_vec eigval = eigs_gen(ham_sp_mat_2, 2);
  cout << eigval;
  
  return 0;
}
