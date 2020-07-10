#include "../hilbertSym.cpp"
#include "../waveFon.cpp"
#include "../generalOp.cpp"
//#include "../opBoseHubbard.cpp"
#include "../opHardCoreBose.cpp"
#include "../opSpinOneHalf.cpp"


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
  HilbertSym<cx_double> hil(localSizes);
  hil.readSymFromFile("input_sym_10_complex.csv");
  hil.createWithFixedQn(Qn);
  const HilbertSym<cx_double>* hil_ptr = hil.getHilPtr();
  cout << "Creating HilSym took " << timer.toc() << endl;
  hil.printSymmetries();

  GeneralOp<cx_double> ham(hil_ptr);
  ham.readFromFile<SpiSmj<cx_double>>("input_ops_2is_complex.csv");
  ham.readFromFile<Szi>("input_ops_1i.csv");
  ham.print();

  GeneralOp<cx_double> ham1(hil_ptr);
  ham1.readFromFile<BdagiBj<cx_double>>("input_ops_2is_complex.csv");
  ham1.readFromFile<Ni>("input_ops_1i.csv");
  ham1.print();

  return 0;
}
