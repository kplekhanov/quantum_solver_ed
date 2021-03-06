#include "../lib/hilbertSym.hpp"
#include "../lib/waveFon.hpp"
#include "../lib/generalOp.hpp"
#include "../lib/opHardCoreBose.hpp"
#include "../lib/opSpinOneHalf.hpp"


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
  hil.readSymFromFile("input/input_sym_10_complex.csv");
  hil.createWithFixedQn(Qn);
  const HilbertBones* hil_ptr = hil.getHilPtr();
  cout << "Creating HilSym took " << timer.toc() << endl;
  hil.printSymmetries();
  cout << endl;

  GeneralOp<cx_double> ham(hil_ptr);
  ham.readFromFile<SpiSmj_2<cx_double>>("input/input_ops_2is_complex.csv");
  ham.readFromFile<Szi_2>("input/input_ops_1i.csv");
  ham.print();
  cout << endl;

  GeneralOp<cx_double> ham1(hil_ptr);
  ham1.readFromFile<BdagiBj_2<cx_double>>("input/input_ops_2is_complex.csv");
  ham1.readFromFile<Ni_2>("input/input_ops_1i.csv");
  ham1.print();
  cout << endl;

  return 0;
}
