#include "../lib/elementaryOp.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ---------------------------- *** //
  // *** class of elemetary operators *** //

  void ElementaryOp::print(ostream& os) const{
	os << this->name;
  }

  void ElementaryOp::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->hil_ptr = hil_ptr;
	cout << "Warning: filling is not implemented for this operator" << endl;
	cout << "Reading indices: ";
	for (auto i=indices_int.begin(); i<indices_int.end(); i++){
	  cout << *i << " ";
	}
	cout << endl;
  }

  // *** functions
  ostream& operator<<(ostream& os, const ElementaryOp& eop){
	eop.print(os);
	return os;
  }
  
  // *** ------------------------ *** //

  // *** ------------------------ *** //
  // *** instatation of template  *** //
  
  template class OffDiagOp<double>;
  template class OffDiagOp<cx_double>;

  // *** ------------------------ *** //

} //* namespace quantum_solver_ed
