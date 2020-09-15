#include "../lib/hilbertBase.hpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** ----------------- *** //
  // *** HilbertBase class *** //

  HilbertBase::HilbertBase(const urowvec& localhs_sizes_int){
	this->localhs_sizes_int = localhs_sizes_int;
	this->num_states = prod(localhs_sizes_int);
	this->num_sites = localhs_sizes_int.n_cols;
	HilbertBase::configureMasks();
  }

  urowvec HilbertBase::confInt2confVec(const uword& confint) const{
	urowvec confvec = urowvec(this->num_sites, fill::zeros);
	for (uword i=0; i<this->num_sites; ++i)
	  confvec[i] = HilbertBase::readQnAtSite(i, confint);
	return confvec;
  }

  uword HilbertBase::confVec2confInt(const urowvec& confvec) const{
	if(this->num_sites != confvec.n_elem)
	  throw logic_error("In Hilbert::confVec2confInt."
							 "The size of confVec is incorrect");
	uword confint = 0;
	for (uword i=0; i<this->num_sites; ++i)
	  confint = confint | (confvec[i] << this->localhs_shifts[i]);
	return confint;
  }

  uword HilbertBase::readQnAtSite(uword i, const uword& confint) const{
	return (confint & this->localhs_masks[i]) >> this->localhs_shifts[i];
  }

  uword HilbertBase::changeQnAtSite(uword i, uword qn, const uword& confint) const{
	return (confint & this->localhs_nonmasks[i]) | (qn << this->localhs_shifts[i]);
  }

  void HilbertBase::configureMasks(){
	this->localhs_sizes_bit = log2(this->localhs_sizes_int + 1);
	uword size_max = accu(this->localhs_sizes_bit);
	if(size_max > 64)
	  throw logic_error("In Hilbert::setConfigurationMask."
							 "Number of sites is too big");
	this->localhs_shifts = urowvec(this->num_sites, fill::zeros);
	for(uword i=1; i<this->num_sites; ++i)
	  this->localhs_shifts[i] = this->localhs_shifts[i-1]
		+ this->localhs_sizes_bit[i];
	uword mask_full = ((uword(1) << (size_max - 1)) - 1) << 1 | 1;
	this->localhs_masks = urowvec(this->num_sites, fill::zeros);
	this->localhs_nonmasks = urowvec(this->num_sites, fill::zeros);
	for(uword i=0; i<this->num_sites; ++i){
	  this->localhs_masks[i] = ((uword(1) << this->localhs_sizes_bit[i]) - 1)
		<< this->localhs_shifts[i];
	  this->localhs_nonmasks[i] = mask_full ^ this->localhs_masks[i];
	}
  }

  // *** ----------------- *** //
  
} //* namespace quantum_solver_ed
