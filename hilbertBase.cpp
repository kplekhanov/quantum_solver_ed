#ifndef HILBERT_BASE_H
#define HILBERT_BASE_H

#include "hilbertBones.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  class HilbertBase: public HilbertBones{
  public:
	HilbertBase(const urowvec& localhs_sizes_int);
	const HilbertBase* getHilPtr() const;
	uword getNumStates() const;
	uword getNumSites() const;
	uword getLocalhsSize(uword i) const;
	virtual bool checkSymmetric() const; // polymorphism for HilbertSym
	urowvec confInt2confVec(const uword& confint) const;
	uword confVec2confInt(const urowvec& confvec) const;
	uword readQnAtSite(uword i, const uword& confint) const;
	uword changeQnAtSite(uword i, uword qn, const uword& confint) const;
	void configureMasks();
	virtual void create() = 0;
	virtual void createWithFixedQn(uword qn_tot) = 0;
	virtual void print() const = 0;
  protected:
	uword num_states, num_sites;
	urowvec localhs_sizes_int;
	urowvec localhs_sizes_bit;
	urowvec localhs_shifts;
	urowvec localhs_masks, localhs_nonmasks;
  };

  HilbertBase::HilbertBase(const urowvec& localhs_sizes_int){
	this->localhs_sizes_int = localhs_sizes_int;
	this->num_states = prod(localhs_sizes_int);
	this->num_sites = localhs_sizes_int.n_cols;
	HilbertBase::configureMasks();
  }

  inline
  const HilbertBase* HilbertBase::getHilPtr() const{
	return this;
  }

  inline
  uword HilbertBase::getNumStates() const{
	return this->num_states;
  }

  inline
  uword HilbertBase::getNumSites() const{
	return this->num_sites;
  }

  inline
  uword HilbertBase::getLocalhsSize(uword i) const{
	return this->localhs_sizes_int[i];
  }

  inline
  bool HilbertBase::checkSymmetric() const{
	return false;
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

  inline
  uword HilbertBase::readQnAtSite(uword i, const uword& confint) const{
	return (confint & this->localhs_masks[i]) >> this->localhs_shifts[i];
  }

  inline
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

} //* namespace quantum_solver_ed

  
#endif //* HILBERT_BASE_H
