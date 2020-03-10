#ifndef HILBERT_H
#define HILBERT_H

#include "hilbertBrut.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  class Hilbert: public HilbertBase{
	// same as HilbertBrut but has two maps to search faster
  public:
	Hilbert(const urowvec& localhs_sizes_int);
	const Hilbert* getHilPtr() const;
	uword getConfInt(uword i) const;
	uword getConfPair(const uword& confint) const;
	void create();
	void createWithFixedQn(uword qn_tot=0);
	void print() const;
  protected:
	uword half_shift, mask_left, mask_rght;
	uvec confint_vec;
	HilMap confint_map_left, confint_map_rght;
  };

  Hilbert::Hilbert(const urowvec& localhs_sizes_int)
	:HilbertBase(localhs_sizes_int){
	uword tot_shift = this->localhs_shifts[this->num_sites - 1];
	this->half_shift = this->localhs_shifts[this->num_sites / 2];
	this->mask_left = (uword(1) << this->half_shift) - 1;
	this->mask_rght = ((uword(1) << tot_shift) - 1) << this->half_shift;
  }

  inline
  const Hilbert* Hilbert::getHilPtr() const{
	return this;
  }

  inline
  uword Hilbert::getConfInt(uword i) const{
	return this->confint_vec[i];
  }

  uword Hilbert::getConfPair(const uword& confint) const{
	uword conf_left = this->mask_left & confint;
	uword conf_rght = (this->mask_rght & confint) >> this->half_shift;
	return this->confint_map_left.at(conf_left) + this->confint_map_rght.at(conf_rght);
  }

  void Hilbert::create(){
	//* creating left and right Hilbert spaces
	uword ns_left = this->num_sites / 2;
	uword ns_rght = this->num_sites - ns_left;
	urowvec localhs_sizes_int_left = urowvec(ns_left, fill::zeros);
	urowvec localhs_sizes_int_rght = urowvec(ns_rght, fill::zeros);
	for(uword i=0; i<this->num_sites; ++i){
	  if(i < ns_left)
		localhs_sizes_int_left[i] = this->localhs_sizes_int[i];
	  else
		localhs_sizes_int_rght[i-ns_left] = this->localhs_sizes_int[i];
	}
	HilbertBrut hil_left(localhs_sizes_int_left);
	HilbertBrut hil_rght(localhs_sizes_int_rght);
	cout << this->num_states << endl;
	hil_left.create();
	hil_rght.create();
	cout << this->num_states << endl;

	//* some variables
	uword conf_left, conf_rght, conf_new, index_rght;
	this->confint_vec.set_size(this->num_states);

	//* making map for states in the left Hilbert space
	for(uword i=0; i<hil_left.getNumStates(); ++i){
	  conf_left = hil_left.getConfInt(i);
	  this->confint_map_left.insert(PairUword(conf_left, i));
	}

	//* combining two Hilbert spaces
	for(uword ir=0; ir<hil_rght.getNumStates(); ++ir){
	  conf_rght = hil_rght.getConfInt(ir);
	  index_rght = ir * hil_left.getNumStates();
	  this->confint_map_rght.insert(PairUword(conf_rght, index_rght));
	  conf_rght = conf_rght << this->half_shift;
	  for(uword il=0; il<hil_left.getNumStates(); ++il){
		conf_new = hil_left.getConfInt(il) | conf_rght;
		this->confint_vec[index_rght + il] = conf_new;
	  }
	}
  }

  void Hilbert::createWithFixedQn(uword qn_tot){
	//* creating right and left HS's
	uword ns_left = this->num_sites / 2;
	uword ns_rght = this->num_sites - ns_left;
	urowvec localhs_sizes_int_left = urowvec(ns_left, fill::zeros);
	urowvec localhs_sizes_int_rght = urowvec(ns_rght, fill::zeros);
	for(uword i=0; i<this->num_sites; ++i){
	  if(i < ns_left)
		localhs_sizes_int_left[i] = this->localhs_sizes_int[i];
	  else
		localhs_sizes_int_rght[i-ns_left] = this->localhs_sizes_int[i];
	}
	HilbertBrut hil_left(localhs_sizes_int_left);
	HilbertBrut hil_rght(localhs_sizes_int_rght);
	hil_left.create();
	hil_rght.create();
    
	//* preparing for sorting right and left HS's by qns
	uvec qns_left = uvec(hil_left.getNumStates(), fill::zeros);
	for(uword i=0; i<hil_left.getNumStates(); ++i)
	  qns_left[i] = accu(hil_left.confInt2confVec(hil_left.getConfInt(i)));
	uvec qns_rght = uvec(hil_rght.getNumStates(), fill::zeros);
	for(uword i=0; i<hil_rght.getNumStates(); ++i)
	  qns_rght[i] = accu(hil_rght.confInt2confVec(hil_rght.getConfInt(i)));
	uvec qns_left_unique = uvec(unique(qns_left));
	uvec qns_rght_unique = uvec(unique(qns_rght));

	//* some variables
	uvec inds_left, inds_rght;
	uword qr, conf_new, conf_left, conf_rght, num_states(0), index(0);
	uvec confint_vec_tmp = uvec(0, fill::zeros);
	PairUword confpair;
  
	//* making the map for states in the left Hilbert space
	for(uvec::iterator qlp=qns_left_unique.begin(); qlp<qns_left_unique.end(); ++qlp){
	  qr = qn_tot - *qlp;
	  inds_rght = find(qns_rght == qr);
	  if (inds_rght.n_elem != 0){
		inds_left = find(qns_left == *qlp);
		for(uword irl=0; irl<inds_left.n_elem; ++irl){
		  conf_left = hil_left.getConfInt(inds_left[irl]);
		  this->confint_map_left.insert(PairUword(conf_left, irl));
		}
		num_states += inds_rght.n_elem * inds_left.n_elem;
	  }
	}
	this->confint_vec.set_size(num_states);
  
	//* combining two HS's
	for(uvec::iterator qlp=qns_left_unique.begin(); qlp<qns_left_unique.end(); ++qlp){
	  qr = qn_tot - *qlp;
	  inds_rght = find(qns_rght == qr);
	  if(inds_rght.n_elem != 0){
		inds_left = find(qns_left == *qlp);
		for(uword irp=0; irp<inds_rght.n_elem; ++irp){
		  conf_rght = hil_rght.getConfInt(inds_rght[irp]);
		  this->confint_map_rght.insert(PairUword(conf_rght,
												  index + irp * inds_left.n_elem));
		  conf_rght = conf_rght << this->half_shift;
		  for(uword irl=0; irl<inds_left.n_elem; ++irl){
			conf_new = hil_left.getConfInt(inds_left[irl]) | conf_rght;
			this->confint_vec[index + irp * inds_left.n_elem + irl] = conf_new;
		  }
		}
		index += inds_rght.n_elem * inds_left.n_elem;
	  }
	}
	this->num_states = num_states;
  }

  void Hilbert::print() const{
	cout << "Printing the Hilbert space" << endl;
	for(uword i=0; i<this->num_states; ++i){
	  cout << i << ") ";
	  for(uword j=0; j<this->num_sites; ++j)
		cout << this->confInt2confVec(this->confint_vec[i])[j] << " ";
	  cout << endl;
	}
  }

} //* namespace quantum_solver_ed


#endif //* HILBERT_H
