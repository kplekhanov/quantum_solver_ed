#ifndef HILBERT_BRUT_H
#define HILBERT_BRUT_H

#include "hilbertBase.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  class HilbertBrut: public HilbertBase{
  public:
	HilbertBrut(const urowvec& localhs_sizes_int);
	const HilbertBrut* getHilPtr() const;
	uword getConfInt(uword i) const;
	uword getConfPair(const uword& confint) const;
	void create();
	void createWithFixedQn(uword qn_tot);
	void print() const;
  protected:
	uvec confint_vec;
	HilMap confint_map;
  };

  HilbertBrut::HilbertBrut(const urowvec& localhs_sizes_int)
	:HilbertBase(localhs_sizes_int){}

  inline
  const HilbertBrut* HilbertBrut::getHilPtr() const{
	return this;
  }

  inline
  uword HilbertBrut::getConfInt(uword i) const{
	return this->confint_vec[i];
  }

  inline
  uword HilbertBrut::getConfPair(const uword& confint) const{
	return this->confint_map.at(confint);
  }

  void HilbertBrut::create(){
	uword localhs_size, num_states, localhs_shift;
	uword conf_local, conf_new;
	uvec confint_vec_tmp = uvec(this->num_states, fill::zeros);
  
	num_states = this->localhs_sizes_int[0];
	for(uword i=0; i<num_states; ++i)
	  confint_vec_tmp[i] = i;
	for(uword i=1; i<num_sites; ++i){
	  localhs_size = this->localhs_sizes_int[i];
	  localhs_shift = this->localhs_shifts[i];
	  this->confint_vec.set_size(localhs_size * num_states);
	  for(uword j=0; j<localhs_size; ++j){
		conf_local = j << localhs_shift;
		for(uword k=0; k<num_states; ++k){
		  conf_new = confint_vec_tmp[k] | conf_local;
		  this->confint_vec[j * num_states + k] = conf_new;
		}
	  }
	  num_states *= localhs_size;
	  if(num_states != this->num_states){
		confint_vec_tmp.set_size(num_states);
		confint_vec_tmp = uvec(this->confint_vec);   
	  }
	}
  
	PairUword confpair;
	for(uword i=0; i<num_states; ++i){
	  confpair.first = this->confint_vec[i];
	  confpair.second = i;
	  this->confint_map.insert(confpair);
	}
  }

  void HilbertBrut::createWithFixedQn(uword qn_tot=0){
	//* creating right and left HS's
	uword ns_left = this->num_sites / 2;
	uword ns_rght = this->num_sites - ns_left;
	uword half_shift = this->localhs_shifts[ns_left];
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
	uvec qns_left = uvec(hil_left.num_states, fill::zeros);
	for(uword i=0; i<hil_left.num_states; ++i)
	  qns_left[i] = accu(hil_left.confInt2confVec(hil_left.confint_vec[i]));
	uvec qns_rght = uvec(hil_rght.num_states, fill::zeros);
	for(uword i=0; i<hil_rght.num_states; ++i)
	  qns_rght[i] = accu(hil_rght.confInt2confVec(hil_rght.confint_vec[i]));
	uvec qns_left_unique = uvec(unique(qns_left));
	uvec qns_rght_unique = uvec(unique(qns_rght));
    
	//* combining two HS's
	uvec inds_left, inds_rght;
	uword qr, conf_new, conf_rght, num_states, index(0);
	uvec confint_vec_tmp = uvec(0, fill::zeros);
	PairUword confpair;
	for(uvec::iterator qlp=qns_left_unique.begin(); qlp<qns_left_unique.end(); ++qlp){
	  qr = qn_tot - *qlp;
	  inds_rght = find(qns_rght == qr);
	  if(inds_rght.n_elem != 0){
		inds_left = find(qns_left == *qlp);
		num_states = index + inds_rght.n_elem * inds_left.n_elem;
		this->confint_vec.set_size(num_states);
		for(uword i=0; i<index; ++i)
		  this->confint_vec[i] = confint_vec_tmp[i];
		for(uword irp=0; irp<inds_rght.n_elem; ++irp){
		  conf_rght = hil_rght.confint_vec[inds_rght[irp]] << half_shift;
		  for(uword irl=0; irl<inds_left.n_elem; ++irl){
			conf_new = hil_left.confint_vec[inds_left[irl]] | conf_rght;
			this->confint_vec[index + irp * inds_left.n_elem + irl] = conf_new;
			confpair.first = conf_new;
			confpair.second = index + irp * inds_left.n_elem + irl;
			this->confint_map.insert(confpair);
		  }
		}
		confint_vec_tmp.set_size(num_states);
		confint_vec_tmp = uvec(this->confint_vec);
		index = num_states;
	  }
	}
	this->num_states = num_states;
  }

  void HilbertBrut::print() const{
	cout << "Printing the Hilbert space" << endl;
	for(uword i=0; i<this->num_states; ++i){
	  cout << i << ") ";
	  for(uword j=0; j<this->num_sites; ++j)
		cout << this->confInt2confVec(this->confint_vec[i])[j] << " ";
	  cout << endl;
	}
  }

} //* namespace quantum_solver_ed


#endif //* HILBERT_BRUT_H
