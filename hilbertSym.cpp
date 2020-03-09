#ifndef HILBERT_SYM_H
#define HILBERT_SYM_H

#include "hilbertBrut.cpp"
#include "symmetryBones.cpp"
#include "word.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  class HilbertSym: public HilbertBase{
	typedef vector<SymmetryBones*> VectorSymPtr;
  public:
	HilbertSym(const urowvec& localhs_sizes_int);
	const HilbertSym* getHilPtr() const;
	bool checkSymmetric() const;
	void addSym(SymmetryBones* symptr);
	uword getConfInt(uword i) const;
	uword getConfPair(const uword& confint) const;
	PairUwordT<cx_double> getSymPair(const uword& confint) const;
	PairUword getRep(const uword& confint) const;
	int testRep(const uword& confint) const;
	void fill(uword index, const uword& confint, HilMap* confint_map_tmp_ptr);
	void create();
	void createWithFixedQn(uword qn_tot);
	void print() const;
	void printSymmetries() const;
  private:
	uword num_states_nonsym;
	uword tot_shift, half_shift, mask_left, mask_rght;
	VectorSymPtr symptr_vec;
	uvec confint_vec, rep_lookup_vec;
	HilMap confint_map_left, confint_map_rght;
  };

  HilbertSym::HilbertSym(const urowvec& localhs_sizes_int)
	:HilbertBase(localhs_sizes_int){
	this->tot_shift = this->localhs_shifts[this->num_sites - 1];
	this->half_shift = this->localhs_shifts[this->num_sites / 2];
	this->mask_left = (uword(1) << this->half_shift) - 1;
	this->mask_rght = ((uword(1) << this->tot_shift) - 1) << this->half_shift;
  }

  inline
  const HilbertSym* HilbertSym::getHilPtr() const{
	return this;
  }

  inline
  bool HilbertSym::checkSymmetric() const{
	return true;
  }

  inline
  void HilbertSym::addSym(SymmetryBones* symptr){
	this->symptr_vec.push_back(symptr);
  }

  inline
  uword HilbertSym::getConfInt(uword i) const{
	return this->confint_vec[i];
  }

  uword HilbertSym::getConfPair(const uword& confint) const{
	uword conf_left = this->mask_left & confint;
	uword conf_rght = (this->mask_rght & confint) >> this->half_shift;
	return this->confint_map_left.at(conf_left) + this->confint_map_rght.at(conf_rght);
  }

  PairUwordT<cx_double> HilbertSym::getSymPair(const uword& confint) const{
	uword index = HilbertSym::getConfPair(confint);
	Word w(this->rep_lookup_vec[index]);
	PairUwordT<cx_double> sympair;
	sympair.first = w.getRep();
	sympair.second = sqrt(double(w.getDeg()))
	  * this->symptr_vec[w.getSym()]->getChi();
	return sympair;
  }

  PairUword HilbertSym::getRep(const uword& confint) const{
	uword index_min(0), confint_min(confint);
	for (uword i=0; i<this->symptr_vec.size(); ++i){
	  SymmetryBones& sym_ref = *(this->symptr_vec[i]);
	  uword confint_new = this->confVec2confInt(
												sym_ref.apply(this->confInt2confVec(confint)));
	  if (confint_new < confint_min){
		confint_min = confint_new;
		index_min = i;
	  }
	}
	return PairUword(confint_min, index_min);
  }

  int HilbertSym::testRep(const uword& confint) const{
	int deg = 0;
	cx_double sum_chi = 0;
	for (uword i=0; i<this->symptr_vec.size(); ++i){
	  SymmetryBones& sym_ref = *(this->symptr_vec[i]);
	  if (this->confVec2confInt(sym_ref.apply(this->confInt2confVec(confint)))
		  == confint){
		deg ++;
		sum_chi += sym_ref.getChi();
	  }
	}
	if (abs(sum_chi) < 1e-8)
	  return -1;
	else
	  return deg;
  }

  void HilbertSym::fill(uword index, const uword& confint,
						HilMap* confint_map_tmp_ptr){
	PairUword rep_confpair = this->getRep(confint);
	int rep_deg_int = this->testRep(rep_confpair.first);
	if (rep_deg_int == -1)
	  this->rep_lookup_vec[index] = uword(0); 
	else{
	  uword rep_deg = uword(rep_deg_int);
	  uword rep_w = Word::getWrd(rep_confpair.first, rep_deg, 0);
	  auto search = confint_map_tmp_ptr->find(rep_w);
	  // if the rep of the orbit is already in here we check rep_index
	  if (search != confint_map_tmp_ptr->end()){
		uword rep_index = search->second;
		this->rep_lookup_vec[index] = Word::getWrd(rep_index, rep_deg,
												   rep_confpair.second);
	  }
	  // if not we add rep to the orbit and
	  // we get rep_index = num_states; num_states ++
	  // ATTENTION: rep_lookup_vec[index] gives index of the rep, not the rep
	  else{
		confint_map_tmp_ptr->insert(PairUword(rep_w, this->num_states));
		this->confint_vec[this->num_states] = rep_w;
		this->rep_lookup_vec[index] = Word::getWrd(this->num_states, rep_deg,
												   rep_confpair.second);
		this->num_states ++;
	  }
	}
  }

  void HilbertSym::create(){
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
	hil_left.create();
	hil_rght.create();
  
	//* some variables
	uword conf_left, conf_rght, conf_new, index_rght;
	HilMap confint_map_tmp;

	//* number of states in the symmetry sector is yet unknown
	this->num_states_nonsym = this->num_states;
	this->rep_lookup_vec.set_size(this->num_states_nonsym);
	this->confint_vec.set_size(this->num_states_nonsym);
                             
	//* making map for states in the left Hilbert space
	for(uword i=0; i<hil_left.getNumStates(); ++i){
	  conf_left = hil_left.getConfInt(i);
	  this->confint_map_left.insert(PairUword(conf_left, i));
	}

	//* combining two Hilbert spaces
	this->num_states = 0;
	for(uword ir=0; ir<hil_rght.getNumStates(); ++ir){
	  conf_rght = hil_rght.getConfInt(ir);
	  index_rght = ir * hil_left.getNumStates();
	  this->confint_map_rght.insert(PairUword(conf_rght, index_rght));
	  conf_rght = conf_rght << this->half_shift;
	  for(uword il=0; il<hil_left.getNumStates(); ++il){
		conf_new = hil_left.getConfInt(il) | conf_rght;
		HilbertSym::fill(index_rght + il, conf_new, &confint_map_tmp);
	  }
	}
	this->confint_vec.resize(this->num_states); // updated inside HilbertSym::fill
	cout << "Hilbert size. Tot size: " << this->num_states_nonsym << ". "
		 << "Sym size: " << this->num_states << endl;
  }

  void HilbertSym::createWithFixedQn(uword qn_tot = 0){
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
	uword qr, index_new, conf_new, conf_left, conf_rght, index(0);
	HilMap confint_map_tmp;
	PairUword confpair;
    
	//* making the map for states in the left Hilbert space
	this->num_states_nonsym = 0;
	for(uvec::iterator qlp=qns_left_unique.begin(); qlp<qns_left_unique.end(); ++qlp){
	  qr = qn_tot - *qlp;
	  inds_rght = find(qns_rght == qr);
	  if (inds_rght.n_elem != 0){
		inds_left = find(qns_left == *qlp);
		for(uword irl=0; irl<inds_left.n_elem; ++irl){
		  conf_left = hil_left.getConfInt(inds_left[irl]);
		  this->confint_map_left.insert(PairUword(conf_left, irl));
		}
		this->num_states_nonsym += inds_rght.n_elem * inds_left.n_elem;
	  }
	}
	this->rep_lookup_vec.set_size(this->num_states_nonsym);
	this->confint_vec.set_size(this->num_states_nonsym);
    
	//* combining two HS's
	this->num_states = 0;
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
			index_new = index + irp * inds_left.n_elem + irl;
			conf_new = hil_left.getConfInt(inds_left[irl]) | conf_rght;
			HilbertSym::fill(index_new, conf_new, &confint_map_tmp);
		  }
		}
		index += inds_rght.n_elem * inds_left.n_elem;
	  }
	}
	this->confint_vec.resize(this->num_states); // updated inside HilbertSym::fill
	cout << "Hilbert size. Tot size: " << this->num_states_nonsym << ". "
		 << "Sym size: " << this->num_states << endl;
  }

  void HilbertSym::print() const{
	cout << "Printing the Hilbert space" << endl;
	for(uword i=0; i<this->num_states; ++i){
	  cout << i << ") ";
	  for(uword j=0; j<this->num_sites; ++j){
		Word w(this->confint_vec[i]);
		cout << this->confInt2confVec(w.getRep())[j] << " ";
	  }
	  cout << endl;
	}
  }

  void HilbertSym::printSymmetries() const{
	cout << "Printing Hilbert space symmetries" << endl;
	for(uword i=0; i<this->symptr_vec.size(); ++i){
	  cout << *(this->symptr_vec[i]) << endl;
	}
  }

} //* namespace quantum_solver_ed


#endif //* HILBERT_SYM_H
