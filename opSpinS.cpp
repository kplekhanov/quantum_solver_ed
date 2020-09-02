#ifndef OP_SPIN_H
#define OP_SPIN_H

#include "elementaryOp.cpp"


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;

  // *** -------------- *** //
  // *** diagonal terms *** //

  // *** Zeeman along the z-axis
  class Szi: public DiagOp{
  public:
	Szi(){};
	Szi(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
  };

  Szi::Szi(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void Szi::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 1){
	  throw logic_error( "In Szi::fill. Vector of indices should have size 1." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->name = "Sz_" + to_string(this->i);
  }

  double Szi::apply(const uword& conf_in) const{
	double szi = this->hil_ptr->readQnAtSite(this->i, conf_in);
	szi += (1 - double(this->hil_ptr->getLocalhsSize(this->i))) / 2;
	return szi;
  }

  // *** Exchange along the z-axis
  class SziSzj: public DiagOp{
  public:
	SziSzj(){};
	SziSzj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	double apply(const uword& conf_in) const;
  private:
	uword i;
	uword j;
  };

  SziSzj::SziSzj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  void SziSzj::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSzj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sz_" + to_string(this->j);
  }

  double SziSzj::apply(const uword& conf_in) const{
	double szi = (this->hil_ptr->readQnAtSite(this->i, conf_in) +
				  (1 - double(this->hil_ptr->getLocalhsSize(this->i))) / 2);
	double szj = (this->hil_ptr->readQnAtSite(this->j, conf_in) +
				  (1 - double(this->hil_ptr->getLocalhsSize(this->j))) / 2);
	return szi * szj;
  }   

  // *** -------------- *** //


  // *** ------------------ *** //
  // *** off-diagonal terms *** //

  // *** swap gates
  template<typename T=double>
  class SpiSmj: public OffDiagOp<T>{
  public:
	SpiSmj(){};
	SpiSmj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SpiSmj<T>::SpiSmj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SpiSmj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SpiSmj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sp_" + to_string(this->i) + ".Sm_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SpiSmj<T>::apply(const uword& conf_in) const{
	uword Ni = this->hil_ptr->getLocalhsSize(this->i);
	uword qni = this->hil_ptr->readQnAtSite(this->i, conf_in);
	uword qnj = this->hil_ptr->readQnAtSite(this->j, conf_in);
	if ((qni == Ni-1) | (qnj == 0))
	  return PairUwordT<T>(0, 0);
	uword Nj = this->hil_ptr->getLocalhsSize(this->j);
	uword conf_new = this->hil_ptr->changeQnAtSite(this->i, qni+1, conf_in);
	conf_new = this->hil_ptr->changeQnAtSite(this->j, qnj-1, conf_new);
	double Si = (double(Ni) - 1) / 2;
	double Sj = (double(Nj) - 1) / 2;
	double szi = qni - Si;
	double szj = qnj - Sj;
	return PairUwordT<T>(conf_new, sqrt(Si*(Si+1)-szi*(szi+1))*sqrt(Sj*(Sj+1)-szj*(szj-1)));
  }


  // *** in-plane DM interactions -- Sz Sp
  template<typename T=double>
  class SziSpj: public OffDiagOp<T>{
  public:
	SziSpj(){};
	SziSpj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SziSpj<T>::SziSpj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSpj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSpj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sp_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSpj<T>::apply(const uword& conf_in) const{
	uword Nj = this->hil_ptr->getLocalhsSize(this->j);
	uword qnj = this->hil_ptr->readQnAtSite(this->j, conf_in);
	if (qnj == Nj-1)
	  return PairUwordT<T>(0, 0);
	uword Ni = this->hil_ptr->getLocalhsSize(this->i);
	uword qni = this->hil_ptr->readQnAtSite(this->i, conf_in);
	double szi = qni + (1 - double(Ni)) / 2;
	double Sj = (double(Nj) - 1) / 2;
	double szj = qnj - Sj;
	uword conf_new = this->hil_ptr->changeQnAtSite(this->j, qnj+1, conf_in);
	return PairUwordT<T>(conf_new, szi*sqrt(Sj*(Sj+1)-szj*(szj+1)));
  }

  
  // *** in-plane DM interactions -- Sz Sm
  template<typename T=double>
  class SziSmj: public OffDiagOp<T>{
  public:
	SziSmj(){};
	SziSmj(const HilbertBones* hil_ptr, vector<uword> indices_int);
	void fill(const HilbertBones* hil_ptr, vector<uword> indices_int);
	PairUwordT<T> apply(const uword& conf_in) const;
  private:
	uword i, j;
  };

  template<typename T>
  SziSmj<T>::SziSmj(const HilbertBones* hil_ptr, vector<uword> indices_int){
	this->fill(hil_ptr, indices_int);
  }

  template<typename T>
  void SziSmj<T>::fill(const HilbertBones* hil_ptr, vector<uword> indices_int){
	if (indices_int.size() != 2){
	  throw logic_error( "In SziSmj::fill. Vector of indices should have size 2." );
	}
	this->hil_ptr = hil_ptr;
	this->i = indices_int.at(0);
	this->j = indices_int.at(1);
	this->name = "Sz_" + to_string(this->i) + ".Sm_" + to_string(this->j);
  }

  template<typename T>
  PairUwordT<T> SziSmj<T>::apply(const uword& conf_in) const{
	uword Nj = this->hil_ptr->getLocalhsSize(this->j);
	uword qnj = this->hil_ptr->readQnAtSite(this->j, conf_in);
	if (qnj == 0)
	  return PairUwordT<T>(0, 0);
	uword Ni = this->hil_ptr->getLocalhsSize(this->i);
	uword qni = this->hil_ptr->readQnAtSite(this->i, conf_in);
	double szi = qni + (1 - double(Ni)) / 2;
	double Sj = (double(Nj) - 1) / 2;
	double szj = qnj - Sj;
	uword conf_new = this->hil_ptr->changeQnAtSite(this->j, qnj-1, conf_in);
	return PairUwordT<T>(conf_new, szi*sqrt(Sj*(Sj+1)-szj*(szj-1)));
  }
  
  // *** ------------------ *** //

} //* namespace quantum_solver_ed


#endif //* OP_SPIN_H
