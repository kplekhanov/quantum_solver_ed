#ifndef HEADERS_H
#define HEADERS_H

#define SP << " " <<
#define OMP_NUM_THREADS std::max(atoi(std::getenv("MKL_NUM_THREADS")), 1)

#include <armadillo>
#include <stdexcept>
#include <iostream>
#include <string> 
#include <unordered_map>
#include <vector>
#include <math.h>
#include <omp.h>


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;
  
  typedef unordered_map<uword,uword> HilMap; // map intconf -> index
  typedef pair<uword,uword> PairUword;
  template <typename T> using PairUwordT = pair<uword,T>;
  template <typename T> using MyVec = Col<T>;
  
} //* namespace quantum_solver_ed


#endif //* HEADERS_H
