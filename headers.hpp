#ifndef HEADERS_H
#define HEADERS_H

#define OMP_NUM_THREADS 32

#include <armadillo>
#include <stdexcept>
#include <iostream>
#include <string>
#include <sstream> // stringstream
#include <unordered_map>
#include <vector>
#include <math.h>
#include <omp.h> // a priori not needed since using only pragma


namespace quantum_solver_ed{
  using namespace arma;
  using namespace std;
  
  typedef unordered_map<uword,uword> HilMap; // map intconf -> index
  typedef pair<uword,uword> PairUword;
  template <typename T> using PairUwordT = pair<uword,T>;
  template <typename T> using MyVec = Col<T>;
} //* namespace quantum_solver_ed


#endif //* HEADERS_H
