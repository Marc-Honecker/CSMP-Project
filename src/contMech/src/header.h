#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <cmath>
#include <complex>
#include <fftw3.h>
#include <algorithm>
#include <bitset>

#include <vector>
#include <stdlib.h>
#include <cstring>

#include<ctime>

using namespace std;

typedef ptrdiff_t Lint;
typedef std::complex<double> Complex;

#ifndef _PI_
  #define _PI_
  const double PI = 4.*atan(1.), TWOPI = 2*PI;
#endif

template <class T> int sign(T a){
  return(a<0 ? -1 : 1);
}

