// MM2All: need to correct reading of data with nx = 256, ny = 1.

#include "header.h"

extern double lengthX, lengthY, areaXY;
extern double fireRedrct;
extern int iTime;

// odd integer with a decent mix of zeroes and ones
uint64_t seed = 16273530647158678139ULL;
// odd integer closest to 2^64/((1+sqrt(5))/2);
// arbitrary; good mix of zeroes and ones;
const uint64_t incr = 0x9e3779b97f4a7c15L;
const double DOUBLE_ULP = 1.0 / (1L << 53);


void termination(const string& notice){
  cerr << "\n### " << notice;
  cerr << "\n### termination !!!" << endl << endl;
  exit(0);
}

// MM2: not sure how to use templates here, i.e.,
//	extern templates appear difficult to define

void termination(const string& notice, int flag){
  cerr << "\n### " << notice << flag;
  cerr << "\n### termination" << endl << endl;
  exit(0);
}

void termination(const string& notice, double flag){
  cerr << "\n### " << notice << flag;
  cerr << "\n### termination" << endl << endl;
  exit(0);
}

double rrand() {return(rand()*1./RAND_MAX);}

double mix64(void){
  uint64_t z = ( seed += incr );
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;   // MurmurHash3 mix constants
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  //return (z ^ (z >> 31)); // in case one wants uint64_t output
  return ((z ^ (z >> 31)) >> 11) * DOUBLE_ULP; // a 64-bit double really only has 53 bits (the rest is exponent)
}

void changeSeedMix64(int mySeed){
  seed += 2*mySeed + 1;
}

// input/output functions

// CS = cross section
void dumpRealCS(vector<double*> fields, vector<string> columnNames, Lint nx, Lint ny, int xStep, int yStep, const string& fileName) {
  double dx = lengthX/nx, dy = lengthY/ny;
  ofstream dumPH(fileName);

  // create file header
  string header = (ny>nx)? "# y" : "# x";
  if (nx > ny) {
    for (string column:columnNames) header += "\t" + column + "(x)";
  }
  else if (ny > nx) {
    for (string column:columnNames) header += "\t" + column + "(y)";
  }
  else {
    for (string column:columnNames) header += "\t" + column + "(x)\t" + column + "(y)";
    header += "\td";
    for (string column:columnNames) header += "\t" + column + "(d)";
  }
  header += "\n";
  dumPH << header;

  // only plot along x
  if (nx > ny) {
    int iy = ny/2;
    for (int ix=0; ix<=nx; ix+=xStep) {
      dumPH << (ix-nx/2)*dx;
      for (double* field:fields) dumPH << "\t" << field[(ix%nx)*ny+iy];
      dumPH << "\n";
    }
  }

  // only plot along y
  else if (ny > nx) {
    int ix = nx/2;
    for (int iy=0; iy<=ny; iy+=yStep) {
      dumPH << (iy-ny/2)*dy;
      for (double* field:fields) dumPH << "\t" << field[ix*ny+(iy%ny)];
      dumPH << "\n";
    }
  }

  // plot along x, y, and diagonal
  else {
    for (int ix=0; ix<=nx; ix+=xStep) {
      double x = (ix-nx/2)*dx, y = (ix-nx/2)*dy, r = sqrt(x*x+y*y);
      dumPH << x;
      for (double* field:fields) dumPH << "\t" << field[(ix%nx)*ny+ny/2] << "\t" << field[nx*ny/2+(ix%nx)];
      dumPH << "\t";
      if(x<0) dumPH << "-";
      dumPH << r;
      for (double* field:fields) dumPH << "\t" << field[(ix%nx)*ny+(ix%nx)];
      dumPH << "\n";
    }
  }

  dumPH.close();
}


void dumpReal(vector<double*> fields, Lint nx, Lint ny, int xStep, int yStep, const string& fileName) {
  //change: if input/output ever becomes advanced enough to not break completely
  //        because of it, we should add column headers here as well.

  double dx = lengthX/nx, dy = lengthY/ny;
  ofstream dumP(fileName);

  // create file header
  int nxNew = nx/xStep, nyNew = ny/yStep;
  dumP << "#" << nxNew << "\t" << nyNew << "\n\n";

  // only plot along y
  if(nx==1) {
    for (Lint iy=0; iy<=ny; iy+=yStep) {
      dumP << iy*dy;
      for (double* field:fields) dumP << "\t" << field[iy];
      dumP << "\n";
    }
  }

  // only plot along x
  else if (ny==1) {
    for (Lint ix=0; ix<=nx; ix+=xStep) {
      dumP << ix*dx;
      for (double* field:fields) dumP << "\t" << field[ix];
      dumP << "\n";
    }
  }

  // plot along x, y, and diagonal
  else {
    for (Lint ix=0; ix<=nx; ix+=xStep) {
      for (Lint iy=0; iy<=ny; iy+=yStep) {
        Lint k = (ix%nx)*ny + (iy%ny);
        dumP << ix*dx << "\t" << iy*dy;
        for (double* field:fields) dumP << "\t" << field[k];
        dumP << "\n";
      }
      dumP << "\n";
    }
  }

  dumP.close();
}


void dumpFour(Complex* array, Lint nx, Lint ny, const string& fileName){
  Lint nyHP1 = (ny/2) + 1;
  ofstream dumP(fileName);
  dumP << "#" << nx << "\t" << ny << "\n\n";
  for (Lint ix=0; ix<nx; ++ix) {
    for (Lint iy=0; iy<nyHP1; ++iy) {
      Lint k = ix*nyHP1+iy;
      dumP << ix << "\t" << iy << "\t" << array[k] << "\n";
    } dumP << endl;
  } dumP.close();
}

double biLinInterpol(double q11, double q12, double q21, double q22, double x1,
       double x2, double y1, double y2, double x, double y){
  double x2x1, y2y1, x2x, y2y, yy1, xx1;
  x2x1 = x2 - x1;
  y2y1 = y2 - y1;
  x2x  = x2 - x;
  y2y  = y2 - y;
  yy1  = y  - y1;
  xx1  = x  - x1;
  return (q11*x2x*y2y + q21*xx1*y2y + q12*x2x*yy1 + q22*xx1*yy1) / (x2x1*y2y1);
}

int readReal(double* array, Lint nx, Lint ny, const string& fileName, int column){

  ifstream configIn(fileName);
  if (!configIn.is_open()) return(0);

  int fSizeUpScale = 0;

  Lint nxOld, nyOld;
  double dummy, dispZ;
  string str;
  char dummyChar; //used to hold the "#" character in array.old
  configIn >> dummyChar >> nxOld >> nyOld;

  // sanity check
  if ( ((nx%nxOld)!=0)&&((nxOld%nx)!=0)) {
    cerr << "### old config file incompatible with nx. Caution!\n";
    return(1);
  }
  if ( ((ny%nyOld)!=0)&&((nyOld%ny)!=0)) {
    cerr << "### old config file incompatible with ny. Caution!\n";
    return(1);
  }

  vector <double> arrayOld(nxOld*nyOld);
  for (Lint ixOld = 0; ixOld < nxOld; ixOld++){
    for (Lint iyOld = 0; iyOld < nyOld; iyOld++){
      for (int i = 1; i < column; ++i) configIn >> dummy; // skip other columns
      configIn >> dispZ; getline(configIn, str);
      Lint ii = ixOld*nyOld+iyOld;
      arrayOld[ii] = dispZ;
    }
    getline(configIn, str); // skip extra line
  }
  configIn.close();


  if ((nx<=nxOld)&&(ny<=nyOld)){ // 4 cases
    int incrX = nxOld/nx, incrY = nyOld/ny;
    for (Lint ix = 0; ix < nx; ix++){
      for (Lint iy = 0; iy < ny; iy++){
        int k = ix*ny + iy;
        Lint ixOld = ix*incrX;
        Lint iyOld = iy*incrY;
        int kOld = ixOld*nyOld + iyOld;
        array[k] = arrayOld[kOld];
      }
    }
  } else if ((nx>nxOld)&&(ny>nyOld)){ // 3 cases remaining
    fSizeUpScale = 1;
    int scalX = nx/nxOld, scalY = ny/nyOld;
    for (Lint ix = 0; ix < nx; ix++){
    for (Lint iy = 0; iy < ny; iy++){
      Lint k = ix*ny+iy;
      Lint ixOldLeft = ((k)/(ny*scalY));
      Lint iyOldLow = (((k/scalX)%(ny))% (ny/scalY));
      Lint ixOldRight =(((k)/(ny*scalY))+1+nxOld)%nxOld;
      Lint iyOldTop =((((k/scalX)%(ny))% (ny/scalY))+1+nyOld)%nyOld;
      double q11 = arrayOld[ixOldLeft*ny/scalY+iyOldLow];
      double q21 = arrayOld[ixOldRight*ny/scalY+iyOldLow];
      double q12 = arrayOld[ixOldLeft*ny/scalY+iyOldTop];
      double q22 = arrayOld[ixOldRight*ny/scalY+iyOldTop];
      ixOldLeft  *= scalX;
      ixOldRight *= scalX;
      iyOldLow   *= scalY;
      iyOldTop   *= scalY;
      array[k] = biLinInterpol(q11, q12, q21, q22, ixOldLeft, ixOldRight,
      iyOldLow, iyOldTop, ix, iy);
    } } // end for's
  } else { // 2 cases remaining, to be coded eventually
    cerr << "### nx, ny incompatible with nxOld, nyOld" << endl;
  }
  return(fSizeUpScale);
}

void spectralAnalysis(Complex * f, vector<double> *derivs, Lint nx, Lint ny) {

  Lint nyHP1 = ny/2 + 1;
  double dqx = 2*PI/lengthX, dqy = 2*PI/lengthY;

  double varHeight = 0, gradX2 = 0, gradY2 = 0, curvX2 = 0, curvY2 = 0;

  for (Lint k = 0; k < nx*ny; ++k) {

    double specLoc = f[k].real()*f[k].real() + f[k].imag()*f[k].imag();

    Lint iqx = k/nyHP1, iqy = k%nyHP1;
    int jqx = abs(iqx-nx);
    iqx = (iqx<jqx) ? iqx : jqx;
    double qxN = iqx*dqx, qyN = iqy*dqy;

    double weight = 2;
    if ( (iqy==0) || (iqy==ny/2) ) weight = 1;
    specLoc *= weight;

    // height variance
    varHeight += specLoc;

    // slope variances
    qxN *= qxN; qyN *= qyN;
    gradX2 += qxN * specLoc;
    gradY2 += qyN * specLoc;

    // curvature variances
    qxN *= qxN; qyN *= qyN;
    curvX2 += qxN * specLoc;
    curvY2 += qyN * specLoc;

  }

  (*derivs)[0] = varHeight;
  (*derivs)[1] = gradX2;
  (*derivs)[2] = gradY2;
  (*derivs)[3] = curvX2;
  (*derivs)[4] = curvY2;

}

void realSpaceAnalysis(double * array, vector<double> *derivs, Lint nx, Lint ny) {
  // Calculates surface stats in real space without assuming periodic boundaries
  double msGradRx = 0., msCurvRx = 0.;
  double msGradRy = 0., msCurvRy = 0.;
  double mHeightR = 0., msHeightR = 0.;
  double minHeight = array[0], maxHeight = array[0];
  for (Lint ix=0; ix < nx; ++ix) {
    for (Lint iy=0; iy<ny; ++iy) {

      double specLoc = array[ix*ny+iy];

      if (specLoc < minHeight) minHeight = specLoc;
      else if (specLoc > maxHeight) maxHeight = specLoc;

      mHeightR += specLoc;
      msHeightR += specLoc*specLoc;
      if (ix > 0){
        msGradRx += pow(specLoc - array[(ix-1)*ny+iy],2); // Norm dx missing
        if (ix < nx-1) {
          msCurvRx += pow(array[(ix+1)*ny+iy] + array[(ix-1)*ny+iy] - 2*specLoc, 2); // Norm dx2 missing
        }
      }
      if (iy > 0){
        msGradRy += pow(specLoc - array[ix*ny+iy-1],2); // Norm dy missing
        if (iy < ny-1) {
          msCurvRy += pow(array[ix*ny+iy+1] + array[ix*ny+iy-1] - 2*specLoc, 2); // Norm dy2 missing
        }
      }
  } } // end for's

  // Normalization
  double dx = lengthX/nx, dy = lengthY/ny;
  double dx2 = dx*dx, dy2 = dy*dy;
  mHeightR  /= nx*ny;
  msHeightR /= nx*ny;
  Lint fac = (nx-1)*ny;
  msGradRx /= fac*dx2; // Norm dx
  fac = nx*(ny-1);
  msGradRy /= fac*dy2; // Norm dy
  fac = (nx-2)*ny;
  msCurvRx /= fac*dx2*dx2; // Norm dx2
  fac = nx*(ny-2);
  msCurvRy /= fac*dy2*dy2; // Norm dy2

  double heightDiff = maxHeight - minHeight;

  (*derivs)[0] = msHeightR - pow(mHeightR,2) ;
  (*derivs)[1] = msGradRx;
  (*derivs)[2] = msGradRy;
  (*derivs)[3] = msCurvRx;
  (*derivs)[4] = msCurvRy;
  (*derivs)[5] = heightDiff;

}


// cannot use template for fireRedirect, as complex and real squares differ

void fireRedirecT(Complex *force, Complex *rN, Complex *rC,
		  Complex *rO, double *mass, Lint nR) {

  // rN = rNew, rC = rCopy (of new), rO = rOld

  if(iTime==1) cerr << "# fireRedirecT not well tested" << endl;

  Complex wvDotWV = 0, fDotF = 0; // wv: weighted velocity
  for(Lint k = 0; k < nR; ++ k){
    Complex velocity = rN[k] - rO[k];
    wvDotWV += mass[k] * 1. * velocity   * conj(velocity);
    fDotF += mass[k] * 1. * force[k] * conj(force[k]);
  }

  double scaleV = (1. - fireRedrct);
  double scaleF = fireRedrct * sqrt( wvDotWV.real() / fDotF.real() );

  for(Lint k = 0; k < nR; ++ k){
    rC[k] = rN[k];
    rN[k]  = rO[k] + scaleV * (rN[k] - rO[k]);
    rN[k] += scaleF * force[k];
  }

}
