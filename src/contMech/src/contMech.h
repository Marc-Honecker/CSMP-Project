// functions called in main

void initParams();
void initSystem();
void writeParams();
void initMeasure();
void propagate();
void constrain();
void getForces();
void fire(); 
void measure();
bool finished();
void outMeasure();
void outSystem();

// additional functions

void initParamsDefault();
void initParamsFile();
void writeParamsDefault();

// global parameters

double lengthX, lengthY, areaXY;
Lint nxGlobal, nyGlobal;

Lint iTime, nTime, nRelax;
double dTime, dTimeInit, dTime2, dampGlobal, mdTime;

int fFire, fFireOn, iTimeLFS; // LFS = last fire stop
double fireRedrct, fireIncrmt, fireDecrmt, fireRedrctInit;
bool fFireConstraint;

int fLangevin;
double temp, tempInit,  tempFinal;
bool fExternalDriving;

int nSheet, nInter, nAtomL;

ifstream readParams;
ofstream moni;

int iFrame, freqFrame;

double timeTotal;

// parameters only needed in condMat

int randSeed;

double tKinGlobal, vPotGlobal;
double tKinOld, vPotOld;

time_t timeLocal;
double tPropagate, tGetForces, tConstrain, tFire, tMeasure;
