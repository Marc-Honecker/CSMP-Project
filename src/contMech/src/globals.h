
// global variables from contMech

extern Lint nxGlobal, nyGlobal;
extern double lengthX, lengthY, areaXY;

extern Lint iTime, nTime, nRelax;
extern double dTime, dTime2, dampGlobal, mdTime;

extern int fFire, fFireOn;
extern double fireRedrct, fireIncrmt, fireDecrmt, fireRedrctInit;

extern int fLangevin;
extern double temp, tempInit,  tempFinal;

extern bool fExternalDriving;

extern int nSheet, nInter;

extern ifstream readParams;

// global functions from auxiliary

extern double rrand();

extern void termination(const string&);
extern void termination(const string&, int);
extern void termination(const string&, double);

extern void dumpRealCS(vector<double*>, vector<string>, Lint, Lint, int, int, const string&);
extern void dumpReal(vector<double*>, Lint, Lint, int, int, const string&);
extern int  readReal(double*, Lint, Lint, const string&, int);

extern void spectralAnalysis(Complex *, vector <double> *, Lint nx, Lint ny);
extern void realSpaceAnalysis(double *, vector<double> *, Lint, Lint);
extern void fireRedirecT(Complex *, Complex *, Complex *, Complex *, double *, Lint);

extern double mix64(void);
extern void changeSeedMix64(int);
