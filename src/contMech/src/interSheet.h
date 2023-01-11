
class interSheet{

  private:

  int ID;
  int fInterType;
  int fFirstCurve;

  int sheetID0, sheetID1;
  Lint nx, ny;
  Lint nxH, nyHP1, nxny, nReal, nFour;
  Lint iqx, iqy;
  void get_iqxiqy(Lint);

  ofstream moni;

  // space-holder fields to be used in FFTW
  double *fieldRFFT;
  Complex *fieldFFFT;
  // fftw_plan
  fftw_plan gapR2F, gapGradF2R;

  void initParamsDefault();

  double areaPP; // area per particle

  void (interSheet::*makeGapX)();
  void makeGap1(), makeGap2(), makeGap3(), makeGap4(), makeGap5();

  void (interSheet::*constraintX)();
  void constraint1(), constraint2(), constraint3(), constraint4(), constraint5();

  int fDumpGap, fDumpLateral, fDumpFrame;
  int resolMovie;

  void makeGradientPotential();
  double (interSheet::*potentialX)(double, double*);
  double potential1(double, double*), potential2(double, double*), potential3(double, double*);
  void potentialInit();
  void potentialTest();

  double fullContactElastEnergy();
  void dumpGap();

  double contactArea, attractArea, meanGap, meanPosGap;
  double surfEnerg, potCurveRel, potRange, potCurve;
  double interStiffMax;
  vector <double> potParam;

  double lateralForce();

  int ixMaxTension; // change
  int nTimeOff;

  public:

  void makeGap(), constraint();
  double addStress();
  void initMeasure();
  void measure();
  void dumpFrame(int);
  void outMeasure();

  int fConstraint, fPotential, fReynolds, fElastic, fShifted;
  int fPotentialTest;

  double vPotTot;

  int nGap;
  double *gap, *gapStress;

  void initParams(int), initSystem();
  void writeParams(), writeParamsDefault();
  void outSystem();
  
  ~interSheet();

};
