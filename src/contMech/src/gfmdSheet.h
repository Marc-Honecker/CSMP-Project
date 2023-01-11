//AW2CM: all the modifications can be tracked with "AW modified:"
class gfmdSheet{

  private:

  int ID;

  double areaPP;
  double zeroModeMass;

  double pressure, pressInit, pressFinal;

  int fConstCOM; // = 1 controlled COM
  int fSteppedRamp, rampSteps, rampRelax;
  double zConstCOM, vzConstCOM, stressCOM, stressInter, dzRamp, ddzRamp, dpRamp, ddpRamp;
  double rVeloTurnStep; // relative number of steps before inversion of ramp
  int nVeloTurnStep;	// absolute number of steps ...
  int nVeloTransition, tStartTransition, tEndTransition, t0Transition;
  double vzInit;

  int fzOpposite;
  double zOpposite;

  double vX, vY;
  double dxOffset, dyOffset;

  int fOnSitePotential;
  double fOnSitePeriod, fOnSiteFreq;
  int fKelvinVoigt;
  double tauKV, scalKV;

  double rXhertz, rYhertz, hertzExpnt;
  double rSphere;

  double hurst, lambdaR, lambdaS, rRoughNorm, peklenik;
  int fRollOff, fRoughNorm, fBoxMuller;

  int fRoughRead, fRoughAdd;

  int fAddSWR, nqxAddSWR, nqyAddSWR; // SWR = single-wavelength roughness
  double heightSWR;

  int fFlatPunch;
  double radiusFlatPunch, heightFlatPunch;

  double dzStep;

  // space-holder fields to be used in FFTW
  double *fieldRFFT;
  Complex *fieldFFFT;

  Lint kFastest, rCentral;

  int rampTime;
  double damping; 

  string konfigName;
  int f3dMovie, resolMovie;

  ofstream moni, ramp;
  ofstream traj;//ENH-Maxwell

  void initParamsDefault();
  void initParamsFile();
  void writeParamsDefault();

  void addConfigR();
  void addConfigF();

  void initConfig();
  void addHertz();
  void addSphere();
  void addSelfAffine(), selfAffineAnalysis(vector<double> *derivs);
  void addSWR();
  void addFlatPunch();

  void dumpConfig();
  void makeSteps();//CM-dzSteps
  double sqrtSpec(double);

  void dispZ_prepare();
  void initLateral();
  void moveLateral();

  double deltaIndenter(), einsteinSolid();
  double stressGFMD();
  double stressKV();
  void thermostat();

  void moveCOM();
  void moveOpposite();

  double computeChi2SS();
  void propagateSteadySlide(); 

  // divisions and modulos could be optimized with bit-shift operations!
  Lint irx, iry, iqx, iqy;
  inline Lint getLinRIndex(Lint ix, Lint iy) {return(ix*ny+iy);};
  inline void get_irxiry(Lint k){ irx = k/ny; iry = k%ny; };
  inline Lint getLinCIndex(Lint ix, Lint iy) {return(ix*nyHP1+iy);};
  inline void get_iqxiqy(Lint k){ iqx = k/nyHP1; iqy = k%nyHP1; };

  double getQ(Lint);
  
  //AW modified:
  int fMaxwell, nMaxwell, fMwRelax;
  int nMuTiSt; //multiple time step option
  int fReadMwElement; //read prony parameters from maxwell.in
  
  double dTimeMw;

  //addition field for spring and damper of maxwell model
  vector<double> kMaxwell, dampMaxwell;

  //functions to propagate the model
  void initMaxwell();
  void addStressMaxwell();
  
  //file stream to read in maxwell.in
  ifstream readProny;

  //AW modified end

  public:

  Lint nx, ny;
  Lint nxH, nyHP1, nxny, nReal, nFour;

  double dx, dqx, dqx2, dy, dqy, dqy2, q2Max;
  
  //addedAnle:
  double massGFMD; 
  
  int fRough, nElast;
  int fMassWeightg;
  int fLateral, fSteadySlide;
  bool fSliding;
  double chi2SteadySlide, chi2SteadySlideOld, weightSS, invStiffSS_COM, invStiffSS_COMRef, chi2SSRestart, prefacSSnxny;
  int oldSSPress;

  int fSheetType; //	if rough: += 1; if nElast: += 2; makes 3 different types

  // energy: v = stiffness * q^elastExpnt * u^2 / 2
  // contactMod = 2 * stiffness; ==> default: stiffness = 0.5
  static const int nLayer=4;
  double stiffness[nLayer], elastExpnt[nLayer];
  double stiffHigh[nLayer];//ENH-Maxwell
  double poisson[nLayer], thickness[nLayer];
  int fThickness[nLayer];
  double stiffMax, stiffMin, contModEff, massScal;

  double tKinHalfStep();
  double tKinetic, vPotTot, tKinetic1, tKinetic2;
  double vOnSitePot, vElastic, vGravit;

  // config fields
  double *equilPos;
  Complex *equilPosF, *equilPos0F;

  double  *dispZR, *stressZR, rKV_LPF;
  Complex *dispZF, *dispZFold, *stressZF, *rhsKV_LPF; 

  short int *weightQ;
  double *massCorr, *stressZFPre;
  Complex *preFacCorrSS;
  
  //ENH-Maxwell
  vector<Complex*> dispZMw, stressZMw; 

  // fftw_plans
  fftw_plan dispZR2F, dispZF2R, stressZR2F;
  fftw_plan equilPosF2R;

  // called from contMech
  void initParams(int);
  void initSystem();
  void writeParams();
  void initMeasure();
  double propagate();
  double propagateKV();
  void zeroStress();
  double addStress();
  void measure();
  void fireHalt();
  void fireRedirect();
  void outMeasure();
  void outSystem();
  void dumpFrame(int);
  ~gfmdSheet();

  // called / set from inter
  bool fRealSpaceStress;
  void dispF2R();
  void dispR2F();

};
