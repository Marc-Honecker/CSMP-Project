class atomicSheet{

  private:

  int ID;
  int nAtom;
  static const int nAtomSheetMax = 2;
  int atomSheet[nAtomSheetMax], nAtomSheet;

  static const int nDim = 3;
  vector<vector<double> > rNow, rOld, rOldCopy, force;

  double mass, sigma, epsAtomWall;

  double tKinetic, tKin1, tKin2;

  void initParamsDefault();
  void writeParamsDefault();

  public:

  void initParams(int);
  void initSystem();
  void writeParams();
  void dumpConfig();
  bool readConfig();
  void initMeasure();
  double propagate();

};
