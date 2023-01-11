#include "header.h"
#include "contMech.h"
#include "gfmdSheet.h"
#include "interSheet.h"
#include "atomicSheet.h"

const int nSheetMax=3, nInterMax=3, nAtomLMax=2;
vector <gfmdSheet> sheet(nSheetMax);
vector <interSheet> inter(nInterMax);
vector <atomicSheet> atomL(nAtomLMax);

extern void termination(const string&);
extern double rrand();
extern void changeSeedMix64(int);

int main(){

  initParams();
  writeParams();
  initSystem();
  initMeasure();

  while(!finished()){
    iTime += 1;
    if(iTime) propagate();
    constrain();
    getForces();
    if(fFire) fire();
    measure();
  }

  outMeasure();
  outSystem();

}

void initParams(){

  // setting parameters
  initParamsDefault();
  readParams.open("params.in");
  if ( readParams.is_open() ) initParamsFile();
  else cerr << "# Using  global default parameters\n";
  readParams.close();
  
  // further initializations
  iTime = -1;
  temp = tempInit;
  srand(randSeed);
  changeSeedMix64(randSeed);

  // sanity checks and parameter post processing

  // http://www.graphics.stanford.edu/~seander/bithacks.html
  if ( ((nxGlobal&(nxGlobal-1))!=0) || ((nyGlobal&(nyGlobal-1))!=0) )
  cerr << "\n### nxGlobal and/or nyGlobal is not a power of 2. Caution!!!\n\n";

  if ( (lengthX<=0) || (lengthY<=0) )
  termination("lengthX and/or lengthY not positive.\n");

  if (fFire&&fLangevin) termination("fFire + fLangevin used simultanseously.\n");

  areaXY = lengthX*lengthY;
  dTimeInit = dTime;
  dTime2 = dTime*dTime;

  vPotGlobal = tKinGlobal = 0; //CM-change: if we ever read in old velocities...

  // initialize object parameters
  for (int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].initParams(iSheet);
  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].initParams(iInter);
  for (int iAtomL = 0; iAtomL < nAtomL; ++iAtomL) atomL[iAtomL].initParams(iAtomL);

}

void initSystem(){
  for (int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].initSystem();
  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].initSystem();
  for (int iAtomL = 0; iAtomL < nAtomL; ++iAtomL) atomL[iAtomL].initSystem();
}

void writeParams(){

  writeParamsDefault();
  ofstream output("params.out");

  output << lengthX     << "\t\t# lengthX #" << "\n";
  if (lengthY!=lengthX) output << lengthY     << "\t\t# lengthY #" << "\n\n";

  output << nxGlobal    << "\t\t# nxGlobal #" << "\n";
  if (nyGlobal!=nxGlobal) output << nyGlobal    << "\t\t# nyGlobal #" << "\n\n";

  output << "\n"<<nRelax << "\t\t# nRelax #" << "\n";
  output << nTime << "\t\t# nTime #" << "\n";
  output << dTime << "\t\t# dTime #" << "\n\n";

  output << dampGlobal  << "\t\t# dampGlobal #" << "\n";
  output << randSeed    << "\t\t# randSeed #" << "\n\n";

  if (freqFrame) output << freqFrame << "\t\t# freqFrame #" << "\n";

  if (fFire) {
    output << fFire       << "\t\t# fFire #" << "\n";
    output << fireRedrct  << "\t\t# fireRedrct #" << "\n";
    output << fireIncrmt  << "\t\t# fireIncrmt #" << "\n";
    output << fireDecrmt  << "\t\t# fireDecrmt #" << "\n\n";
  }

  if (fLangevin) {
    output << fLangevin << "\t\t# fLangevin #" << "\n";
    output << tempInit  << "\t\t# tempInit #"  << "\n";
    if (tempInit!=tempFinal) output << tempFinal  << "\t\t# tempFinal #" << "\n\n";
  }

  output << nSheet      << "\t\t# nSheet #" << "\n";
  output << nInter      << "\t\t# nInter #" << "\n";
  output << nAtomL      << "\t\t# nAtomL #" << "\n\n";
  output << 0           << "\t\t# end global parameters \n\n";
  output.close();

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].writeParams();
  for (int iInter=0; iInter<nInter; ++iInter) inter[iInter].writeParams();
  for (int iAtomL=0; iAtomL<nAtomL; ++iAtomL) atomL[iAtomL].writeParams();

} // writeParams

void initMeasure(){

  mdTime = -dTime*nRelax;

  moni.open("gMoni.dat", ofstream::out);
  timeTotal = clock();

  tPropagate = tGetForces = tConstrain = tMeasure = 0;

  for(int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].initMeasure();
  for(int iInter = 0; iInter < nInter; ++iInter) inter[iInter].initMeasure();
  for(int iAtomL = 0; iAtomL < nAtomL; ++iAtomL) atomL[iAtomL].initMeasure();

  if(freqFrame) {
    cerr << "# Animating " << nTime/freqFrame << " frames.\n";
    system("mkdir -p frames");
    system("rm -f frames/*");
  }

}

void propagate(){

  timeLocal = (double) clock();

  if (fFire) {
    if (vPotGlobal>vPotOld) {
      dTime *= fireDecrmt;
      for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].fireHalt();
      
    } else {
      dTime *= fireIncrmt;
      if (fFireConstraint && (dTime>2*dTimeInit) ) dTime = 2*dTimeInit;
      if (dTime>2*dTimeInit) dTime = 2*dTimeInit; // MM change
    }
  }

  mdTime += dTime;

  temp = tempInit;
  if (iTime>nRelax) temp += (iTime-nRelax)*(tempFinal-tempInit)/nTime;

  tKinGlobal = 0;

  for (int iSheet=0; iSheet<nSheet; ++iSheet) 
  tKinGlobal += sheet[iSheet].propagate();

  for (int iAtomL=0; iAtomL<nAtomL; ++iAtomL) 
  tKinGlobal += atomL[iAtomL].propagate();

  tPropagate += (double) clock() - timeLocal;

}

void constrain(){
  timeLocal = (double) clock();

  for (int iInter = 0; iInter < nInter; ++iInter)
  if(inter[iInter].fConstraint) inter[iInter].constraint();

  tConstrain += (double) clock() - timeLocal;
}

void getForces(){ 

  timeLocal = (double) clock();

  if(fFire) fFireOn=1;
  vPotGlobal = 0;

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].zeroStress();

  for (int iInter=0; iInter<nInter; ++iInter) 
  vPotGlobal += inter[iInter].addStress();

  for (int iSheet=0; iSheet<nSheet; ++iSheet) 
  vPotGlobal += sheet[iSheet].addStress();

  tGetForces += (double) clock() - timeLocal;

}

void fire(){
  // recompute kinetic energy on half steps
  tKinGlobal = 0;
  for (int iSheet = 0; iSheet<nSheet; ++iSheet) 
  tKinGlobal += sheet[iSheet].tKinHalfStep();

  if(tKinOld>0.999*tKinGlobal) {
    dTime *= fireDecrmt;
    for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].fireHalt();
    tKinGlobal = 0;
    vPotGlobal = vPotOld;
  }

}

void measure(){

  timeLocal = (double) clock();

  if (iTime==0) moni << "# mdTime\ttKinGlobal\tvPotGlobal\n";
  moni << mdTime << "\t" << tKinGlobal << "\t" << vPotGlobal;
  if(fFire) moni << "\t" << dTime;
  moni << endl;

  vPotOld = vPotGlobal;
  tKinOld = tKinGlobal;

  if(iTime<nRelax) return;

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].measure();

  if( (freqFrame) && !(iTime%freqFrame) ) {
    ++iFrame;
    for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].dumpFrame(iFrame);
    for (int iInter=0; iInter<nInter; ++iInter) inter[iInter].dumpFrame(iFrame);
  }

  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].measure();

  tMeasure += (double) clock() - timeLocal;
  
}

bool finished(){ 

  // planned stop
  if (iTime==(nTime+nRelax)) return(true);

  // fFire stop
  if ( fFire && (dTime<1e-10*dTimeInit) ) {
    cerr << "# Fire: dTime small at iTime = " << iTime << "\n";
    getForces();
    return(true);
  } else if ( fFire && (dTime>1e+100*dTimeInit) ) {
    cerr << "# Fire: dTime large at iTime = " << iTime << "\n";
    getForces();
    return(true);
  }

  // hard stop
  ifstream testStop1("stop");
  if (testStop1.is_open()) {
    system("rm stop");
    testStop1.close();
    termination("stop requested.");
  } testStop1.close();

  // early stop
  ifstream testStop2("finish");
  if (testStop2.is_open()) {
    system("rm finish");
    testStop2.close();
    cerr << "\n# finish requested at iTime = " << iTime << "\n\n";
    return(true);
  } testStop2.close();

  return(false);
}

void outMeasure(){

  timeTotal = (double) clock() - timeTotal;
  double timeKnown = tPropagate + tGetForces + tConstrain + tFire + tMeasure; 

  ofstream timing("params.out",ofstream::app);
  timing << timeTotal/1.e6/60 << "\t# absolute computing time (min)\n";
  timing << timeKnown/1.e6/60 << "\t# measured computing time (min)\n\n";

  timing.precision(3);
  timing << "0 \t# relative computing times:\n" << fixed;

  timing << tPropagate/timeKnown << "\t\t# tPropagate\n";
  timing << tGetForces/timeKnown << "\t\t# tGetForces\n";

  tConstrain /= timeKnown;
  if(tConstrain>0.001) timing << tConstrain << "\t\t# tConstrain\n";

  tFire /= timeKnown;
  if(tFire>0.001) timing << tFire << "\t\t# tFire\n";

  tMeasure /= timeKnown;
  if(tMeasure>0.001) timing << tMeasure << "\t\t# tMeasure\n";

  timing << 0 << "\t# relative computing times:\n\n";
  timing.close();

  // up-date of real-space stress field required
  getForces();  //MM2XY: Why is the elastic stress (sometimes) destroyed in konfig-D.dat?

  for (int iSheet=0; iSheet<nSheet; ++iSheet) sheet[iSheet].outMeasure();
  for (int iInter=0; iInter<nInter; ++iInter) inter[iInter].outMeasure();

} // outMeasure

void outSystem(){

  for (int iSheet = 0; iSheet < nSheet; ++iSheet) sheet[iSheet].outSystem(); 
  for (int iInter = 0; iInter < nInter; ++iInter) inter[iInter].outSystem(); 

  ofstream output("params.out",ofstream::app);
  output << 1 << "\t\t# D O N E";
  output.close(); 

  moni.close();
}

void initParamsDefault(){

  lengthX = 1.; lengthY = 1.;
  nxGlobal = 128; nyGlobal = 128;

  nTime = 100; dTime = dTimeInit = 0.25; dampGlobal = 1.5;
  nRelax = 0;
  randSeed = 4711;

  fFire = 0;
  fireRedrct = 0.0;
  fireIncrmt = 1.2;
  fireDecrmt = 0.5;
  fFireConstraint = false;

  fLangevin = 0;
  tempInit  = 0.01;
  tempFinal = tempInit;

  freqFrame = 0;
  iFrame = 0;

  fExternalDriving = false;

  nSheet = 2;
  nInter = 1;
  nAtomL = 0;

}

void initParamsFile(){

  double param;
  std::string ROL; // rest of line
  std::size_t NIS = std::string::npos; // NIS == Not In String
  int fGlobalRead = 1;
  int iCount = 0; 
  while ( !readParams.eof() && fGlobalRead) {
    readParams >> param; getline(readParams,ROL);
    if(++iCount>40) break;
    if      (ROL.find("# lengthX #")  !=NIS){lengthX = param;  lengthY = lengthX;}
    else if (ROL.find("# lengthY #")  !=NIS) lengthY = param;  
    else if (ROL.find("# nxGlobal #") !=NIS){nxGlobal = param; nyGlobal = nxGlobal;}
    else if (ROL.find("# nyGlobal #") !=NIS) nyGlobal = param; 
    else if (ROL.find("# nRelax #")   !=NIS) nRelax = param;
    else if (ROL.find("# nTime #")    !=NIS) nTime = param;
    else if (ROL.find("# dTime #")    !=NIS) dTime = param;
    else if (ROL.find("# dampGlobal #") !=NIS) dampGlobal = param;
    else if (ROL.find("# randSeed #") !=NIS) randSeed = param;
    else if (ROL.find("# fFire #")    !=NIS) fFire = param;
    else if (ROL.find("# fireRedrct #") !=NIS) fireRedrct = param;
    else if (ROL.find("# fireIncrmt #") !=NIS) fireIncrmt = param;
    else if (ROL.find("# fireDecrmt #") !=NIS) fireDecrmt = param;
    else if (ROL.find("# fLangevin #")!=NIS) fLangevin = param;
    else if (ROL.find("# tempInit #") !=NIS) {tempInit = param; tempFinal = tempInit;}
    else if (ROL.find("# tempFinal #")!=NIS) tempFinal = param;
    else if (ROL.find("# freqFrame #")!=NIS) freqFrame = param;
    else if (ROL.find("# nSheet #")   !=NIS) nSheet = param;
    else if (ROL.find("# nInter #")   !=NIS) nInter = param;
    else if (ROL.find("# nAtomL #")   !=NIS) nAtomL = param;
    else if (ROL.find("# end global") !=NIS) fGlobalRead = param;
  }

}

void writeParamsDefault(){

  ofstream output("params.def");
  output << "! ! ! ! ! ! ! ! ! global (default) values\n\n";
  
  output << 1 << "\t\t# dampGlobal #" << "\n";
  output << 4711 << "\t\t# randSeed #" << "\n\n";

  output << 1 << "\t\t# fFire #" << "\n";
  output << 1e-3  << "\t\t# fireRedrct #" << "\n";
  output << 1.2  << "\t\t# fireIncrmt #" << "\n";
  output << 0.5  << "\t\t# fireDecrmt #" << "\n\n";

  output << 1 << "\t\t# fLangevin #" << "\n";
  output << 0.01  << "\t\t# tempInit #" << "\n";
  output << 0.01 << "\t\t# tempFinal #" << "\n\n";

  output << 100 << "\t\t# freqFrame #" << "\n\n";

  output.close();

}
