// MM2CM: http://gnuplot.sourceforge.net/docs_4.2/node103.html
// ... zu bedenken, wenn wieder (2D statt 1D) Movies gemacht werden

// MM2CM: flag for C33 or E_Young coupling of COM mode

// MM2All: Think about nVeloStopStep, nVeloStartStep

#include "header.h"
#include "gfmdSheet.h"
#include "globals.h"

// functions called from main

void gfmdSheet::initParams(int newID){

  if (newID>1) cerr << "# ID might be too large in gfmdSheet.";
  ID = newID;

  // setting parameters
  initParamsDefault();
  readParams.open("params.in");
  if ( readParams.is_open() ) initParamsFile();
  readParams.close();

  //AW modified:
  if (fMaxwell) initMaxwell(); 

  // sanity checks and post-processing of parameters

  nxH = nx/2; nyHP1 = (ny/2) + 1;
  nFour = nx*nyHP1;
  nReal = 2*nFour;
  nxny = nx*ny;
  areaPP = areaXY/nxny;

  dx = lengthX / nx; dqx = TWOPI / lengthX; dqx2 = dqx*dqx;
  dy = lengthY / ny; dqy = TWOPI / lengthY; dqy2 = dqy*dqy;

  double q2Max = 0;
  if (nx>1) q2Max += pow(PI/dx,2);
  if (ny>1) q2Max += pow(PI/dy,2);
  if (q2Max==0) termination("q2Max not supposed to be zero.");
  double qMax = sqrt(q2Max);

  // do post processing of roughness parameters first

  fRough = 0;

  // add 4 to fRoughAdd, if this bit is not set in fRoughAdd
  if( fAddSWR && !(fRoughAdd&4) ) fRoughAdd += 4;

  if( fRoughAdd || fRoughRead ) fRough = 1;


  fSheetType = 0;
  if(fRough) fSheetType += 1;
  if(nElast) fSheetType += 2;

  if ( rRoughNorm<=0 ) {
    rRoughNorm = 1;
    cerr << "# Forcing rRoughNorm=1" << endl;
  }

  if (peklenik<=0) {
    cerr << "# Forcing peklenik=1" << endl;
    peklenik = 1;
  }

  if ( fSteppedRamp && (vzConstCOM!=0) ) {
    cerr << "# Setting vzConstCOM=0 due to fSteppedRamp" << endl;
    vzConstCOM = 0;
  }

  if (fFire&&fKelvinVoigt) termination("fFire AND fKelvinVoigt not allowed!");

  if(fKelvinVoigt) {
    if(fMassWeightg) termination("fMassWeightg AND fKelvinVoigt not allowed!");
    if(rKV_LPF < 0) rKV_LPF = 0;
  }

  konfigName = "konfig" + to_string(ID);
  if(fRough) konfigName = konfigName + "E"; 
  if(nElast) konfigName = konfigName + "D"; 
  konfigName = konfigName + ".";

  // finish function with processing of elasticity-related parametrizations
  if(!nElast) return;

  if ( (fSteppedRamp==1) && (!fConstCOM) ) fConstCOM = 1; 

  if (zOpposite) zOpposite *= nxny;

  if ( (fConstCOM&&fFire) && (vzConstCOM!=0) ) {
    cerr << "# Switching off FIRE becaucse of vzConstCOM!=0\n";
    fFire = 0;
  }

  if ( (fSteppedRamp)&&(fFire) ) {
    if ( (fireIncrmt!=1)||(fireDecrmt!=1)||(fireRedrct!=0) )
    cerr << "# Resetting FIRE variables because of ramp, which may be unnecessary.\n";
    fireIncrmt = 1;
    fireDecrmt = 1;
    fireRedrct = 0;
  }

  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    if ( (elastExpnt[iLayer]!=1) && (fThickness[iLayer]!=0) )
    termination("Exponent of layer violates thickness.");
  }

  if (fSteppedRamp==1) {
    ddzRamp = nxny * dzRamp / rampSteps;
    pressInit = pressFinal = 0;
  } else if (fSteppedRamp==2) {
    ddpRamp = dpRamp / rampSteps;
  }

  pressure = pressInit;

  if (rVeloTurnStep>0) nVeloTurnStep = rVeloTurnStep*nTime;
  if ((nVeloTurnStep>0) && (nVeloTransition>0)) {
    tStartTransition = nVeloTurnStep - nVeloTransition/2;
    tEndTransition = nVeloTurnStep + 3*nVeloTransition/2;
    t0Transition = nVeloTurnStep + nVeloTransition/2;
  }
  else {
    nVeloTransition = -1;
    tStartTransition = -1;
    tEndTransition = -1;
  }
  vzInit = vzConstCOM;


  if (fOnSitePotential!=0) fRealSpaceStress = 1;
  if ((fOnSitePeriod!=0)&&(fOnSitePotential==1)) fOnSiteFreq = 2*PI/fOnSitePeriod;

  if(scalKV<1+1e-5) scalKV = 1+1e-5;

  if(tauKV<=1e-10) cerr << "# Check units of tauKV!\n";

  double q2Min = min(dqx2,dqx2);
  if(nx==1) q2Min = dqy2;
  if(ny==1) q2Min = dqx2;
  double qMin = sqrt(q2Min);

  // stiffMin is smallest non-zero elastic stiffness
  stiffMax = contModEff = stiffMin = 0;
  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    stiffMax += stiffHigh[iLayer]*pow(qMax,elastExpnt[iLayer]);//ENH-Maxwell
    stiffMin += stiffness[iLayer]*pow(qMin,elastExpnt[iLayer]);
  }
  contModEff = 2*stiffMax / sqrt(q2Max);
  //change: stiffMax might be for q=0 if fThickness==2.

  // assign masses and damping (unit of damping is inverse seconds)
  massScal = stiffMax;
  double stiffEff = sqrt(0*pressure*pressure/areaXY+stiffMin*stiffMin); 
  damping = dampGlobal*sqrt(stiffEff/massScal); 


} // initParams

// function called from within object

void gfmdSheet::initParamsDefault(){

  nx = nxGlobal;
  ny = nyGlobal;

  fRoughAdd = fRoughRead = nElast = 0;

  if(ID==0) fRoughAdd = 1;
  if(ID==1) nElast = 1;
  if (nSheet==1) {fRoughAdd = 0;}
  else if (nSheet!=2) cerr << "# Be aware of default setting for nSheet = " << nSheet << endl;

  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    stiffness[iLayer] = stiffHigh[iLayer] = 0;
    elastExpnt[iLayer] = 1;
    poisson[iLayer] = 0.25;
    fThickness[iLayer] = 0;
    thickness[iLayer] = lengthY;
  } stiffness[0] = 0.5; // makes contact modulus of regular sheet = 1

  fMassWeightg = 0;
  zeroModeMass = 1;
  massGFMD = 1;

  // default values for elastomer
  pressure = pressInit = pressFinal = 0.01;

  fConstCOM = 0;	// 1 (2) : constr velocity of top (bottom) surface
  zConstCOM = 0; 
  vzConstCOM = 0;
  fSteppedRamp = 0;
  dzRamp = 1.e-2;
  rampSteps = 130;
  rampRelax = 170;
  rVeloTurnStep = -1;
  nVeloTurnStep = -1;
  nVeloTransition = 0;

  fzOpposite = 0;  
  zOpposite = 0;  

  fLateral = vX = vY = 0;
  fSteadySlide = false;

  fOnSitePotential = 0;
  fOnSitePeriod = 0;

  fKelvinVoigt = 0;
  tauKV = 1.;
  scalKV = 1000.;
  rKV_LPF = 0;
  
  //AW modified:
  fMaxwell = 0;
  nMaxwell = 2;
  fMwRelax = 1;
  nMuTiSt = 1;
  fReadMwElement = 0; //not needed now
  //dTimeMw = 0.025;//ENH-Maxwell

  // default roughness values

  // Hertz roughness 	fRoughAdd&1
  rXhertz = rYhertz = 1;
  hertzExpnt = 2;

  // random roughness	fRoughAdd&2
  hurst = 0.8;
  lambdaR = 0.5;
  lambdaS = 0.02;
  fRollOff = 1;
  fRoughNorm = 1;	// 1: normalize to rms gradient; 2: ... to rms height
  rRoughNorm = 1.;
  peklenik = 1.;
  fBoxMuller = 0;

  // SWR: single-wavelength roughness	fRoughAdd&4
  fAddSWR = 0;
  nqxAddSWR = 2; nqyAddSWR = nqxAddSWR;
  heightSWR = 0.1;

  // flat punch:	fRoughAdd&8
  radiusFlatPunch = 0.2;
  heightFlatPunch = 1;

  // hemisphere:	fRoughAdd&16
  rSphere = 0.2;

  // discrete steps
  dzStep = -1./200;

  fRealSpaceStress = 0;

  // movie parameters
  f3dMovie = 0;
  resolMovie = 512; // max number of points per line
}

void gfmdSheet::initParamsFile(){

  // overwrite following default setting
  nElast = fRoughAdd = fAddSWR = 0;

  int fReadParams = 0;
  while ( !readParams.eof() ) {

    double param;
    std::string ROL; // rest of line
    std::size_t NIS = std::string::npos; // NIS == Not In String
    if (readParams.eof()) break;
    readParams >> param; getline(readParams,ROL);

    if ( (param==ID) && (ROL.find("# sheet start") !=NIS) ) fReadParams = 1;
    if (fReadParams==0) continue;
    if (ROL.find("# sheet end")   !=NIS) break;

    if (ROL.find("# fRoughAdd #") !=NIS) fRoughAdd = param;
    if (ROL.find("# fRoughRead #") !=NIS) fRoughRead = param;
    if (ROL.find("# rXhertz #") !=NIS) {rXhertz = param; rYhertz = param;}
    if (ROL.find("# rYhertz #") !=NIS) rYhertz = param;
    if (ROL.find("# hertzExpnt #") !=NIS) hertzExpnt = param;
    if (ROL.find("# rSphere #") !=NIS) rSphere = param;
    if (ROL.find("# hurst #") !=NIS) hurst = param;
    if (ROL.find("# lambdaR #") !=NIS) lambdaR = param;
    if (ROL.find("# lambdaS #") !=NIS) lambdaS = param;
    if (ROL.find("# rRoughNorm #") !=NIS) rRoughNorm = param;
    if (ROL.find("# peklenik #") !=NIS) peklenik = param;
    if (ROL.find("# fRollOff #") !=NIS) fRollOff = param;
    if (ROL.find("# fRoughNorm #") !=NIS) fRoughNorm = param;
    if (ROL.find("# fBoxMuller #") !=NIS) fBoxMuller = param;
    if (ROL.find("# fAddSWR #") !=NIS) fAddSWR = param;
    if (ROL.find("# nqxAddSWR #") !=NIS) {nqxAddSWR = param; nqyAddSWR = nqxAddSWR;}
    if (ROL.find("# nqyAddSWR #") !=NIS) nqyAddSWR = param;
    if (ROL.find("# heightSWR #") !=NIS) heightSWR = param;
    if (ROL.find("# radiusFlatPunch #") !=NIS) radiusFlatPunch = param;
    if (ROL.find("# heightFlatPunch #") !=NIS) heightFlatPunch = param;
    if (ROL.find("# dzStep #") !=NIS) dzStep = param;
    if (ROL.find("# nElast #") !=NIS) nElast = param;
    for (int iLayer=0; iLayer<nLayer; ++iLayer) {
      std::string text = "# stiffness" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) stiffness[iLayer] = param;
      text = "# stiffHigh" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) stiffHigh[iLayer] = param; //ENH-Maxwell
      text = "# elastExpnt" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) elastExpnt[iLayer] = param; 
      text = "# fThickness" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) fThickness[iLayer] = param; 
      text = "# poisson" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) poisson[iLayer] = param; 
      text = "# thickness" + to_string(iLayer) + " #";
      if (ROL.find(text) != NIS) thickness[iLayer] = param; 
    }
    if (ROL.find("# fMassWeightg #") !=NIS) fMassWeightg = param;
    if (ROL.find("# zeroModeMass #") !=NIS) zeroModeMass = param;
    if (ROL.find("# fKelvinVoigt #") !=NIS) fKelvinVoigt = param;
    
    //AW modified:
    if (ROL.find("# fMaxwell #") !=NIS) fMaxwell = param;
    if (ROL.find("# nMaxwell #") !=NIS) {
      nMaxwell = param;
      kMaxwell.resize(nMaxwell,0); dampMaxwell.resize(nMaxwell,0);
    }
    if (ROL.find("# massGFMD #") !=NIS) massGFMD = param;
    if (ROL.find("# fMwRelax #") !=NIS) fMwRelax = param; 
    if (ROL.find("# nMuTiSt #")  !=NIS)  nMuTiSt = param;
    //only read in the maxwell.in when we set the fReadMwElemnt = 1;
    if (ROL.find("# fReadMwElement #")!=NIS) fReadMwElement = param;
    //if (ROL.find("# dTimeMw #")!=NIS) dTimeMw = param;//ENH-Maxwell
    //AW modified end;

    if (ROL.find("# tauKV #") !=NIS) tauKV = param;
    if (ROL.find("# scalKV #") !=NIS) scalKV = param;
    if (ROL.find("# rKV_LPF #") !=NIS) rKV_LPF = param;
    if (ROL.find("# pressInit #") !=NIS) {pressInit  = param; pressFinal=param;}
    if (ROL.find("# pressFinal #") !=NIS) pressFinal = param;
    if (ROL.find("# fOnSitePotential #") !=NIS) fOnSitePotential = param;
    if (ROL.find("# fOnSitePeriod #") !=NIS) fOnSitePeriod = param;
    if (ROL.find("# fConstCOM #") !=NIS) fConstCOM = param;
    if (ROL.find("# zConstCOM #") !=NIS) zConstCOM = param;
    if (ROL.find("# vzConstCOM #") !=NIS) vzConstCOM = param;
    if (ROL.find("# fSteppedRamp #") !=NIS) fSteppedRamp = param;
    if (ROL.find("# rampSteps #") !=NIS) rampSteps = param;
    if (ROL.find("# rampRelax #") !=NIS) rampRelax = param;
    if (ROL.find("# dzRamp #") !=NIS) dzRamp = param;
    if (ROL.find("# dpRamp #") !=NIS) dpRamp = param;
    if (ROL.find("# rVeloTurnStep #") !=NIS) rVeloTurnStep = param;
    if (ROL.find("# nVeloTurnStep #") !=NIS) nVeloTurnStep = param;
    if (ROL.find("# nVeloTransition #") !=NIS) nVeloTransition = param;
    if (ROL.find("# fzOpposite #") !=NIS) fzOpposite = param;
    if (ROL.find("# zOpposite #") !=NIS)  zOpposite = param;
    if (ROL.find("# fLateral #") !=NIS) fLateral = param; 
    if (ROL.find("# vX #") !=NIS) vX = param;
    if (ROL.find("# vY #") !=NIS) vY = param;
    if (ROL.find("# fSteadySlide #") !=NIS) fSteadySlide = param; 
    if (ROL.find("# f3dMovie #") !=NIS) f3dMovie = param;
    if (ROL.find("# resolMovie #") !=NIS) resolMovie = param;
  } 

  //ENH-Maxwell
  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    if (stiffHigh[iLayer] == 0) stiffHigh[iLayer] = stiffness[iLayer];
  }

}

//AW modified:
void gfmdSheet::initMaxwell(){
  ifstream input("maxwell.in");
  int nMaxwell_test;
  double dTime_test=0, kLow_test=0, kHigh_test=0, damping_test=0, mass_test=0;

  if(!input.is_open()) termination("maxwell.in does not exist.");
  else {
    cerr << "Reading maxwell.in.\n";

    size_t pos; // current position in file stream
    std::string ROL; //REST OF LINE
    std::size_t NIS = std::string::npos; //NIS = NOT IN STRING
    
    while (input.peek()!=EOF){

      // skip empty lines and lines starting with '#'
      char firstChar;
      pos = input.tellg();
      input >> firstChar;
      input.seekg(pos,input.beg);
      if (firstChar=='#' || firstChar=='\n') {
        getline(input, ROL);
        //cout << "skipping: _" << ROL << "_\n";//DEBUG
        continue;
      }

      // read kLow, kHigh, nMaxwell
      pos = input.tellg();
      double param;
      input >> param;
      getline(input, ROL);
      if (ROL.find("# kLow #") !=NIS) kLow_test = param;
      if (ROL.find("# kHigh #") !=NIS) kHigh_test = param;
      if (ROL.find("# nMaxwell #") !=NIS) nMaxwell_test = param;
      if (ROL.find("# dTime #") !=NIS) dTime_test = param;
      if (ROL.find("# massGFMD #") !=NIS) mass_test = param;
      if (ROL.find("# dampGlobal #") !=NIS) damping_test = param;
      if (ROL.find("# MWelement #") !=NIS) {
        double iMw=-1, kMw=-1, etaMw=-1;

        input.seekg(pos,input.beg);
        input >> iMw >> kMw >> etaMw;
        getline(input,ROL);//TODO: something necessary to skip to \n?
        kMaxwell[iMw] = kMw;
        dampMaxwell[iMw] = etaMw;
        //cout << iMw << " " << kMw << " " << etaMw << endl;//DEBUG
      }
    }
  }

  //DEBUG
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    cout << "#" << iMw << " " << kMaxwell[iMw] << " " << dampMaxwell[iMw] << endl;
  }

  // sanity checks
  if (fMaxwell && !fMassWeightg) termination("fMaxwell only works with fMassWeightg.");
  if (kLow_test != stiffness[0]) termination("Lowest stiffness incompatible with maxwell.in.");
  if (kHigh_test != stiffHigh[0]) termination("Highest stiffness incompatible with maxwell.in.");
  if (nMaxwell_test != nMaxwell) termination("nMaxwell incompatible with maxwell.in.");
  if (damping_test != dampGlobal) termination("dampGlobal incompatible with maxwell.in.");
  if (mass_test != massGFMD) termination("massGFMD incompatible with maxwell.in.");
  if (dTime_test < 0.99*dTime) {
    cerr << dTime_test << "\t\t# dTime #\n";
    termination("dTime too large. Please use the above or smaller.");
  } 
  else if (dTime < dTime_test) cerr << "CAUTION! You are using a smaller dTime than suggested by maxwell.cpp!\n";

}


void gfmdSheet::writeParams(){
  if(!ID) writeParamsDefault();

  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet start\n";

  if(nElast) {
    output << nElast << "\t\t# nElast #\n";
    output << pressInit << "\t\t# pressInit #\n";
    if(pressInit!=pressFinal) output << pressFinal << "\t\t# pressFinal #\n";
    if(fConstCOM) {
      output << fConstCOM << "\t\t# fConstCOM #\n";
      output << zConstCOM << "\t\t# zConstCOM #\n";
      output << vzConstCOM << "\t\t# vzConstCOM #\n";
    }

    if(fSteppedRamp) {
      output << fSteppedRamp << "\t\t# fSteppedRamp #\n";
      output << rampSteps << "\t\t# rampSteps #\n";
      output << rampRelax << "\t\t# rampRelax #\n";
      if(fSteppedRamp==1) output << dzRamp << "\t\t# dzRamp #\n";
      else if(fSteppedRamp==2) output << dpRamp << "\t\t# dpRamp #\n";
    }

    if(rVeloTurnStep>0) output << rVeloTurnStep << "\t\t# rVeloTurnStep #\n";
    if(nVeloTurnStep>0) output << nVeloTurnStep <<"\t\t# nVeloTurnStep #\n";
    if(nVeloTransition>0) output << nVeloTransition <<"\t\t# nVeloTransition #\n";
    if(fzOpposite) {
      output << fzOpposite <<"\t\t# fzOpposite #\n";
      output << zOpposite/nxny <<"\t\t# zOpposite #\n";
      if( (!fThickness[0])&&nElast ) 
      output << thickness[0] << "\t\t# thickness0 #\n";
    }

    if(fOnSitePotential) {
      output << fOnSitePotential <<"\t\t# fOnSitePotential #\n";
      output << fOnSitePeriod <<"\t\t# fOnSitePeriod #\n";
    }

    if(fKelvinVoigt) {
      output << fKelvinVoigt <<"\t\t# fKelvinVoigt #\n";
      output << tauKV <<"\t\t# tauKV #\n";
      output << scalKV <<"\t\t# scalKV #\n";
      output << rKV_LPF <<"\t\t# rKV_LPF #\n";
    }

   //AW modified:
    if (fMaxwell) {
      output << fMaxwell << "\t\t# fMaxwell #\n";
      output << nMaxwell << "\t\t# nMaxwell #\n";
      output << fMwRelax << "\t\t# fMwRelax #\n";
      output << nMuTiSt << "\t\t# nMuTiSt #\n";
      output << fReadMwElement << "\t\t# fReadMwElement #\n";
      //output << dTimeMw << "\t\t# dTimeMw #\n";//ENH-Maxwell
    }

    for (int iLayer=0; iLayer<nLayer; ++iLayer) {
      if (stiffness[iLayer] == 0)  continue;
      output << stiffness[iLayer] << "\t\t# stiffness" << iLayer << " #\n";
      if (stiffness[iLayer] != stiffHigh[iLayer])
        output << stiffHigh[iLayer] << "\t\t# stiffHigh" << iLayer << " #\n";//ENH-Maxwell
      if (fThickness[iLayer]) {
        output << fThickness[iLayer] << "\t\t# fThickness" << iLayer << " #\n";
        output << poisson[iLayer]    << "\t\t# poisson"    << iLayer << " #\n";
        output << thickness[iLayer]  << "\t\t# thickness"  << iLayer << " #\n";
      }
      if (elastExpnt[iLayer] != 1)
        output << elastExpnt[iLayer] << "\t\t# elastExpnt" << iLayer << " #\n";
    }
    if (fMassWeightg) {
      output << fMassWeightg << "\t\t# fMassWeightg" << " #\n";
      output << zeroModeMass << "\t\t# zeroModeMass" << " #\n";
    }
    if (massGFMD != 1) output << massGFMD << "\t\t# massGFMD #\n";
      
  } // end if(nElast)

  if(fRoughRead) output << fRoughRead << "\t\t# fRoughRead #\n";

  if(fRoughAdd) {
    output << fRoughAdd << "\t\t# fRoughAdd #\n";
    if (fRoughAdd&1) {
      output << rXhertz << "\t\t# rXhertz #\n";
      if(rYhertz!=rXhertz) output << rYhertz << "\t\t# rYhertz #\n";
      if(hertzExpnt!=2) output << hertzExpnt << "\t\t# hertzExpnt #\n";
    }
    if (fRoughAdd&2) {
      output << hurst << "\t\t# hurst #\n";
      output << lambdaR << "\t\t# lambdaR #\n";
      output << lambdaS << "\t\t# lambdaS #\n";
      if(fRoughNorm!=1) output << fRoughNorm << "\t\t# fRoughNorm #\n";
      if(rRoughNorm!=1) output << rRoughNorm << "\t\t# rRoughNorm #\n";
      if(peklenik!=1)output << peklenik << "\t\t# peklenik #\n";
      output << fRollOff << "\t\t# fRollOff #\n";
      if(fBoxMuller) output << fBoxMuller << "\t\t# fBoxMuller #\n";
    }
    if (fRoughAdd&4) {
      output << fAddSWR << "\t\t# fAddSWR #\n";
      output << nqxAddSWR << "\t\t# nqxAddSWR #\n";
      if(nqxAddSWR!=nqyAddSWR) output << nqyAddSWR << "\t\t# nqyAddSWR #\n";
      output << heightSWR << "\t\t# heightSWR #\n";
    }
    if (fRoughAdd&8) {
      output << radiusFlatPunch << "\t\t# radiusFlatPunch #\n";
      output << heightFlatPunch << "\t\t# heightFlatPunch #\n";
    }
    if (fRoughAdd&16) {
      output << rSphere << "\t\t# rSphere #\n";
    }

    if(dzStep>0) output << dzStep << "\t\t# dzStep #\n";
  }

  if (resolMovie != 512) output << resolMovie << "\t\t# resolMovie #\n";
  if (f3dMovie) output << f3dMovie << "\t\t# f3dMovie #\n";

  if(fLateral) {
    output << fLateral << "\t\t# fLateral #\n";
    if(vX!=0) output << vX << "\t\t# vX #\n";
    if(vY!=0) output << vY << "\t\t# vY #\n";
  }

  if(fSteadySlide) output << fSteadySlide << "\t\t# fSteadySlide #\n";


  output << ID << "\t# sheet end" << endl << endl;
  output.close();

} // writeParams

void gfmdSheet::writeParamsDefault(){
  ofstream output("params.def",ofstream::app);
  output << "! ! ! ! ! ! ! ! ! sheet (default) values\n\n";

  output << 0    << "\t\t# nElast #\n";
  output << 0.5  << "\t\t# stiffnessN #\n";
  output << 1    << "\t\t# elastExpntN #\n\n";
  output << 0    << "\t\t# fThicknessN #\n";
  output << 0.25 << "\t\t# poissonN #\n";
  output << 1    << "\t\t# thicknessN #\n\n";

  output << 1    << "\t\t# fKelvinVoigt #\n";
  output << 1    << "\t\t# tauKV #\n";
  output << 1000 << "\t\t# scalKV #\n";
  output << 1    << "\t\t# rKV_LPF #\n\n";
  
  //AW modified:
  output << 1    << "\t\t# fMaxwell #\n";
  output << 2    << "\t\t# nMaxwell #\n";
  output << 1    << "\t\t# fMwRelax #\n";
  output << 1    << "\t\t# nMuTiSt #\n";
  output << 1    << "\t\t# fReadMwElement #\n";
  //output << 0.025 << "\t\t# dTimeMw #\n\n";//ENH-Maxwell

  output << 1 << "\t\t# fMassWeightg #\n";
  output << 1 << "\t\t# zeroModeMass #\n\n";

  output << 0.1 << "\t\t# pressInit #\n";
  output << 0.2 << "\t\t# pressFinal #\n\n";

  output << 1 << "\t\t# fOnSitePotential #\n";
  output << 1 << "\t\t# fOnSitePeriod #\n\n";

  output << 1     << "\t\t# fConstCOM #\n";
  output << 1     << "\t\t# zConstCOM #\n";
  output << 0.001 << "\t\t# vzConstCOM #\n\n";
  output << 1     << "\t\t# fSteppedRamp #\n";
  output << 0.01  << "\t\t# dzRamp #\n";
  output << 0.001 << "\t\t# dpRamp #\n";
  output << 130   << "\t\t# rampSteps #\n";
  output << 170   << "\t\t# rampRelax #\n\n";
  output << 0.5   << "\t\t# rVeloTurnStep #\n";
  output << 1000  << "\t\t# nVeloTurnStep #\n";
  output << 50    << "\t\t# nVeloTransition #\n";
  output << 1	    << "\t\t# fzOpposite #\n";
  output << 0	    << "\t\t# zOpposite #\n\n";

  output << 1 << "\t\t# fLateral #\n";
  output << 1 << "\t\t# fSteadySlide #\n";
  output << 1 << "\t\t# vX #\n";
  output << 1 << "\t\t# vY #\n\n";

  output << 1 << "\t\t# fRoughAdd #\n";
  output << 1 << "\t\t# fRoughRead #\n\n";
  
  output << 1 << "\t\t# rXhertz #\n";
  output << 1 << "\t\t# rYhertz #\n";
  output << 2 << "\t\t# hertzExpnt #\n\n";

  output << 0.8 << "\t\t# hurst #\n";
  output << 0.5 << "\t\t# lambdaR #\n";
  output << .05 << "\t\t# lambdaS #\n";
  output << 1.0 << "\t\t# rRoughNorm #\n";
  output << 1.0 << "\t\t# peklenik #\n";
  output << 1   << "\t\t# fRollOff #\n";
  output << 1   << "\t\t# fRoughNorm #\n";
  output << 1   << "\t\t# fBoxMuller #\n\n";

  output << 1  << "\t\t# fAddSWR #\n";
  output << 2  << "\t\t# nqxAddSWR #\n";
  output << 2  << "\t\t# nqyAddSWR #\n";
  output << .1 << "\t\t# heightSWR #\n\n";

  output << 0.2 << "\t\t# radiusFlatPunch #\n";
  output << 1.  << "\t\t# heightFlatPunch #\n";

  output << 0.2 << "\t\t# rSphere #\n\n";

  output << .01 << "\t\t# dzStep #\n\n";
  
  output << 0   << "\t\t# f3dMovie #\n";
  output << 512 << "\t\t# resolMovie #\n\n";

  output.close();
}

void gfmdSheet::initSystem(){

  // allocate all fields needed by FFTW, first "commodity fields";

  // fieldRFFT not always needed, but for ease of coding always allocated
  fieldRFFT = (double *)  fftw_malloc( nReal*sizeof(double) );
  fieldFFFT = (Complex *) fftw_malloc( nFour*sizeof(Complex));

  if (fRough) initConfig();
  if (nElast) dispZ_prepare(); 
  if (fLateral) initLateral();

  kFastest = getLinCIndex(nx/2, ny/2);
  rCentral = getLinRIndex(nx/2, ny/2);

  if (fRough) dumpConfig();
}

void gfmdSheet::initConfig(){
  //cout << iTime << "\tsheet" << ID << ".initConfig()\n";//DEBUG-FLOW

  equilPos = (double *) fftw_malloc( nReal * sizeof(double) );
  for (int k = 0; k < nReal; ++k) equilPos[k] = 0;

  if (fRoughRead&1) addConfigR();
  if (fRoughRead&2) addConfigF();

  if (fRoughAdd&1) addHertz();
  if (fRoughAdd&2) addSelfAffine();
  if (fRoughAdd&4) addSWR();
  if (fRoughAdd&8) addFlatPunch();
  if (fRoughAdd&16) addSphere();

  if (dzStep > 0) makeSteps();

  if (ID) {// (rigid indenter might come from below if ID>0)
    for (int k = 0; k < nReal; ++k) equilPos[k] *= -1;
  }
}

void gfmdSheet::addConfigR(){

  string fileName = konfigName + "real";
  ifstream test(fileName);
  if (!test.is_open()) termination("File "+fileName+ " does not exist.");
  test.close();
  readReal(fieldRFFT, nx, ny, fileName, 3);
  for(Lint k = 0; k < nReal; ++k) equilPos[k] += fieldRFFT[k];

  // make surface analysis
  vector<double> props(6);
  realSpaceAnalysis(fieldRFFT, &props, nx, ny);
  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet [in] roughness info start\n";
  output << sqrt(props[0]) << "\t" << props[5] << "\t# rms and max height\n";
  output << sqrt(props[1]) << "\t" << sqrt(props[2]) << "\t# x and y rms grad\n";
  output << sqrt(props[3]) << "\t" << sqrt(props[4]) << "\t# x and y rms curv\n";
  output << sqrt(props[1]+props[2]) << "\t"  << 1. / sqrt( (props[3]+props[4])/2 )
         << "\t# rms grad and radOfCurv\n";
  output << ID << "\t# sheet [in] roughness info end\n\n";
  output.close();
  
}

void gfmdSheet::addConfigF(){
  return;

  // MM2Change: infrastructure for reading Fourier
  string fileName = konfigName + "four";
  ifstream test(fileName);
  if (!test.is_open()) termination("File "+fileName+ " does not exist.");
  test.close();
  // readFour(fieldFFFT, nx, ny, fileName, 3);
  for(Lint k = 0; k < nReal; ++k) equilPos[k] += fieldRFFT[k];

}

void gfmdSheet::addHertz() {
  for (Lint k=0; k<nxny; ++k) {
    get_irxiry(k);
    double deltaX = (irx-nxH)*dx;
    double deltaY = (iry-ny/2)*dy;
    double height = deltaX*deltaX/rXhertz + deltaY*deltaY/rYhertz;
    height = pow(height,hertzExpnt/2);
    height /= hertzExpnt;
    equilPos[k] += height;
  }
}

void gfmdSheet::addSphere() {
  for (Lint k=0; k<nxny; ++k) {
    get_irxiry(k);
    double deltaX2 = (irx-nxH)*dx; deltaX2 *= deltaX2;
    double deltaY2 = (iry-ny/2)*dy; deltaY2 *= deltaY2;
    double rSphere2 = rSphere*rSphere;
    if (deltaX2 + deltaY2 < rSphere2)
    equilPos[k] += rSphere - sqrt(rSphere2 - deltaX2 - deltaY2);
    else equilPos[k] += rSphere;
  }
}

void gfmdSheet::addSelfAffine() {

  if (lambdaS>lengthX) cerr << "# lambdaS>lengthX\n";
  if (lambdaS>lengthY) cerr << "# lambdaS>lengthY\n";

  // NX, NY of a system (just) large enough to accommodate complete spectrum
  Lint NX = 0.1+2*lengthX/lambdaS;
  Lint NY = 0.1+2*lengthY/lambdaS;

  // redefine hurst for one-dimensional interface
  if(nx==1) {NX = 1; hurst -= 0.5;}
  if(ny==1) {NY = 1; hurst -= 0.5;}

  Lint NYHP1 = NY/2 + 1;
  Lint nFOUR = NX * NYHP1;

  // fftw_plan for heightF (fieldFFFT) --> heightR (fieldRFFT)
  fftw_plan heightF2R =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, fieldRFFT, FFTW_ESTIMATE);
  for(Lint k = 0; k<nFour; ++k) fieldFFFT[k] = 0;

  double qRoll = TWOPI / lambdaR;
  double qMax  = TWOPI / lambdaS;

  double varHeight = 0, gradX2 = 0, gradY2 = 0, curvX2 = 0, curvY2 = 0;

  for (int iWarm = 0; iWarm<100; ++iWarm) rrand();

  for (Lint K = 1; K<nFOUR; ++K) {

    // compute (hypothetical) indices
    Lint IQX = K/NYHP1, IQY = K%NYHP1;
    Lint JQX = abs(NX-IQX);
    IQX = (IQX<JQX) ? IQX : JQX;

    // check if wavenumber within cutoff
    double q = sqrt( peklenik*pow(IQX,2)*dqx2 + pow(IQY,2)*dqy2/peklenik );
    if (q>qMax) continue;

    // draw random numbers first
    double phase = rrand()*TWOPI;
    double bm = 1;
    if (fBoxMuller) bm = sqrt(-2.*log(rrand()))*cos(2*PI*rrand());
    
    // compute contribution to spectrum 
    double absHeight = sqrtSpec(q);
    double specLocal = absHeight*absHeight;
    double weight=2;
    if ( (IQY==0) || (IQY==NY/2) ) weight = 1;
    specLocal *= weight;

    // compute height characteristics

    // height variance
    varHeight += specLocal;

    // slope variance
    double qxN = IQX*dqx2*IQX;
    double qyN = IQY*dqx2*IQY;
    gradX2 += qxN * specLocal;
    gradY2 += qyN * specLocal;

    // curvature variance
    qxN *= qxN; qyN *= qyN;
    curvX2 += qxN * specLocal;
    curvY2 += qyN * specLocal;

    // assign to FT(height) if it fits
    if( (IQX>nx/2) || (IQY>ny/2) ) continue;

    if(IQX!=JQX) iqx = IQX; 
    else iqx = nx-IQX;
    iqy = IQY;

    Lint k = iqx*nyHP1 + iqy;

    fieldFFFT[k] = bm * absHeight * exp(Complex(0,1)*phase);

    // enforce fieldFFFT({\bf q}) = fieldFFFT*(-{\bf q})
    if ( (iqy==0) || (iqy==ny/2) ) {
      if ( (iqx==0) || (iqx==nxH) ) {   // Fourier transform is purely real
        phase = (phase>PI) ? -1 : 1;
        fieldFFFT[k] = absHeight * phase;
      } else if (iqx>nxH) {             // complex conjugate already drawn
        Lint kConj = getLinCIndex(nx-iqx,iqy);
        fieldFFFT[k]  = conj(fieldFFFT[kConj]);
      }
    }

  } // loop over k

  // normalize surface in desired way (height or rms gradient)
  double scaleWeight = rRoughNorm*rRoughNorm, grad2 = gradX2 + gradY2;
  if(fRoughNorm==1) scaleWeight /= grad2;
  else if (fRoughNorm==2) scaleWeight /= varHeight;
  else cerr << "# no rescaling in addSelfAffine of sheet " << ID << endl;

  varHeight *= scaleWeight;
  grad2 *= scaleWeight; gradX2 *= scaleWeight; gradY2 *= scaleWeight;
  curvX2 *= scaleWeight; curvY2 *= scaleWeight;

  scaleWeight = sqrt(scaleWeight);
  for (int k=0; k<nFour; ++k) fieldFFFT[k] *= scaleWeight;
  
  // transform to real space
  fftw_execute(heightF2R);

  // shift bottom layer to zero
  double minHeight = fieldRFFT[0], maxHeight = fieldRFFT[0];
  double height1 = 0, height2 = 0;
  for (Lint k = 1; k < nxny; ++k) {
    minHeight = min(minHeight, fieldRFFT[k]);
    maxHeight = max(maxHeight, fieldRFFT[k]);
    height1 += fieldRFFT[k];
    height2 += fieldRFFT[k]*fieldRFFT[k];
  } maxHeight -= minHeight;
  height1 /= nxny;
  height2 /= nxny; height2 -= height1*height1;

  for (Lint k = 0; k < nxny; ++k) {
    double dummy = fieldRFFT[k] - minHeight;
    equilPos[k] += dummy;
  }

  // output measurements
  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet [self-affine] roughness info start\n";
  output << sqrt(height2)   << "\t\t# rms real-space height\n";
  output << sqrt(varHeight) << " " << maxHeight << "\t# rms and max height (Fourier)\n";
  output << sqrt(gradX2) << "  " << sqrt(gradY2) << "\t# x and y rms slope\n";
  output << sqrt(curvX2) << "  " << sqrt(curvY2) << "\t\t# x and y rms curve\n";
  output << sqrt(grad2) << "\t"  << 1. / sqrt( (curvX2+curvY2)/2 )
         << "\t# rms grad and radOfCurv\n";
  output << ID << "\t# sheet [self-affine] roughness info end\n\n";
  output.close();

  fftw_destroy_plan(heightF2R);

} // addSelfAffine

void gfmdSheet::addSWR() {

// fAddSWR = ...
// 1 single wave parallel to y
// 2 single wave parallel to x
// 3 square roughness, see Eq. (1) in Dapp and Muser, EPL 109, 44001 (2015)
// 4 hexagonal, see Eq. (2) in same paper
// 5 triangular, see Eq. (3)
// 6 hexagonal  again but x <--> y
// 7 triangular again but x <--> y
// 4 ... 7, correct lengthX/lengthY ratio must be imposed to yield correct results

   const double qx = nqxAddSWR*2*PI/lengthX, qy = nqyAddSWR*2*PI/lengthY;

   if( ((fAddSWR==4)||(fAddSWR==5)) && (nqyAddSWR%2) ) {
     cerr << "### commensurability issue in sheet " << ID << ".with nqy";
//   termination("");
   }

   if( ((fAddSWR==6)||(fAddSWR==7)) && (nqxAddSWR%2) ) {
     cerr << "### commensurability issue in sheet " << ID << ".with nqx";
//   termination("");
   }

   for (Lint ix=0; ix<nx; ++ix) { double x = ix*lengthX/nx;
   for (Lint iy=0; iy<ny; ++iy) { double y = iy*lengthY/ny;
     double height;
     if(fAddSWR==1) height =  1 + cos(qx*x) ;
     if(fAddSWR==2) height =  1 + cos(qy*y) ;
     if(fAddSWR==3) height =  2 + cos(qx*x) + cos(qy*y);
     if(fAddSWR==4) height = sqrt(2./3) * ( 1.5 + 2*cos(qx*x)*cos(qy*y/2) + cos(qy*y));
     if(fAddSWR==5) height = sqrt(13.5) 
			   - sqrt(2./3) * ( 1.5 + 2*cos(qx*x)*cos(qy*y/2) + cos(qy*y));
     if(fAddSWR==6) height = sqrt(2./3) * ( 1.5 + 2*cos(qy*y)*cos(qx*x/2) + cos(qx*x));
     if(fAddSWR==7) height = sqrt(13.5)
                           - sqrt(2./3) * ( 1.5 + 2*cos(qy*y)*cos(qx*x/2) + cos(qx*x));
     equilPos[ix*ny+iy] += heightSWR * height;
   } }

  if(fAddSWR>3) { // report relative deviations of wavelengths in x and y
    ofstream output("params.out", ofstream::app);
    output << ID << "\t# sheet [addSWR] roughness info start\n";
    double lambdaRatio = (lengthX/nqxAddSWR) / (lengthY/nqyAddSWR);
    if( (fAddSWR==4)||(fAddSWR==5) ) lambdaRatio /= (2./sqrt(3.));
    if( (fAddSWR==6)||(fAddSWR==7) ) lambdaRatio *= (2./sqrt(3.));
    double violation = 1 - lambdaRatio;
    output << violation << "\t\t# wavelength mismatch\n";
    if(abs(violation)>0.1) 
    cerr << "\n# " << violation << " wavelength mismatch in ID "
	 << ID << endl;
    output << ID << "\t# sheet [addSWR] roughness info end\n\n";
    output.close();
  }

} // addSWR

void gfmdSheet::addFlatPunch(){

  int fFlatPunch = 1; // MM2change (allow for read in and dump)
  double r2cut = (1.+1e-6)*radiusFlatPunch*radiusFlatPunch;
  for (Lint k=0; k<nReal; ++k) {
    get_irxiry(k);
    double deltaX = (irx-nxH)*dx;
    double deltaY = (iry-ny/2)*dy;
    double r2 = deltaX*deltaX+deltaY*deltaY;
//  if (fIndMode==0) r2 += deltaX*deltaX;
    if (r2>2*r2cut) equilPos[k] += heightFlatPunch;
    else if(r2>r2cut) {
      if (fFlatPunch==1) equilPos[k] += heightFlatPunch;
      else equilPos[k] += 
      heightFlatPunch*(cos(PI*(r2-2*r2cut)/(r2cut))+1)/2;
    }
  }


}


void gfmdSheet::makeSteps(){
  for (Lint k = 0; k < nReal; ++k) {
    double zLocal = equilPos[k];
    zLocal /= dzStep;
    zLocal = dzStep*round(zLocal);
    equilPos[k] = zLocal;
  }
}


double gfmdSheet::sqrtSpec(double q){
  double qRoll = TWOPI / lambdaR;
  double qRed = q/qRoll;
  if (fRollOff==1) return( pow(1./sqrt(1.+qRed*qRed),1+hurst) ); // smooth roll-off
  if (qRed<1.-1e-9) {
    if (fRollOff==0) return(0.);        // cut-off
    else if (fRollOff==2) return(1.);   // roll-off with cusp
  }
  return(pow(1./qRed,1+hurst));
}

void gfmdSheet::dispZ_prepare(){

  // initialize fields
  dispZR = (double *)  fftw_malloc(nReal*sizeof(double) );
  dispZF = (Complex *) fftw_malloc(nFour*sizeof(Complex));
  int iDummy = 0;
  if(fConstCOM) iDummy = 1;
  memset(dispZR, iDummy*zConstCOM, sizeof(double)*nReal);

  stressZR = (double *)  fftw_malloc( nReal*sizeof(double) );
  memset(stressZR, 0, sizeof(double)*nReal);

  dispZF    = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  dispZFold = (Complex *) fftw_malloc( nFour*sizeof(Complex) );

  for (Lint k = 1; k < nFour; ++k) dispZFold[k] = dispZF[k] = 0.;

  stressZF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
  if(fKelvinVoigt) rhsKV_LPF = (Complex *) fftw_malloc( nFour*sizeof(Complex) );


  weightQ = (short int *) fftw_malloc( nFour*sizeof(short int) );

  if(!fKelvinVoigt) massCorr = (double *) fftw_malloc( nFour*sizeof(double) );
  if(fSteadySlide) preFacCorrSS = (Complex *) fftw_malloc( nFour*sizeof(Complex) );

  stressZFPre = (double *) fftw_malloc( nFour*sizeof(double) );
  
  //AW modified: Maxwell field initialization
  if (fMaxwell){
    dispZMw.resize(nMaxwell);
    stressZMw.resize(nMaxwell);
    for (int iMaxwell = 0; iMaxwell < nMaxwell; ++iMaxwell){
      dispZMw[iMaxwell] = (Complex *) fftw_malloc( nFour*sizeof(Complex) );
      stressZMw[iMaxwell] = (Complex *) fftw_malloc( nFour*sizeof(Complex));
    }
    for (Lint k = 0; k < nFour; ++k){
      for (int iMaxwell = 0; iMaxwell < nMaxwell; ++iMaxwell){
        dispZMw[iMaxwell][k] = 0.;
        stressZMw[iMaxwell][k] = 0.;
      }
    }
  } //prepare for the Maxwell displacement and stress field

  // read displacemenmt

  string fileName = konfigName + "old";
  ifstream test(fileName);
  if(test.is_open()) {
    test.close();
    if ( (nx!=1)&&(ny!=1) ) {
      if (fRough) readReal(dispZR, nx, ny, fileName, 4);
      else readReal(dispZR, nx, ny, fileName, 3);
    } else {
      if (fRough) readReal(dispZR, nx, ny, fileName, 3);
      else readReal(dispZR, nx, ny, fileName, 2);
    }
  }

  // make all plans

  dispZR2F =
  fftw_plan_dft_r2c_2d(nx, ny, dispZR, (fftw_complex*) dispZF, FFTW_ESTIMATE);

  dispZF2R =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, dispZR, FFTW_ESTIMATE);

  stressZR2F =
  fftw_plan_dft_r2c_2d(nx, ny, stressZR, (fftw_complex*) stressZF, FFTW_ESTIMATE);

  // initialize (Fourier) displacements from file

  fftw_execute(dispZR2F);
  for (int k = 0; k < nFour; ++k) dispZFold[k] = dispZF[k];
  if(fConstCOM) {
    dispZFold[0] = dispZF[0] = zConstCOM*nxny;
    dispF2R();
  }

  if(fSteadySlide){
    fConstCOM = 1;
    zConstCOM = dispZF[0].real()/nxny;
  }

  // initialize stress-constant prefactors, weights, and Q-dependent inertia

  weightQ[0] = 1;
  stressZFPre[0] = 0;
  // fThickness==2 is the only case, for which stressZFPre[0]!=0.
  for (int iLayer=0; iLayer<nLayer; ++iLayer) {
    if (fThickness[iLayer] == 2) { //MM2XYZ double check if correct 
      double nu = poisson[iLayer], thick = thickness[iLayer], stiff = stiffness[iLayer];
      stressZFPre[0] += (stiff/thick) * (2*(1.-nu)*(1.-nu))/(1.-2*nu); // Carbone
      // stressZFPre[0] += (stiff/thick) * (1.-nu*nu); // E_Young
    }
  }

  // change MM2SS: Better way to initialize?
  if(fKelvinVoigt) {
    for (Lint k = 0; k < nFour; ++k) rhsKV_LPF[k] = 0;
  }

  for (Lint k = 1; k < nFour; ++k) {
    double q = getQ(k);
    for (int iLayer=0; iLayer<nLayer; ++iLayer) {
      double preFac = 1;
      double width = thickness[iLayer]*q;
      double width2 = width*width;
      // from Carbone et al. Eur. Phys. J. E 29, 275-284 (2009)
      const double POISSON = poisson[iLayer];
      if (fThickness[iLayer]==1) {              // bottom boundary: constant stress
        preFac = cosh(2*width)-2*width2-1;
        preFac /= sinh(2*width)+2*width;
        if ( width>40 ) preFac = 1;
      } else if (fThickness[iLayer]==2) {       //bottom boundary: constant strain
        preFac  = (3-4*POISSON)*cosh(2*width)+2*width2-4*POISSON*(3-2*POISSON)+5;
        preFac /= (3-4*POISSON)*sinh(2*width)-2*width;
        if (width>40) preFac = 1;
      }

      // Green's function with width-dependet q-Factor
      stressZFPre[k] += preFac * pow(q,elastExpnt[iLayer]) * stiffness[iLayer];
    } // iLayer end

    weightQ[k] = 2;
    if ( (iqy==0) || (iqy==ny/2) ) weightQ[k] = 1;

  } // loop over q-vectors */

  if(!fKelvinVoigt) {
    for (Lint k=0; k<nFour; ++k) massCorr[k] = 1;
  }

  if(fMassWeightg) {
    if (fMaxwell || zeroModeMass==0) {
      double stiffRatio = stiffHigh[0]/stiffness[0];
      for (int k = 1; k < nFour; ++k) {
        massCorr[k] = massGFMD*stressZFPre[k]/stiffness[0] / stiffMax;
        massCorr[k] = 1./massCorr[k] ;
      }
      massCorr[0] = 2*massCorr[1];
    }
    else {
      for (int k = 0; k < nFour; ++k) {
        massCorr[k] = pow(zeroModeMass*stiffMin,2.) + pow(stressZFPre[k],2);
        massCorr[k] = sqrt(massCorr[k]) / stiffMax;
        massCorr[k] = 1./massCorr[k] ;
      }
    }
    damping = dampGlobal; 
  } 
  kFastest = getLinCIndex(nx/2, ny/2);

  if(fSteadySlide) {
    for (int k = 0; k < nFour; ++k) {
      // compute frequency of mode
      get_iqxiqy(k);
      if(iqx>nxH)   iqx -= nxH;
      if(iqy>nyHP1) iqy -= nyHP1;
      double qx = iqx*dqx, qy = iqy*dqy;
      double omega = qx*vX + qy*vY;
      if(fKelvinVoigt) {
        preFacCorrSS[k]  = (1.+Complex(0,1)*omega*tauKV/scalKV);
        preFacCorrSS[k] /= (1.+Complex(0,1)*omega*tauKV);
      } else {
        termination("# Steady-state sliding only implemented for KV.");
      }
    }
  }

  ofstream output("params.out", ofstream::app);
  output << ID << "\t# sheet elasticity info start\n";
  if(!fConstCOM) output << pressure*areaXY << "\t\t# total force\n";
  output << stiffMin << "\t\t# stiffMin\n";
  output << stiffMax << "\t\t# stiffMax\n";
  output << stressZFPre[1] << "\t\t# stressZFPre[1]\n";
  output << stressZFPre[kFastest] << "\t\t# stressZFPre[kFastest]\n";
  if (!fKelvinVoigt){
    output << massScal/massCorr[0] << "\t\t# massCOM\n";
    output << massScal/massCorr[kFastest] << "\t\t# massFast\n";
  }
  output << damping  << "\t\t# damping\n";
  if(stressZFPre[0]!=0) output << stressZFPre[0] << "\t\t# stressZFPre[0]\n";
  output << ID << "\t# sheet elasticity info end\n\n";
  output.close();


  //DEBUG
  if (!fKelvinVoigt) {
    cerr << "\n4kLow/m = " << 4*stressZFPre[1] * massCorr[1]/stiffMax;
    cerr << "\n4kHigh/m = " << 4*stressZFPre[1]*stiffHigh[0]/stiffness[0] * massCorr[1]/stiffMax;
    cerr << "\ng^2 = " << damping*damping << "\n";
  }

} // dispZ_prepare

void gfmdSheet::initLateral(){

  if ( (vX!=0) || (vY!=0) ) {
    fSliding = true;		// reset class variable
    fExternalDriving = true;	// reset global variable
  }

  if (fRough) {
    equilPos0F = (Complex *) fftw_malloc( nFour * sizeof(Complex) );
    for(int k=0; k<nReal; ++k) fieldRFFT[k] = equilPos[k]; 
    fftw_plan equilPosR2F =
    fftw_plan_dft_r2c_2d(nx, ny, fieldRFFT, (fftw_complex*) equilPos0F, FFTW_ESTIMATE);
    fftw_execute(equilPosR2F);
    fftw_destroy_plan(equilPosR2F);
    equilPosF2R =
    fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) fieldFFFT, equilPos , FFTW_ESTIMATE);
  }
}

void gfmdSheet::moveLateral() {

  double xShiftAbs = vX*mdTime, xShiftRel = vX*dTime;
  double yShiftAbs = vY*mdTime, yShiftRel = vY*dTime;

  if (fRough) {
    for (Lint k = 0; k < nFour; ++k) {
      get_iqxiqy(k);
      if (iqx >= nx/2) iqx-=nx;
      double qx = iqx*dqx, qy = iqy*dqy;
      fieldFFFT[k] = exp(Complex(0,1)*(qx*xShiftAbs+qy*yShiftAbs)) * equilPos0F[k];
      fieldFFFT[k] /= (double) nxny;
    }
    fftw_execute(equilPosF2R);
  }

  if (nElast) {
    termination("dispZFold might have to be shifted too!"); 
    for (Lint k = 0; k < nFour; ++k) {
      get_iqxiqy(k);
      if (iqx >= nx/2) iqx-=nx;
      double qx = iqx*dqx, qy = iqy*dqy;
      fieldFFFT[k] = exp(Complex(0,1)*(qx*xShiftRel+qy*yShiftRel)) * dispZF[k];
      fieldFFFT[k] /= (double) nxny;
    }
    fftw_execute(dispZF2R);
  }

  return;

}

double gfmdSheet::getQ(Lint k){
  get_iqxiqy(k);
  int jqx = abs(iqx-nx);
  jqx = (iqx<jqx) ? iqx : jqx;
  double q2 = jqx*jqx*dqx2 + iqy*iqy*dqy2;
  return (sqrt(q2));
}

void gfmdSheet::dumpConfig(){
  //cout << iTime << "\tsheet" << ID << ".dumpConfig()\n";//DEBUG-FLOW

  string fileName = konfigName + "dat";

  if ( (!fRoughAdd) && (!nElast) ) { // dump cross section of input surface
    vector<double*> arrays = {equilPos};
    vector<string> arrayNames = {"Z_eq"};
    dumpRealCS(arrays, arrayNames, nx, ny, 1, 1, fileName+"H");
    return;
  }

  // fields to be dumped
  vector<double*> arrays;
  vector<string> arrayNames;
  if (fSheetType&1) {
    arrays.push_back(equilPos);
    arrayNames.push_back("Z_eq");
  }
  if (fSheetType&2) {
    arrays.push_back(dispZR);
    arrayNames.push_back("dispZ");
    arrays.push_back(stressZR);
    arrayNames.push_back("stressZ");
  }

  if(fSheetType==1) {
    dumpReal(arrays, nx, ny, 1, 1, fileName);
    dumpRealCS(arrays, arrayNames, nx, ny, 1, 1, fileName+"H");
    return;
  }

  // compute elastic stress, start ...

  for(Lint k = 0; k < nFour; ++k) stressZF[k] = 0.;
  for(Lint k = 0; k < nReal; ++k) stressZR[k] = 0.;

  if(fKelvinVoigt) stressKV();
  else stressGFMD();

  if(!fConstCOM)
  stressZF[0] += pressure * nxny;  // external stress as part of elastic stress

  // Fourier based filtering of stress field would go here

  // transform elastic force to real space
  fftw_plan stressZF2R  =
  fftw_plan_dft_c2r_2d(nx, ny, (fftw_complex*) stressZF, stressZR, FFTW_ESTIMATE);
  fftw_execute(stressZF2R);
  for (Lint k = 0; k < nxny; ++k) stressZR[k] /= nxny;

  // ... end, compute elastic stress

  dumpReal(arrays, nx, ny, 1, 1, fileName);
  dumpRealCS(arrays, arrayNames, nx, ny, 1, 1, fileName+"H");

}

void gfmdSheet::dumpFrame(int iFrame){

  // no need to dump frame of static equilPos
  if ((fSheetType<2)&&(!fLateral)) return; 

  int xStep = (nx>resolMovie) ? (nx-1)/resolMovie + 1: 1;
  int yStep = (ny>resolMovie) ? (ny-1)/resolMovie + 1: 1;

  string fileName = "frames/" + konfigName + to_string(iFrame) + ".dat";

  for (int k = 0; k < nReal; ++k) {
    fieldRFFT[k] = 0;
    if (nElast) fieldRFFT[k]  = dispZR[k];
    if (fRough) fieldRFFT[k] += equilPos[k];
  }

  // fields to be dumped
  vector<double*> arrays = {fieldRFFT};
  vector<string> arrayNames = {"Z"};
  if (fSheetType>1) {
    arrays.push_back(stressZR);
    arrayNames.push_back("stressZ");
  }

  if (f3dMovie) dumpReal(arrays, nx, ny, xStep, yStep, fileName);
  dumpRealCS(arrays, arrayNames, nx, ny, xStep, yStep, fileName+"H");

}

void gfmdSheet::zeroStress(){
  // memset better? 
  if(!nElast) return;
  for (int k = 0; k < nFour; ++k) stressZF[k] = 0;
  if(fRealSpaceStress) for (int k = 0; k < nReal; ++k) stressZR[k] = 0;
}

double gfmdSheet::addStress(){
  //if (iTime <= 2) cout << iTime << "\tsheet" << ID << ".addStress()\n";//DEBUG-FLOW

  if (!nElast) return(0.);

  vOnSitePot = 0.;

  // complete real-space stress with on-site potential
  if(fOnSitePotential&1) vOnSitePot += deltaIndenter();
  if(fOnSitePotential&2) vOnSitePot += einsteinSolid();

  // initialize Fourier stress with FFTW from real-space stress
  memset(stressZF, 0, sizeof(Complex) * nFour);
  if (fRealSpaceStress) fftw_execute(stressZR2F);
  stressInter = stressZF[0].real()/nxny;

  // add elastic and thermal stresses
  vElastic = stressGFMD(); 
  if ( fLangevin && (temp!=0) ) thermostat();
  
  //AW modified:
  if (fMaxwell) addStressMaxwell();
  //cout << mdTime << "\t" << stressZF[1].real() << "\t" << abs(stressZF[kFastest]) << "\n";//TEMP-stressRelax

  // add external pressure
  if(fConstCOM) vGravit = 0;
  else {
    stressZF[0] += pressure * nxny;
    vGravit = - pressure * dispZF[0].real() * areaPP;
  }
  stressCOM = stressZF[0].real()/nxny;

  if(fSteadySlide) chi2SteadySlide = computeChi2SS();
  if(iTime==0){
    chi2SteadySlideOld = chi2SteadySlide;
    chi2SSRestart = 1e+50*chi2SteadySlide;
    if (chi2SteadySlide<1e-10) chi2SSRestart = 1e+50*lengthX*lengthY;
    weightSS = 1e-5;
    invStiffSS_COMRef = 1./stressZFPre[1];
    invStiffSS_COM = invStiffSS_COMRef;
    oldSSPress = 1;
    prefacSSnxny = pow(nxny,0.9);
  }

  vPotTot = vOnSitePot+vElastic+vGravit;
  return(vPotTot);
}

//ENH-Maxwell
void gfmdSheet::addStressMaxwell(){
  //if (iTime <= 2) cout << iTime << "\tsheet" << ID << ".addStressMaxwell()\n";//DEBUG-FLOW

  // compute new Maxwell stress if necessary
  int kStart = 0;
  if (fThickness[0]==0) kStart = 1;//TODO: other cases necessary?
  for (Lint k = kStart; k < nFour; ++k){
    for (int iMaxwell = 0; iMaxwell < nMaxwell; ++iMaxwell){
      //TODO: nMuTiSt!
      double stiffZZ = stressZFPre[k]*kMaxwell[iMaxwell]/stiffness[0];
      double dampZZ = dampMaxwell[iMaxwell]*stressZFPre[k]/stiffness[0];
      stressZMw[iMaxwell][k] = stiffZZ*(dispZF[k] - dispZMw[iMaxwell][k]);
      dispZMw[iMaxwell][k] += dTime*stressZMw[iMaxwell][k] / dampZZ;//ENH-Maxwell
    }
  }
  /*//DEBUG-stressRelax start
  double stressReal = 0, stressImag = 0;
  for (int iMw = 0; iMw < nMaxwell; ++iMw) {
    stressReal += stressZMw[iMw][1].real();
    stressImag += stressZMw[iMw][1].imag();
  }
  cout << mdTime << "\t" << stressZF[1].real() << "\t" << stressReal << "\n";
  //DEBUG-stressRelax end*/

  // add Maxwell stress
  for (Lint k = 0; k < nFour; ++k) {
    for (int iMw=0; iMw<nMaxwell; ++iMw) stressZF[k] -= stressZMw[iMw][k];
  }//TEMP-check*/
}

void gfmdSheet::fireHalt(){
  if (!nElast) return;
//memcpy(dispZFcopy, dispZFold, nFour*sizeof(Complex));
  memcpy(dispZF, dispZFold,  nFour*sizeof(Complex));
  dispF2R();
}

void gfmdSheet::fireRedirect(){

  if ( (!nElast) || (fireRedrct==0) ) return;

  Complex vDotV = 0., vDotF = 0., fDotF = 0;
  for(Lint k = 0; k < nFour; ++k){
    Complex velocity = dispZF[k] - dispZFold[k];
    vDotV += weightQ[k] * 1. * velocity   * conj(velocity);
    vDotF += weightQ[k] * 1. * velocity   * conj(stressZF[k]);
    fDotF += weightQ[k] * 1. * stressZF[k] * conj(stressZF[k]);
  }

  double delta_PH = (1. - fireRedrct) * vDotF.real() / fDotF.real();
  double delta_Q  = fireRedrct*(fireRedrct - 2) * vDotV.real() / fDotF.real();

  double scaleV = (1. - fireRedrct);
  double scaleF = -delta_PH + sqrt(delta_PH*delta_PH - delta_Q);

  for(Lint k = 0; k < nFour; ++ k){
    dispZF[k]  = dispZFold[k] + scaleV * (dispZF[k] - dispZFold[k]);
    dispZF[k] += scaleF * stressZF[k];
  }

  if(fRealSpaceStress) dispF2R();
}

double gfmdSheet::deltaIndenter(){
  const int iPoint = 0*rCentral;
  double stress = -pressure * nxny;
  if(fOnSitePeriod!=0) stress *= sin(fOnSiteFreq*mdTime);
  stressZR[iPoint] += stress;
  return (-stress*dispZR[iPoint]*areaPP);
}

double gfmdSheet::einsteinSolid(){
  double vPot = 0;
  double zEq = 1, kEinstein = sqrt(stiffMax*stiffMin);
  if (iTime==0) cerr << "# kEinstein = " << kEinstein
		     << " (stress)\t" << kEinstein*areaPP << " (per CGA)" << endl;
  for (Lint k = 0; k < nxny; ++k) {
    double dz = dispZR[k]-zEq;
    stressZR[k] -= kEinstein*dz;
    vPot += dz*dz;
  } 
  vPot *= kEinstein*areaPP/2;
  return(vPot); 
}

double gfmdSheet::stressGFMD(){
  //if (iTime <= 2) cout << iTime << "\tsheet" << ID << ".stressGFMD()\n";//DEBUG-FLOW
  double vElastic = 0;

  // treat 0 mode separately 
  if (fzOpposite) {
    Complex stressLoc = stressZFPre[0] * (dispZF[0] - zOpposite);
    stressZF[0] -= stressLoc;
    Complex dummy = stressLoc * conj(dispZF[0] - zOpposite);
    vElastic += weightQ[0]*dummy.real();
  }

  if(fKelvinVoigt) return(vElastic);

  for (Lint k = 1; k<nFour; ++k) {
    Complex stressLoc = stressZFPre[k] * dispZF[k];
    stressZF[k] -= stressLoc;
    Complex dummy = stressLoc * conj(dispZF[k]);
    vElastic += weightQ[k]*dummy.real();
  }

  vElastic *= areaPP / 2;
  vElastic /= nxny;

  return(vElastic);
}

double gfmdSheet::stressKV(){
  double vElastic = 0;

  for (Lint k = 1; k != nFour; ++k) {
    stressZF[k] = -stressZFPre[k]*dispZF[k];
    Complex dummy = -stressZF[k]*conj(dispZF[k]);
    vElastic += weightQ[k]*dummy.real();
  }
  vElastic *= areaPP/2;
  vElastic /= nxny;
  return(vElastic);
    
}

void gfmdSheet::thermostat() {

  double preFac = sqrt(1.5*nxny*damping*massScal*temp*0.5/areaPP) / dTime; 
  for (Lint k = 0; k != nFour; ++k) {
    double randZR = 2*mix64() - 1;
    double randZI = 2*mix64() - 1;
    Complex randZC = Complex(randZR, randZI);
    stressZF[k] += preFac*randZC / sqrt(2*massCorr[k]);
    get_iqxiqy(k);
  }

}

void gfmdSheet::initMeasure(){
  //cout << iTime << "\tsheet" << ID << ".initMeasure()\n";//DEBUG-FLOW
  if(!nElast) return;

  tKinetic1 = tKinetic2 = 0;

  string moniName = "moni" + to_string(ID) + "-";
  string rampName = "ramp" + to_string(ID) + "-";
  if (ny<1000) {moniName = moniName + "0"; rampName = rampName + "0";}
  if (ny<100)  {moniName = moniName + "0"; rampName = rampName + "0";}
  if (ny<10)   {moniName = moniName + "0"; rampName = rampName + "0";}
  moniName = moniName + to_string(ny) + ".dat";
  rampName = rampName + to_string(ny) + ".dat";
  moni.open(moniName, ofstream::out);
  if (fSteppedRamp) {
    ramp.open(rampName, ofstream::out);
    ramp << "# displacement  stress  time\n";
  }
  moni << "# mdTime"
  << "\tdispZF[0].real()"
  << "  dispZF[1].real()"
  << "  dispZF[kFastest].real()"
  << "  dispZR[0]"
  << "  dizpZR[central]";
  if (fConstCOM) moni << "  stressCOM"; 
  if ( (pressFinal!=pressInit) || (fSteppedRamp==2) ) moni << "  pressure";
  if (fzOpposite) moni << "  zOpposite";
  moni << endl;
}

void gfmdSheet::measure(){
  //if (iTime <= 2) cout << iTime << "\tsheet" << ID << ".measure()\n";//DEBUG-FLOW
  if(!nElast) return;

  tKinetic1 += tKinetic;
  tKinetic2 += tKinetic*tKinetic;

  if (iTime) {
    moni << mdTime << "\t";
    moni << dispZF[0].real()/nxny
         << "\t" << dispZF[1].real()/nxny
         << "\t" << dispZF[kFastest].real()/nxny
         << "\t" << dispZR[0]
         << "\t" << dispZR[rCentral];
    if (fConstCOM) {
      if (fzOpposite) moni << "\t" << stressInter;
      else moni << "\t" << stressCOM;
    }
    if ( (pressFinal!=pressInit) || (fSteppedRamp==2) ) moni << "\t" << pressure;
    if (fzOpposite) moni << "\t" << zOpposite/nxny;
    moni << endl;
  }

}

double gfmdSheet::tKinHalfStep(){

  if(!nElast) return(0);

  // initialize prefactors
  double dt2 = dTime2 / stiffMax;
  double weightNow = 2, weightOld = 1;

  tKinetic = 0;
  int kStart = 0;
  if (fConstCOM) kStart = 1;
  for (Lint k = kStart; k < nFour; ++k) {
    Complex vDiffN = dispZF[k] - dispZFold[k];
    vDiffN = vDiffN * conj(vDiffN);
    tKinetic += weightQ[k] * vDiffN.real() / massCorr[k];
  }
  tKinetic *= areaPP / (2*dt2*nxny);
  return(tKinetic);
}

double gfmdSheet::propagate(){
  //if (iTime <= 2) cout << iTime << "\tsheet" << ID << ".propagate()\n";//DEBUG-FLOW

  // MM2MM: need to add moveLateral on displacements (for nElast==1 sheets)
  if(fSliding&&!fSteadySlide) moveLateral(); 

  if(!nElast) return(0.);

  //switch off FIRE while being ramped
  if(fSteppedRamp) {
    rampTime = iTime % (rampSteps+rampRelax);
    if(fFire && (rampTime<rampSteps)) fFireOn = 0;
  }

  // update pressure
  if (fSteppedRamp != 2) {
    pressure = pressInit;
    if(iTime>nRelax) pressure += (iTime-nRelax)*(pressFinal-pressInit)/(nTime-nRelax);
  } else if (iTime>nRelax) {
    if ( (rampTime<rampSteps) && (iTime>rampSteps) ) {
      double prefacRamp = sin(rampTime*PI/rampSteps);
      prefacRamp *= 2*prefacRamp;
      pressure += ddpRamp*prefacRamp;
    } else if (rampTime+1==rampSteps+rampRelax) { // dump info at last relax step
      ramp << dispZF[0].real()/nxny << "\t";
      ramp << stressInter << "\t" << mdTime << endl;
    }
  }

  // viscoelasticity is treated by propagateKV instead
  if(fKelvinVoigt) {
    if(fSteadySlide) {
      propagateSteadySlide();
      return(0); 
    } 
    propagateKV();
    return(0); 
  }

  // initialize prefactors
  double dt2 = dTime2 / stiffMax;
  double weightNow = 2, weightOld = 1;

  if(!fFireOn) {
    double dampDL = damping*dTime; // damping in units of dTime
    dt2 /= (1.0 + dampDL);
    double weightDiff = 1-exp(-dampDL);
    weightNow -= weightDiff;
    weightOld -= weightDiff;
  }
  
  // external propagation of COM/Opposite
  int kStart = 0;
  if (fConstCOM) {
    if ( (fzOpposite) && (iTime>nRelax) ) moveOpposite();
    else {
      kStart = 1;
      if(iTime>nRelax) moveCOM();
    }
  }

  // propagate internal q modes with Verlet
  tKinetic = 0;
  for (Lint k = kStart; k < nFour; ++k) {
    Complex dispFnew = weightNow*dispZF[k] - weightOld*dispZFold[k]
		     + massCorr[k]*stressZF[k]*dt2;
    Complex vDiffN = dispFnew - dispZFold[k];
    vDiffN = vDiffN * conj(vDiffN);
    tKinetic += weightQ[k] * vDiffN.real() / massCorr[k]; 
    dispZFold[k] = dispZF[k];
    dispZF[k] = dispFnew;//TEMP-stressRelax
  }
  tKinetic *= areaPP / (8*dt2*nxny);

  if(fFire&&(fireRedrct!=0)) fireRedirect();

  // turn around movement if necessary
  if ((nVeloTransition>0) && (iTime>=tStartTransition) && (iTime<tEndTransition)) { 
    if (iTime == t0Transition) t0Transition -= 1; // skip time step where vz would be 0
    vzConstCOM = vzInit*(t0Transition - iTime)/nVeloTransition;
  }
  else if (iTime==nVeloTurnStep) {
    vzConstCOM *= -1;
    ddzRamp *= -1;
    ddpRamp *= -1;
  }

  dispF2R();
  return(tKinetic);
}

void gfmdSheet::moveCOM(){
  //cout << iTime << "\tsheet" << ID << ".moveCOM()\n";//DEBUG-FLOW

  dispZF[0] += vzConstCOM*dTime*nxny;

  if (fSteppedRamp) {
    if ( (rampTime<rampSteps) && (iTime>rampSteps) ) {
      double prefacRamp = sin(rampTime*PI/rampSteps);
      prefacRamp *= 2*prefacRamp;
      dispZF[0] += ddzRamp*prefacRamp;
    } else if (rampTime+1==rampSteps+rampRelax) {
      ramp << dispZF[0].real()/nxny << "\t";
      ramp << stressInter << "\t" << mdTime << endl;
    }
    dispZFold[0] = dispZF[0];
  }

}

void gfmdSheet::moveOpposite(){
  zOpposite += vzConstCOM*dTime*nxny;

  if (fSteppedRamp) {
    if ( (rampTime<rampSteps) && (iTime>rampSteps) ) {
      double prefacRamp = sin(rampTime*PI/rampSteps);
      prefacRamp *= 2*prefacRamp;
      zOpposite += ddzRamp*prefacRamp;
    } else if (rampTime+1==rampSteps+rampRelax) {
      ramp << zOpposite/nxny << "\t";
      ramp << stressInter << "\t" << mdTime << endl;
    }
  }
}


double gfmdSheet::propagateKV(){

  // initialize prefactors
  double weightNew = 1./(1+rKV_LPF);
  double weightOld = 1-weightNew;
  double prefacDot = weightNew*tauKV/scalKV;

  // external propagation of COM/Opposite
  int kStart = 0;
  if (fConstCOM) {
    if ( (fzOpposite) && (iTime>nRelax) ) moveOpposite();
    else {
      kStart = 1;
      if(iTime>nRelax) moveCOM();
    }
  }

  // propagate internal q modes with KV-SLS
  double dissipKV = 0.; // MM2fix: dissipKV not yet computed
  for (Lint k = kStart; k < nFour; ++k) {

    // NOTE: dispZFold is actually stressZFold here.
    Complex extStressDot = (stressZF[k]-dispZFold[k])/dTime; 
    Complex dispLocOld = dispZF[k];

    rhsKV_LPF[k] *= weightOld;
    rhsKV_LPF[k] += weightNew*stressZF[k];
    if(iTime>1) rhsKV_LPF[k] += prefacDot*extStressDot;

    //0 mode fix
    Complex effectiveStress;
    if (k==0) effectiveStress = rhsKV_LPF[k];
    else effectiveStress = rhsKV_LPF[k] - stressZFPre[k] * dispZF[k];


    Complex  velocity = effectiveStress/ ( stressZFPre[k]*tauKV);
    if(stressZFPre[k]==0) velocity = effectiveStress/ ( stressZFPre[1]*tauKV); //0 mode fix
    dispZF[k] += velocity * dTime;

    dispZFold[k] = stressZF[k];

  }

  // turn around movement if necessary
  if ((nVeloTransition>0) && (iTime>=tStartTransition) && (iTime<tEndTransition)) {
    if (iTime == t0Transition) t0Transition -= 1; // skip time step where vz would be 0
    vzConstCOM = vzInit*(t0Transition - iTime)/nVeloTransition;
  }
  else if (iTime==nVeloTurnStep) {
    vzConstCOM *= -1;
    ddzRamp *= -1;
    ddpRamp *= -1;
  }

  dispF2R();
  
  return(dissipKV);

}

double gfmdSheet::computeChi2SS(){
  Complex chi2 = 0;
  for (int k = 1; k < nFour; ++k) {
    stressZF[k] *= preFacCorrSS[k]/stressZFPre[k];
    Complex diff = stressZF[k] - dispZF[k];
    chi2 += diff*conj(diff);
  }
  return(chi2.real()); 
}

void gfmdSheet::propagateSteadySlide(){

  if (chi2SteadySlide >= chi2SteadySlideOld){
    chi2SteadySlide = chi2SteadySlideOld;
    weightSS /= 4;
    if(weightSS<1e-10) {
      stressCOM = stressZF[0].real()/nxny;
      double pressDiff = pressInit + stressCOM;
      if (abs(pressDiff) < pressInit*1.e-5){
        dispF2R();
        cerr << "# relative pressure error<1e-5 induces finish";
        ofstream finish("finish");
        finish.close();
      } else {
        // adjusting the effective mass
        if (oldSSPress*pressDiff > 0) invStiffSS_COM*=pow(2.,1./4);
        else {
	  if (invStiffSS_COM > invStiffSS_COMRef) invStiffSS_COM = invStiffSS_COMRef;
	  invStiffSS_COM /= 2;
          oldSSPress *= -1;
        }
        dispZF[0] += invStiffSS_COM * pressDiff * prefacSSnxny;
        
        chi2SteadySlideOld = chi2SSRestart;
        weightSS = 1.e-5; // restart the weight
        for (int k = 1; k < nFour; ++k) dispZFold[k] = dispZF[k];
      }
    }
  } else {
    chi2SteadySlideOld = chi2SteadySlide;
    for (int k = 1; k < nFour; ++k){
      Complex newDispZF;
      newDispZF = weightSS*stressZF[k] + (1.-weightSS)*dispZFold[k];
      dispZFold[k] = dispZF[k];
      dispZF[k] = newDispZF;
    }
    weightSS *= pow(2.,1./4);
    if(weightSS>4) weightSS = 4;
  }
 
  dispF2R();
}

void gfmdSheet::dispF2R(){
  // transform from Fourier to real space with back-up copy
  memcpy(fieldFFFT, dispZF, nFour*sizeof(Complex));
  fftw_execute(dispZF2R);
  for (Lint k = 0; k < nxny; ++k) dispZR[k] /= nxny;
}

void gfmdSheet::dispR2F() {
  fftw_execute(dispZR2F);
  // MM2fix: if we ever start with finite velocities, then we need
  //		a real copy of dispZF at iTime==0. 
  if (iTime==0) memcpy(dispZFold, dispZF, nFour*sizeof(Complex));
}

void gfmdSheet::outMeasure(){

  if(!nElast) return;

  ofstream output("params.out",ofstream::app);
  output << ID << "\t\t# sheet measure start\n";

  if(vElastic!=0) output << vElastic << "\t\t# vElastic\n"; 
  if(vOnSitePot!=0) output << vOnSitePot << "\t\t# vOnSitePot\n"; 
  if(vGravit!=0) output << vGravit << "\t\t# vGravit\n"; 
  if (areaXY!=1) {
    if(vElastic!=0) output << vElastic/areaXY << "\t\t# vElastic/areaXY\n"; 
    if(vOnSitePot!=0) output << vOnSitePot/areaXY << "\t\t# vOnSitePot/areaXY\n"; 
    if(vGravit!=0) output << vGravit/areaXY << "\t\t# vGravit/areaXY\n"; 
  }

  tKinetic1 /= (iTime-nRelax);
  tKinetic2 /= (iTime-nRelax);
  tKinetic2 -= tKinetic1*tKinetic1;
  if(fLangevin || fExternalDriving) {
    output << tKinetic1 << "\t" << tKinetic2 << "\t#  tKinetic1, tKinetic2\n";
  }

  output << dispZF[0].real()/nxny << "\t\t# zConstCOM #\n";
  if (fzOpposite) output << zOpposite/nxny << "\t\t# zOpposite #\n";

  output << ID << "\t\t# sheet measure end\n\n";
  output.close();
}

void gfmdSheet::outSystem(){
  if(nElast||fLateral) dumpConfig();
}

gfmdSheet::~gfmdSheet(){
  //DEBUG
  if (!fKelvinVoigt) {
    ofstream gfmd("GFMD.out");
    gfmd << "# q\tk(q)/m(q)\n";
    for (Lint k = 0; k<nFour; ++k) {
      gfmd << getQ(k) << "\t" << stressZFPre[k]*massCorr[k] << "\n";
    }
    gfmd.close();
  }

  if(nElast) {
    moni.close();
    if (fSteppedRamp) ramp.close();
  }

  fftw_free(fieldRFFT);
  fftw_free(fieldFFFT);
  fftw_free(equilPos);
  fftw_free(dispZR);
  fftw_free(dispZF);
  fftw_free(dispZFold);
  fftw_free(stressZR);
  fftw_free(stressZF);
  fftw_free(weightQ);
  fftw_free(massCorr);
  fftw_free(stressZFPre);
  fftw_free(equilPos0F);

}

