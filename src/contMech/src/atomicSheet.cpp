#include "header.h"
#include "atomicSheet.h"
#include "gfmdSheet.h"
#include "globals.h"

// additional global variables
extern vector <gfmdSheet> sheet;

void atomicSheet::initParams(int newID){

  ID = newID;

  initParamsDefault();

  ifstream input("params.in");
  if (input.is_open()) nAtom = 0;
  int fReadParams = 0;
  while (input.is_open()) {
    double param;
    std::string ROL; // rest of line
    std::size_t NIS = std::string::npos; // NIS == Not In String
    if (input.eof()) break;
    input >> param; getline(input,ROL);
    if ( (param==ID) && (ROL.find("# atomL start") !=NIS) ) fReadParams = 1;
    if ( fReadParams==0 ) continue;
    if (ROL.find("# atomL end") !=NIS) break;
    if (ROL.find("# nAtom #")   !=NIS) nAtom = param;
    if (ROL.find("# nAtomSheet #") !=NIS) nAtomSheet = param;
    if (ROL.find("# atomSheet0 #") !=NIS) atomSheet[0] = param;
    if (ROL.find("# atomSheet1 #") !=NIS) atomSheet[1] = param;
    if (ROL.find("# mass #")   !=NIS) mass = param;
    if (ROL.find("# epsAtomWall #")!=NIS) epsAtomWall = param;
    if (ROL.find("# sigma #")   !=NIS) sigma = param;
  }
  input.close();

}

void atomicSheet::initParamsDefault(){

  nAtom = nxGlobal * nyGlobal;

  nAtomSheet = 2;
  atomSheet[0] = 0;
  atomSheet[1] = 1;

  mass = 1;
  epsAtomWall = 1./(13*12-2*6*7);

  if(nyGlobal==1) sigma = lengthX / nxGlobal;
  else sigma = lengthY / nyGlobal;

}

void atomicSheet::initSystem(){

  // allocate rNow, rOld, force seperately to generate consecutive memory
  rNow.resize(nAtom);
  for (int iAtom=0; iAtom<nAtom; ++iAtom) rNow[iAtom].resize(nDim);
  rOld.resize(nAtom);
  for (int iAtom=0; iAtom<nAtom; ++iAtom) rOld[iAtom].resize(nDim);
  force.resize(nAtom);
  for (int iAtom=0; iAtom<nAtom; ++iAtom) force[iAtom].resize(nDim);
  if(fFire&&(fireRedrct!=0)) {
    rOldCopy.resize(nAtom);
    for (int iAtom=0; iAtom<nAtom; ++iAtom) rOldCopy[iAtom].resize(nDim);
  }

  // initialize rNow
   bool fRead = readConfig();
  if(!fRead) { // set-up random structure
    for(int iAtom=0; iAtom<nAtom; ++iAtom) {
      rNow[iAtom][0] = mix64()*lengthX;
      rNow[iAtom][1] = mix64()*lengthY;
      rNow[iAtom][2] = mix64();
    }
  }

  for (int iAtom=0; iAtom<nAtom; ++iAtom)
  for (int iDim=0; iDim<nDim; ++iDim) rOld[iAtom][iDim] = rNow[iAtom][iDim];

  dumpConfig();

}

void atomicSheet::writeParams(){

  ofstream output("params.out",ofstream::app);
  output << ID << "\t# atomL start\n";

  output << nAtom << "\t\t# nAtom #\n";
  if(nAtomSheet!=2)   output << nAtomSheet   << "\t\t# nAtomSheet #\n";
  if(atomSheet[0]!=0) output << atomSheet[0] << "\t\t# atomSheet0 #\n";
  if(atomSheet[1]!=1) output << atomSheet[1] << "\t\t# atomSheet1 #\n";
  output << mass << "\t# mass #\n";
  output << epsAtomWall << "\t# epsAtomWall #\n";
  output << sigma << "\t# sigma #\n";

  output << ID << "\t# atomL end\n\n";
  output.close();

  if(!ID) writeParamsDefault();

}

void atomicSheet::writeParamsDefault(){
  ofstream output("params.def",ofstream::app);
  output << "! ! ! ! ! ! ! ! ! selected atomicSheet (default) values\n\n";
  output << nAtom << "\t\t# nAtom #\n";
  output << 2   << "\t\t# nAtomSheet #\n";
  output << 0 << "\t\t# atomSheet0 #\n";
  output << 1 << "\t\t# atomSheet1 #\n";
  output << 1. << "\t\t# mass #\n";
  output << 0.0138889 << "\t# epsAtomWall #\n";
  output << 1 << "\t\t# sigma #\n\n";
  output.close();

}

void atomicSheet::dumpConfig(){
  string fileName = "konfigA" + to_string(ID) + ".dat";
  ofstream konfig(fileName);
  konfig << nAtom << "\n\n"; 
  for (int iAtom=0; iAtom<nAtom; ++iAtom) {
    konfig << "Ar";
    for (int iDim=0; iDim<nDim; ++iDim) konfig << "\t" << rNow[iAtom][iDim];
    konfig << "\n";
  }
  konfig.close();
}

bool atomicSheet::readConfig(){
  string fileName = "konfigA" + to_string(ID) + ".old";
  ifstream konfig(fileName);
  if(!konfig.is_open()) return(false);
  Lint nAtomOld;
  konfig >> nAtomOld;
  if(nAtom!=nAtomOld) {
    cerr << "# nAtom and nAtomOld inconsistent in " << ID << endl;
    return(false);
  }
  for (int iAtom=0; iAtom<nAtom; ++iAtom) {
    konfig >> fileName;
    for (int iDim=0; iDim<nDim; ++iDim) konfig >> rNow[iAtom][iDim];
  }
  konfig.close();
  return(true);
}

void atomicSheet::initMeasure(){
  tKin1 = tKin2 = 0;
}

double atomicSheet::propagate(){
  tKinetic = 0;

  double dt2 = dTime2 / mass;

  double dampDL = dampGlobal*dTime; // damping in units of dTime
  if(fFireOn) dampDL = 0;

  dt2 /= (1.0 + dampDL);
  double weightNow = 2, weightOld = 1;
  double weightDiff = 2*dampDL / (1+dampDL);
  weightNow -= weightDiff;
  weightOld -= weightDiff;

  for (int iAtom = 0; iAtom < nAtom; ++iAtom) {
  for (int iDim = 0; iDim < nDim; ++iDim) {
    double rNew = weightNow*rNow[iAtom][iDim] - weightOld*rOld[iAtom][iDim]
		+ force[iAtom][iDim]*dt2;
    double dr = rNew - rOld[iAtom][iDim];
    tKinetic += dr*dr;
    rOld[iAtom][iDim] = rNow[iAtom][iDim];
    rNow[iAtom][iDim] = rNew;
  } }

  tKinetic *= mass/(8*dTime2);
  
  if( fFireOn && (fireRedrct!=0) ) {
    for (int iAtom = 0; iAtom < nAtom; ++iAtom) {
    for (int iDim = 0; iDim < nDim; ++iDim) {
      rOldCopy[iAtom][iDim] = rNow[iAtom][iDim];
    } }
  }

  return(tKinetic);
}
