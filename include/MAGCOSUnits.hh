// Updated on 9/16/2014, Hexc, Olesya: added CLHEP namespace for g4.10 version
#ifndef MAGCOSUNITS_HH
#define MAGCOSUNITS_HH 
// DESCRIPTION
// -----------
//
// This interface define new units of length and magtnetic field used in
// magnetocosmics. 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#include "globals.hh"
static const G4double re=6371.2*km;
//static const G4double re=2439.7*km;
static const G4double Re=6371.2*km;
//static const G4double Re=2439.7*km;
static const G4double nT=1.e-9*tesla;
static const G4double nanotesla=1.e-9*tesla;
static const G4double GV=GeV;
static const G4double MV=MeV;




#endif
