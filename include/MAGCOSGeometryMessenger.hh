#ifndef MAGCOSGeometryMessenger_h
#define MAGCOSGeometryMessenger_h 1

// DESCRIPTION
// -----------
//
// This class defines the geometry messenger for MAGNETOCOSMICS. 
// It defines UI commond that allows to control interactively 
// the MAGCOSGeometryContruction object.
//  
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
// 

#include "globals.hh"
#include "G4UImessenger.hh"

class MAGCOSGeometryConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class MAGCOSGeometryMessenger: public G4UImessenger
{
  public:
    MAGCOSGeometryMessenger(MAGCOSGeometryConstruction * myDet);
    ~MAGCOSGeometryMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    MAGCOSGeometryConstruction* myGeometry;
    G4UIdirectory*             myUserLimitDir;
    
   
    G4UIcmdWithADoubleAndUnit* MaxStepLengthCmd;
    G4UIcmdWithADoubleAndUnit* MaxTrackLengthCmd;
    G4UIcmdWithADoubleAndUnit* MaxTimeCmd;    
    G4UIcmdWithADoubleAndUnit* SetStopAltitudeCmd;
    G4UIcmdWithoutParameter* RemoveEarthCmd;
    G4UIcmdWithoutParameter* ConsiderEarthCmd;
    
    
    
    
 };

#endif

