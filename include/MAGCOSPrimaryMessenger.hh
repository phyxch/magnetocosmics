#ifndef MAGCOSPrimaryMessenger_h
#define MAGCOSPrimaryMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MAGCOSPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class MAGCOSPrimaryMessenger: public G4UImessenger
{
  public:
    MAGCOSPrimaryMessenger(MAGCOSPrimaryGeneratorAction * aGenerator);
    ~MAGCOSPrimaryMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    MAGCOSPrimaryGeneratorAction* myGenerator;
    G4UIdirectory*             myStartGeoDir;
    G4UIdirectory*             myRigidityFilterDir;
   
    
    //Commands for defining primaries 
    G4UIcmdWithADoubleAndUnit* SetRigidityCmd;
    G4UIcommand* SetPositionVectorCmd;
    G4UIcommand* SetPositionCmd;
    G4UIcommand* SetPositionOnDipoleMagneticShellCmd;
    G4UIcommand* SetDirectionCmd;
    G4UIcommand* SetDirectionVectorCmd;
    G4UIcommand* SetPositionAndDirectionVectorCmd;
    G4UIcommand* SetPositionAndDirectionCmd;
    G4UIcommand* SetDirectionFromPitchAngleCmd;
    
    //Commands for defining rigidity vector
    G4UIcommand* AddValuesToRigidityVectorCmd;
    G4UIcmdWithoutParameter* SetDefaultRigidityVectorCmd;
    G4UIcmdWithoutParameter* ResetRigidityVectorCmd;

    //Verbosity
    G4UIcmdWithAnInteger*    SetVerbosityCmd;
    G4UIcmdWithoutParameter* PrintBfieldAtPrimaryCmd;
    
 };

#endif

