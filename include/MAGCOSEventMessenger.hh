#ifndef MAGCOSEventMessenger_h
#define MAGCOSEventMessenger_h 1
// DESCRIPTION
// -----------
//
// This class defines the event messenger for MAGNETOCOSMICS. 
// It defines UI commond that allows to control interactively 
// the MAGCOSEventAction object. 
//  
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
// 

#include "globals.hh"
#include "G4UImessenger.hh"

class MAGCOSEventAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class MAGCOSEventMessenger: public G4UImessenger
{
  public:
    MAGCOSEventMessenger(MAGCOSEventAction * myAct);
    ~MAGCOSEventMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    MAGCOSEventAction* myAction;
    G4UIdirectory*     myDrawDir;
   
 
    
    G4UIcmdWith3Vector* SetDrawColourCmd;
    G4UIcmdWithABool*   SetDrawTrajectoryCmd;    
    G4UIcmdWithAString* SetDrawingCoordinateSystemCmd;
    G4UIcmdWithABool*  SetDrawPointsCmd;
    G4UIcmdWithADouble* SetPointSizeCmd;
    G4UIcmdWithoutParameter* DrawCmd;
    G4UIcmdWithoutParameter* ResetCmd;
    G4UIcmdWithADoubleAndUnit* TraceMagnetopauseLineCmd; 
    
    //commands for registring trajectories and field lines in a file
    G4UIcmdWithABool* SetRegisterTrajectoryCmd;
    G4UIcmdWithAString* SaveTrajectoriesCmd;
   
    
 };

#endif

