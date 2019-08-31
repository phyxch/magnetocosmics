#include "G4RunManager.hh"
#include "MAGCOSGeometryMessenger.hh"
#include "MAGCOSGeometryConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
// #include "G4VViewer.hh"

MAGCOSGeometryMessenger::MAGCOSGeometryMessenger(MAGCOSGeometryConstruction * myDet)
:myGeometry(myDet)
{ 

  myUserLimitDir = new G4UIdirectory("/MAGCOS/USERLIMIT/");
  myUserLimitDir->SetGuidance("Define user limit");
  
 
  MaxStepLengthCmd = new
             G4UIcmdWithADoubleAndUnit("/MAGCOS/USERLIMIT/SetMaxStepLength",this);
  MaxStepLengthCmd->SetGuidance("Set the maximum tracking step length");
  MaxStepLengthCmd->SetParameterName("MaxStepLength",false);
  MaxStepLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  MaxTrackLengthCmd = new
             G4UIcmdWithADoubleAndUnit("/MAGCOS/USERLIMIT/SetMaxTrajectoryLength",this);
  MaxTrackLengthCmd
    ->SetGuidance("Set the maximum trajectory length");
  MaxTrackLengthCmd->SetParameterName("MaxTrackLength",false);
  MaxTrackLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  MaxTrackLengthCmd->SetUnitCategory("Length");
 
  MaxTimeCmd = new
             G4UIcmdWithADoubleAndUnit("/MAGCOS/USERLIMIT/SetMaxTrajectoryTime",this);
  G4String guidance = "Set an upper limit to the time it takes ";
  guidance += "for a particle to follow the computed trajectory. \n";
  guidance +="When this time is reached the tracing of the particle ";
  guidance +="trajectory is stopped."; 
  MaxTimeCmd->SetGuidance(guidance);
  MaxTimeCmd->SetParameterName("MaxTime",false);
  MaxTimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  MaxTimeCmd->SetUnitCategory("Time");
  
  SetStopAltitudeCmd = new
             G4UIcmdWithADoubleAndUnit("/MAGCOS/USERLIMIT/SetStopAltitude",this);
  guidance = "Set stop altitude ";
  SetStopAltitudeCmd->SetGuidance(guidance);
  SetStopAltitudeCmd->SetParameterName("StopAltitude",false);
  SetStopAltitudeCmd->SetRange("StopAltitude>0.");
  SetStopAltitudeCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  SetStopAltitudeCmd->SetUnitCategory("Length");
  
 /* RemoveEarthCmd = new
             G4UIcmdWithoutParameter("/MAGCOS/USERLIMIT/RemoveEarth",this);
  guidance = "Remove the Earth";
  RemoveEarthCmd->SetGuidance(guidance);
  RemoveEarthCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 

  
  ConsiderEarthCmd = new
             G4UIcmdWithoutParameter("/MAGCOS/USERLIMIT/ConsiderEarth",this);
  guidance = "The Earth is taken into accountg";
  ConsiderEarthCmd->SetGuidance(guidance);
  ConsiderEarthCmd->AvailableForStates(G4State_PreInit,G4State_Idle); */

  
  
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSGeometryMessenger::~MAGCOSGeometryMessenger()
{ delete MaxStepLengthCmd;
  delete myUserLimitDir;
  delete MaxTrackLengthCmd;
  delete MaxTimeCmd; 
  delete SetStopAltitudeCmd;
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  
  if( command == MaxStepLengthCmd ) 
            myGeometry->SetMaxStepLength(MaxStepLengthCmd
                              ->GetNewDoubleValue(newValues));
  else if( command == MaxTrackLengthCmd ) 
            myGeometry->SetMaxTrackLength(MaxTrackLengthCmd
                              ->GetNewDoubleValue(newValues));
  else if( command == MaxTimeCmd ) 
            myGeometry->SetMaxTimeOfFlight(MaxTimeCmd
                              ->GetNewDoubleValue(newValues));				      		      
  else if( command == SetStopAltitudeCmd ) 
            myGeometry->SetStopAltitude(SetStopAltitudeCmd
                              ->GetNewDoubleValue(newValues));	
 /* else if( command == RemoveEarthCmd ) 
            myGeometry->SetRemoveEarth(true);
  else if( command == RemoveEarthCmd ) 
            myGeometry->SetRemoveEarth(false);  */                        			      
}

