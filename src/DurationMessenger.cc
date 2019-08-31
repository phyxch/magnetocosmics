// 9/22/2014: Hexc & Olesya - Added CLHEP namespace
//
#include"DurationManager.hh"
#include"DurationMessenger.hh"
#include"G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//units
#include "G4UnitsTable.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

DurationMessenger::DurationMessenger (DurationManager* aDurationManager )
{ theDurationManager = aDurationManager;

  
 //directories
  DurationDir = new G4UIdirectory("/MAGCOS/TIME/");
  DurationDir->SetGuidance("Duration control");
  
 
// Duration commands

  SetMaxTotalDurationCmd = new
               G4UIcmdWithADoubleAndUnit("/MAGCOS/TIME/SetMaxTotalDuration",this);
  SetMaxTotalDurationCmd->SetGuidance("Set the maximum duration of the execution of the program");
  SetMaxTotalDurationCmd->SetParameterName("MaxTotalDuration",false);
  SetMaxTotalDurationCmd->SetRange("MaxTotalDuration >0.");
  SetMaxTotalDurationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetMaxTotalDurationCmd->SetUnitCategory("Time");
  
  
  SetMaxRunDurationCmd = new
               G4UIcmdWithADoubleAndUnit("/MAGCOS/TIME/SetMaxRunDuration",this);
  SetMaxRunDurationCmd->SetGuidance("Set the maximum duration of a run");
  SetMaxRunDurationCmd->SetParameterName("MaxRunDuration",false);
  SetMaxRunDurationCmd->SetRange("MaxRunDuration >0.");
  SetMaxRunDurationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetMaxRunDurationCmd->SetUnitCategory("Time");
  
  SetMaxEventDurationCmd = new
               G4UIcmdWithADoubleAndUnit("/MAGCOS/TIME/SetMaxEventDuration",this);
  SetMaxEventDurationCmd->SetGuidance("Set the maximum duration of an event");
  SetMaxEventDurationCmd->SetParameterName("MaxEventDuration",false);
  SetMaxEventDurationCmd->SetRange("MaxEventDuration >0.");
  SetMaxEventDurationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetMaxEventDurationCmd->SetUnitCategory("Time");
  
  
    
}
////////////////////////////////////////////////////////////////////////////////
//
DurationMessenger::~DurationMessenger()
{ delete SetMaxTotalDurationCmd;
  delete SetMaxRunDurationCmd;
  delete SetMaxEventDurationCmd;
}		  
////////////////////////////////////////////////////////////////////////////////
//
void DurationMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{  
  
  if( command == SetMaxTotalDurationCmd)
            theDurationManager->SetMaxTotalDuration(SetMaxTotalDurationCmd->GetNewDoubleValue(newValues)/s);
   
  else if ( command == SetMaxRunDurationCmd)
            theDurationManager->SetMaxRunDuration(SetMaxRunDurationCmd->GetNewDoubleValue(newValues)/s);
  
  else if ( command == SetMaxEventDurationCmd)
            theDurationManager->SetMaxEventDuration(SetMaxEventDurationCmd->GetNewDoubleValue(newValues)/s);
   
  	   
}
