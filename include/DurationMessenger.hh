#ifndef DurationMessenger_h
#define DurationMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DurationManager;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

class DurationMessenger: public G4UImessenger
{public:
   DurationMessenger(DurationManager* aDurationManager);
  ~DurationMessenger();
  
   void SetNewValue(G4UIcommand * command,G4String newValues);
 
 private:
     DurationManager* theDurationManager;
     G4UIdirectory*   DurationDir;
   
  //set time duration cmd
     G4UIcmdWithADoubleAndUnit* SetMaxTotalDurationCmd; 
     G4UIcmdWithADoubleAndUnit* SetMaxEventDurationCmd;
     G4UIcmdWithADoubleAndUnit* SetMaxRunDurationCmd;     
};
#endif
