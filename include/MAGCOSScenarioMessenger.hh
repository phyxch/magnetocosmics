#ifndef MAGCOSScenarioMessenger_h
#define MAGCOSScenarioMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MAGCOSApplicationScenario;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

class MAGCOSScenarioMessenger : public G4UImessenger
{public:
   MAGCOSScenarioMessenger(MAGCOSApplicationScenario* ApplicationScenario);
  ~MAGCOSScenarioMessenger();
  
   void SetNewValue(G4UIcommand * command,G4String newValues);
 
 private:
   MAGCOSApplicationScenario* theApplicationScenario;
   G4UIdirectory*             ScenarioDir;
   
   
 
   // Scenario command
   
   G4UIcmdWithoutParameter* BlineCmd;
   G4UIcmdWithoutParameter* ParticleTrajectoryCmd;
   G4UIcmdWithoutParameter* ReverseParticleTrajectoryCmd;
   G4UIcmdWithAString* ComputeRigidityFilterCmd;
   G4UIcommand* ComputeDirectionFilterCmd;
   G4UIcommand* RCutoffVsPositionCmd;
   G4UIcommand* RCutoffVsPositionOnLShellCmd;
   G4UIcommand* RCutoffVsSpenvisPositionGridCmd;
   G4UIcommand* RCutoffVsSpenvisTrajectoryCmd;
   G4UIcommand* RCutoffVsDirectionCmd;
   G4UIcommand* RCutoffVsTimeCmd;
   G4UIcmdWithABool* SetAutoDetectionOfPenumbra;
   
   //SpenvisCSV command
   G4UIcmdWithABool* SetRegisterResultsInSpenvisCSVFileCmd;
   G4UIcmdWithAString* SetSpenvisCSVFileNameCmd;
   
   
   
   
   
   
    

};












#endif
