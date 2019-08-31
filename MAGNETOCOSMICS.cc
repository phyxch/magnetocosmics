// Updated on 8/6/2014, Hexc, Olesya
// Rebuild magnetocosmics with Geant4.9.6
// Updated on 9/16/2014, Hexc, Olesya: added CLHEP namespace for g4.10 version
// Updated on 9/29/2014, Hexc, Olesya: Fixing the OPENGl display problem following
//      the dnaphysics.cc example. We simply replaced all of the display-related header
//      files based on the dnaphysics example.
// Updated on 7/22/2015, Hexc, Olesya - update the main code for using newer visualization manager
//
#include "G4RunManager.hh"

#include "MAGCOSGeometryConstruction.hh"
#include "MAGCOSPhysicsList.hh"
#include "MAGCOSApplicationScenario.hh"
#include "MAGCOSPrimaryGeneratorAction.hh"
#include "MAGCOSEventAction.hh"
#include "MAGCOSTrackingAction.hh"
#include "MAGCOSSteppingAction.hh"
#include "MAGCOSMagneticField.hh"
#include "DurationManager.hh"
#include "G4UImanager.hh"

#include "MAGCOSUnits.hh"
#include "G4UnitsTable.hh"

#include "SpaceCoordinateConvertor.hh"
#include "G4ProcessTable.hh"

#include "G4UImanager.hh"
//#include "G4UIterminal.hh"
//#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

using namespace CLHEP;

int main(int argc,char** argv) 
{
  
  // Units
  // Definition of new units in the unit table should be defined 
  // at the beginning before the
  // instantiation of the runManager and should be followed by 
  // G4UnitDefinition::BuildUnitsTable()  
  
  new G4UnitDefinition("earth radii","re","Length",re);
  new G4UnitDefinition("earth radii 1","Re","Length",re);
  new G4UnitDefinition("earth radii 2","RE","Length",re);
  new G4UnitDefinition("hour","hour","Time",3600.*s);
  new G4UnitDefinition("minute","minute","Time",60.*s);
  new G4UnitDefinition("day","day","Time",24.*3600.*s);
  new G4UnitDefinition("nanotesla","nT","Magnetic flux density",nT);
  new G4UnitDefinition("gigavolt","GV","Electric potential",1000.*megavolt);
  
  new G4UnitDefinition("/cm3","#/cm3","Density",1./cm3);
  new G4UnitDefinition("/m3","#/m3","Density",1./m3);
  new G4UnitDefinition("/km3","#/km3","Density",1./km3);
  new G4UnitDefinition("/dm3","#/dm3","Density",1.e-3/cm3);
  
  new G4UnitDefinition("km/s","km/s","Speed",km/s);
  new G4UnitDefinition("m/s","m/s","Speed",m/s);
   
  G4UnitDefinition::BuildUnitsTable();

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  MAGCOSGeometryConstruction* MAGCOSgeometry = new MAGCOSGeometryConstruction;
  runManager->SetUserInitialization(MAGCOSgeometry);
  runManager->SetUserInitialization(new MAGCOSPhysicsList);


  //Duration manager
  
  DurationManager* theDurationManager;
  theDurationManager = DurationManager::getInstance();
  
  // UserAction classes
  
  MAGCOSPrimaryGeneratorAction* 
    thePrimaryAction = new MAGCOSPrimaryGeneratorAction();	     
  runManager->SetUserAction(new MAGCOSApplicationScenario); 	     
  runManager->SetUserAction(thePrimaryAction);
  runManager->SetUserAction(new MAGCOSTrackingAction);
  runManager->SetUserAction(new MAGCOSEventAction);
  MAGCOSSteppingAction* theSteppingAction =new MAGCOSSteppingAction();
  MAGCOSgeometry->SetSteppingAction(theSteppingAction);
  runManager->SetUserAction(theSteppingAction);

  //Initialize G4 kernel  
  runManager->Initialize(); 
  
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  //  G4UIsession* session = 0;
  G4UIExecutive* session = 0;

  // Get the pointer to the User Interface manager   
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  
  if (argc==1)   // Define UI session for interactive mode.
  { 
    /*
#ifdef _WIN32
    G4UIsession * session = new G4UIterminal();
#else
    G4UIsession * session = new G4UIterminal(new G4UItcsh);
#endif
    UI->ApplyCommand("/control/execute prerun.g4mac");    
    */
    
    session = new G4UIExecutive(argc, argv);
    session->SessionStart();
    delete session;
  }
  else           // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
    
    session = new G4UIExecutive(argc, argv);
    delete session;
  }

  
  G4ProcessTable::GetProcessTable()->SetProcessActivation("MYTransportation",false);
  G4ProcessTable::GetProcessTable()->SetProcessActivation("Transportation",true);
  
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager; 
  
  return 0;
}
