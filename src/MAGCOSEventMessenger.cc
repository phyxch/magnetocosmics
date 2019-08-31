#include "G4RunManager.hh"
#include "MAGCOSEventMessenger.hh"
#include "MAGCOSEventAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "MAGCOSSpenvisManager.hh"

////////////////////////////////////////////////////////////////////////////////
//
MAGCOSEventMessenger::MAGCOSEventMessenger(MAGCOSEventAction * myAct)
:myAction(myAct)
{ 
  G4String cmd_name;
  G4String guidance;
  myDrawDir = new G4UIdirectory("/MAGCOS/DRAW/");
  myDrawDir->SetGuidance("Define parameters for drawing trajectories");
  
 //Unfortunately line width and line style are not yet consider 
 // in G4SceneHandler when  drawing G4Polyline 
 
 /* SetDrawLineWidthCmd = new
             G4UIcmdWithADouble("/MAGCOS/DRAW/SetLineWidth",this);
  SetDrawLineWidthCmd->SetGuidance("Set the  line width for drawing trajectories  and magnetic field lines");
  SetDrawLineWidthCmd->SetParameterName("LineWidth",false);
  SetDrawLineWidthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawLineStyleCmd = new
             G4UIcmdWithAnInteger("/MAGCOS/DRAW/SetLineStyle",this);
  SetDrawLineStyleCmd->SetGuidance("Set the  line Style for drawing trajectories  and magnetic field lines");
  SetDrawLineStyleCmd->SetParameterName("LineStyle",false);
  SetDrawLineStyleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);*/
  
  
  SetDrawColourCmd = new
             G4UIcmdWith3Vector("/MAGCOS/DRAW/SetColour",this);
  SetDrawColourCmd->SetGuidance("Set the colour for drawing the next computed trajectories  and magnetic field lines");
  SetDrawColourCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  SetDrawingCoordinateSystemCmd = new
             G4UIcmdWithAString("/MAGCOS/DRAW/SetCoordinateSystem",this);
  SetDrawingCoordinateSystemCmd->SetGuidance("Set the reference coordinate system for drawing trajectories and field lines");
  SetDrawingCoordinateSystemCmd->SetParameterName("LineWidth",false);
  SetDrawingCoordinateSystemCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  
  SetDrawTrajectoryCmd = new
             G4UIcmdWithABool("/MAGCOS/DRAW/DrawTrajectory",this);
   	     
  
  guidance ="If true the next computed trajectories and/or field lines are stocked in a vector ";
  guidance +="of curves that can be plotted at any time";
  SetDrawTrajectoryCmd->SetGuidance(guidance);
  SetDrawTrajectoryCmd->SetParameterName("DrawTrajectory",false);
  SetDrawTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetDrawPointsCmd = new
             G4UIcmdWithABool("/MAGCOS/DRAW/DrawPoints",this);
  guidance ="If true the step points of the next computed ";
  guidance +="trajectories and/or field lines  are stoked in a vector of points ";
  guidance +="that can be plotted at any time"; 	     
  SetDrawPointsCmd->SetGuidance(guidance);
  SetDrawPointsCmd->SetParameterName("DrawPoints",false);
  SetDrawPointsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SetPointSizeCmd = new
             G4UIcmdWithADouble("/MAGCOS/DRAW/SetPointSize",this);
  SetPointSizeCmd->SetGuidance("Set the size of step points for plotting purpose");
  SetPointSizeCmd->SetParameterName("StepSize",false);
  SetPointSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DrawCmd = new
             G4UIcmdWithoutParameter("/MAGCOS/DRAW/Show",this);
 
  guidance ="Show/Plot the curves and step positions that have been ";
  guidance +="computed and registered during previous computations";
  DrawCmd->SetGuidance(guidance);
  DrawCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ResetCmd = new
             G4UIcmdWithoutParameter("/MAGCOS/DRAW/Reset",this);
  ResetCmd->SetGuidance("Reset the vector of stored particle trajectories and magnetic field  lines ");
  ResetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  //commands for registring trajectories and field lines in a file
  
   
  SetRegisterTrajectoryCmd = new
             G4UIcmdWithABool("/MAGCOS/DRAW/RegisterTrajectory",this);
  guidance ="If true the trajectory or field lines that should be drawn ";
  guidance +="are also registered as  a SpenvisCVS object for later saving in a file" ; 	     
  SetRegisterTrajectoryCmd->SetGuidance(guidance);
  SetRegisterTrajectoryCmd->SetParameterName("RegisterTrajectory",false);
  SetRegisterTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  SaveTrajectoriesCmd = new
             G4UIcmdWithAString("/MAGCOS/DRAW/SaveTrajectories",this);
  SaveTrajectoriesCmd->SetGuidance("The trajectories that have been regsiterd in SencvisCVS objects are saved in a file ");
  SaveTrajectoriesCmd->SetParameterName("LineWidth",false);
  SaveTrajectoriesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
 /* cmd_name = "/MAGCOS/DRAW/MagnetopauseLine";
  TraceMagnetopauseLineCmd = new G4UIcmdWithADoubleAndUnit(cmd_name,this);
  TraceMagnetopauseLineCmd->SetGuidance("Trace a line on the magnetopause");
  TraceMagnetopauseLineCmd->SetParameterName("theta",false);
  TraceMagnetopauseLineCmd->SetDefaultValue(0.);
  TraceMagnetopauseLineCmd->SetUnitCategory("Angle");
  TraceMagnetopauseLineCmd->SetDefaultUnit("degree");*/
  
}
////////////////////////////////////////////////////////////////////////////////
//

MAGCOSEventMessenger::~MAGCOSEventMessenger()
{ delete   myDrawDir;
  delete   SetDrawTrajectoryCmd;
  delete   SetDrawColourCmd;
  delete   SetDrawingCoordinateSystemCmd;
  delete   SetDrawPointsCmd;
  delete   SetPointSizeCmd;
  delete   ResetCmd;
  delete   DrawCmd;
  delete   TraceMagnetopauseLineCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSEventMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if( command == SetDrawTrajectoryCmd ) 
            myAction->SetDrawTrajectory(SetDrawTrajectoryCmd
                              ->GetNewBoolValue(newValues));
   
  else if( command == SetDrawingCoordinateSystemCmd ) 
            myAction->SetDrawingCoordinateSystem(newValues);
   
   
  else if( command == SetDrawColourCmd ){ 
  	G4ThreeVector vec=SetDrawColourCmd->GetNew3VectorValue(newValues);
	myAction->SetDrawColour(G4Colour(vec.x(),vec.y(),vec.z()));
  }			        	      
   
  else if( command == SetDrawPointsCmd ) 
            myAction->SetDrawPoints(SetDrawPointsCmd
                              ->GetNewBoolValue(newValues));
   
  else if( command == SetPointSizeCmd ) 
            myAction->SetPointSize(SetPointSizeCmd
                              ->GetNewDoubleValue(newValues));			      
  else if( command == DrawCmd ) 
                 myAction->DrawTrajectoriesAndFieldLines(.5,45.,45.);
  else if( command == ResetCmd ) {
  	myAction->ResetVectorObjectToBeDrawn();
   	MAGCOSSpenvisManager::getInstance()
			->ResetTheParticleTrajectoryCSVBlocks();
  }
  else if( command == SetRegisterTrajectoryCmd ) 
            MAGCOSSpenvisManager::getInstance()->SetRegisterParticleTrajectory(SetRegisterTrajectoryCmd
                              ->GetNewBoolValue(newValues));
  else if( command == SaveTrajectoriesCmd ) 
            MAGCOSSpenvisManager::getInstance()->SaveTrajectoryCSVBlocks(newValues);
   
   			      	
  /*else if (command == TraceMagnetopauseLineCmd){
  	myAction->TraceMagnetopauseLine(TraceMagnetopauseLineCmd
	                                      ->GetNewDoubleValue(newValues)); 
  }*/
 		       
}

