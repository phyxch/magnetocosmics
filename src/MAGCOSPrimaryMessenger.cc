#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "MAGCOSPrimaryMessenger.hh"
#include "MAGCOSPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4VSensitiveDetector.hh" 
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "SpaceCoordinateConvertor.hh"


/////////////////////////////////////////////////////////////////////////////////
//
MAGCOSPrimaryMessenger::MAGCOSPrimaryMessenger(MAGCOSPrimaryGeneratorAction
* aGenerator)
:myGenerator(aGenerator)
{ G4UIparameter* param;
  G4String cmd_name;
  G4String guidance;

  myStartGeoDir = new G4UIdirectory("/MAGCOS/SOURCE/");
  guidance ="Definition of the start position and direction ";
  guidance += "of a particle in different space coordinate systems";
  myStartGeoDir->SetGuidance(guidance);
   
  myRigidityFilterDir = new G4UIdirectory("/MAGCOS/RIGIDITYVECTOR/");
  guidance ="Definition  of the rigidity values at which bacward ";
  guidance += "trajectories will be integrated in order to compute ";
  guidance += "asymptotic directions and rigidity cutoff"; 
  myRigidityFilterDir->SetGuidance(guidance);
  
  
  // parameters
  //------------
  
  G4UIparameter* coord_sys_param = 
               new G4UIparameter("Coordinate system ",'s',false);
  
  coord_sys_param->SetParameterCandidates("GEO GEOID GSM SM MAG GEI GSE");
  
  G4UIparameter* direction_coord_sys_param = 
               new G4UIparameter("Coordinate system ",'s',false);
  
  direction_coord_sys_param->SetParameterCandidates("GEO GEOID GSM SM MAG GEI GSE");
  
  G4UIparameter* mag_shell_coord_sys_param = 
               new G4UIparameter("Magnetic coordinate system ",'s',false);
  
  mag_shell_coord_sys_param->SetParameterCandidates("SM MAG"); 
  
  G4UIparameter* X_param= new G4UIparameter("X",'d',false);
  G4UIparameter* Y_param= new G4UIparameter("Y",'d',false);
  G4UIparameter* Z_param= new G4UIparameter("Z",'d',false);
  
  G4UIparameter* length_unit_param = new G4UIparameter("Length Unit",'s',false);
  length_unit_param->SetParameterCandidates("Re re RE km m");
  
  G4UIparameter* angle_unit_param = new G4UIparameter("Angle Unit",'s',false);
  angle_unit_param->SetParameterCandidates("degree deg rad radian mrad milliradian");
  
  G4UIparameter* altitude_param = new G4UIparameter("Altitude",'d',false);
  G4UIparameter* longitude_param = new G4UIparameter("Longitude",'d',false);
  G4UIparameter* latitude_param = new G4UIparameter("Latitude",'d',false);
  G4UIparameter* zenith_param = new G4UIparameter("Zenith",'d',false);
  G4UIparameter* azimuth_param = new G4UIparameter("Azimuth",'d',false);
  
  G4UIparameter* L_param = new G4UIparameter("L",'d',false);
  L_param->SetParameterRange(" L > 0.");
  L_param->SetGuidance("L shell parameter in Re");
  
  G4UIparameter* pitch_angle_param = new G4UIparameter("pitch_angle",'d',false);
  G4UIparameter* phi_angle_param = new G4UIparameter("phi_angle",'t',false);
  G4UIparameter* angle_unit_param1 = new G4UIparameter("Angle Unit",'s',false);
  angle_unit_param1->SetParameterCandidates("degree deg rad radian mrad milliradian");
  angle_unit_param1->SetDefaultValue("degree");
  
  // commands
  //----------
   

  SetPositionVectorCmd = new
             G4UIcommand("/MAGCOS/SOURCE/SetPositionVector",this);
  guidance ="Set the start position as a vector  in your selected ";
  guidance +="coordinate system";
  SetPositionVectorCmd->SetGuidance(guidance); 
  SetPositionVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetPositionVectorCmd->SetParameter(coord_sys_param);
  SetPositionVectorCmd->SetParameter(X_param);
  SetPositionVectorCmd->SetParameter(Y_param);
  SetPositionVectorCmd->SetParameter(Z_param);
  SetPositionVectorCmd->SetParameter(length_unit_param);
  
  SetPositionCmd = new
             G4UIcommand("/MAGCOS/SOURCE/SetPosition",this);
  guidance = "Set the altitude, latitude and longitude defining the start position ";
  guidance += "in your selected coordinate system";	     
  SetPositionCmd->SetGuidance(guidance);
  guidance ="[usage] /MAGCOS/SOURCE/SetPosition ";
  guidance +="CoordSys altitude length_ unit latitude longitude angle_unit";  
  SetPositionCmd->SetGuidance(guidance); 
  SetPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetPositionCmd->SetParameter(coord_sys_param);
  SetPositionCmd->SetParameter(altitude_param);
  SetPositionCmd->SetParameter(length_unit_param);
  SetPositionCmd->SetParameter(latitude_param);
  SetPositionCmd->SetParameter(longitude_param);
  SetPositionCmd->SetParameter(angle_unit_param);
  
  
  
  
  SetPositionOnDipoleMagneticShellCmd = new
             G4UIcommand("/MAGCOS/SOURCE/SetPositionOnDipoleMagneticShell",this);
  guidance="Set the L parameter (in Re), the latitude and longitude ";
  guidance += "defining the start position on a dipole magnetic shell "; 
  guidance += "in the SM or MAG coordinate system";
  SetPositionOnDipoleMagneticShellCmd
    ->SetGuidance(guidance);
  SetPositionOnDipoleMagneticShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(mag_shell_coord_sys_param);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(L_param);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(latitude_param);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(longitude_param);
  SetPositionOnDipoleMagneticShellCmd->SetParameter(angle_unit_param);
  
  
  SetDirectionVectorCmd = new
             G4UIcommand("/MAGCOS/SOURCE/SetDirectionVector",this);
  guidance  = "Set the start direction of the pirmary particle as a vector ";
  guidance  +="in your selected coordinate system";	     
  SetDirectionVectorCmd->SetGuidance(guidance);
  SetDirectionVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetDirectionVectorCmd->SetParameter(coord_sys_param);
  SetDirectionVectorCmd->SetParameter(X_param);
  SetDirectionVectorCmd->SetParameter(Y_param);
  SetDirectionVectorCmd->SetParameter(Z_param);
  
  
  SetDirectionCmd = new G4UIcommand("/MAGCOS/SOURCE/SetDirection",this);
  guidance ="The direction of the primary particle is defined ";
  guidance +="by a zenith and azimuth angle ";
  guidance +="from the vertical direction in your selected coordinate system. \n";
  guidance +="As the vertical direction depends on the ";
  guidance += "position, you should always defined the start position ";
  guidance += "before invoking this command";
  SetDirectionCmd->SetGuidance(guidance);
  SetDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  SetDirectionCmd->SetParameter(direction_coord_sys_param);
  SetDirectionCmd->SetParameter(zenith_param);
  SetDirectionCmd->SetParameter(azimuth_param);
  SetDirectionCmd->SetParameter(angle_unit_param);
  
  
  SetDirectionFromPitchAngleCmd= new
           G4UIcommand("/MAGCOS/SOURCE/SetDirectionFromPitchAngle",this);
  guidance ="Defines the direction of the particle from the pitch angle";
  SetDirectionFromPitchAngleCmd->SetGuidance(guidance );
  SetDirectionFromPitchAngleCmd->SetParameter(pitch_angle_param);
  SetDirectionFromPitchAngleCmd->SetParameter(phi_angle_param);
  SetDirectionFromPitchAngleCmd->SetParameter(angle_unit_param1);
  
  
  
  
  SetRigidityCmd =
    new G4UIcmdWithADoubleAndUnit("/MAGCOS/SOURCE/SetRigidity",this);
  SetRigidityCmd->SetGuidance("Define the rigidity of the primary particle");
  SetRigidityCmd->SetParameterName("Rigidity",false);
  SetRigidityCmd->SetUnitCategory("Electric potential");
  SetRigidityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
   
  AddValuesToRigidityVectorCmd = new
             G4UIcommand("/MAGCOS/RIGIDITYVECTOR/AddValues",this);
  AddValuesToRigidityVectorCmd
    ->SetGuidance("Increase the rigidity vector with new values in GV");
  AddValuesToRigidityVectorCmd
    ->SetGuidance("[usage] /MAGCOS/RIGIDITYVECTOR/AddValues   val0 dval nvalues" ); 
  AddValuesToRigidityVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  param = new G4UIparameter("val0",'d',false);
  AddValuesToRigidityVectorCmd->SetParameter(param);
  
  param = new G4UIparameter("dval",'d',false);
  param->SetParameterRange("dval<0.");
  AddValuesToRigidityVectorCmd->SetParameter(param);
  
  
  param = new G4UIparameter("nvalues",'i',true);
  AddValuesToRigidityVectorCmd->SetParameter(param);
  
  
  SetDefaultRigidityVectorCmd= 
   new G4UIcmdWithoutParameter("/MAGCOS/RIGIDITYVECTOR/SetDefault",this);
  SetDefaultRigidityVectorCmd->SetGuidance("Set default rigidity vector");
  SetDefaultRigidityVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  ResetRigidityVectorCmd= 
   new G4UIcmdWithoutParameter("/MAGCOS/RIGIDITYVECTOR/Reset",this);
  ResetRigidityVectorCmd->SetGuidance("Reset rigidity vector");
  ResetRigidityVectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
  
  
  cmd_name = "/MAGCOS/SOURCE/verbose";
  SetVerbosityCmd = new G4UIcmdWithAnInteger(cmd_name,this);
  guidance = "If verbose >0 the position, direction, and energy ";
  guidance += "of the primary particles are printed"; 
  SetVerbosityCmd->SetGuidance(guidance);
  SetVerbosityCmd->SetParameterName("verbosity", false);
  SetVerbosityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cmd_name = "/MAGCOS/SOURCE/BfieldAtPrimaryPosition";
  PrintBfieldAtPrimaryCmd = new G4UIcmdWithoutParameter(cmd_name,this);
  PrintBfieldAtPrimaryCmd->SetGuidance("Print the Bfield at the  primary position");
  PrintBfieldAtPrimaryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  
 
  
  
}
/////////////////////////////////////////////////////////////////////////////////
//
MAGCOSPrimaryMessenger::~MAGCOSPrimaryMessenger()
{ delete myStartGeoDir;
  delete myRigidityFilterDir;
 
  delete SetRigidityCmd;
  delete SetPositionCmd;
  delete SetPositionOnDipoleMagneticShellCmd;
  delete SetDirectionCmd;
  delete SetDirectionVectorCmd;
  delete SetPositionAndDirectionVectorCmd;
  delete SetPositionAndDirectionCmd;
  delete SetDirectionFromPitchAngleCmd;
  
  delete AddValuesToRigidityVectorCmd;
  delete SetDefaultRigidityVectorCmd;
  delete ResetRigidityVectorCmd;
  
  delete SetVerbosityCmd;
  delete PrintBfieldAtPrimaryCmd;
  

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ if ( command == SetPositionVectorCmd ){
  	const char* paramString=newValues;
        G4double x,y,z;
	G4String  coord_sys,unit_str;
	std::istringstream is((char*)paramString);
        is >> coord_sys >> x >> y >> z >> unit_str ;
        G4ThreeVector GEOposition =
	     G4ThreeVector(x,y,z)*G4UnitDefinition::GetValueOf(unit_str);
	G4bool test;
	test= myGenerator->SetPosition(coord_sys,GEOposition);
  }
  else if ( command == SetPositionCmd ){
  	const char* paramString=newValues;
        G4double altitude,latitude,longitude;
	G4String  coord_sys,length_unit,angle_unit;
	std::istringstream is((char*)paramString);
        is >> coord_sys >> altitude >> length_unit >> latitude >> 
	       longitude >> angle_unit ;
        altitude *=G4UnitDefinition::GetValueOf(length_unit);
	longitude *=G4UnitDefinition::GetValueOf(angle_unit);
	latitude *=G4UnitDefinition::GetValueOf(angle_unit);
	G4bool test;
	test= myGenerator->SetPosition(coord_sys,altitude,longitude,latitude);
  }  
  else if ( command == SetPositionOnDipoleMagneticShellCmd ){
  	const char* paramString=newValues;
        G4double L,latitude,longitude;
	G4String  coord_sys,angle_unit;
	std::istringstream is((char*)paramString);
        is >> coord_sys >> L >> latitude >> 
	       longitude >> angle_unit ;
	longitude *=G4UnitDefinition::GetValueOf(angle_unit);
	latitude *=G4UnitDefinition::GetValueOf(angle_unit);
	G4bool test;
	test=
            myGenerator->SetPositionOnDipoleMagneticShell(coord_sys,L,latitude,longitude);
  }     
        
  else if ( command == SetDirectionVectorCmd ){
  	const char* paramString=newValues;
        G4double x,y,z;
	G4String  coord_sys;
	std::istringstream is((char*)paramString);
        is >> coord_sys >> x >> y >> z ;
        G4ThreeVector GEODirection = G4ThreeVector(x,y,z);
	GEODirection/=GEODirection.mag();
	G4bool test;
	test= myGenerator->SetDirection(coord_sys,GEODirection);
  }
       
  else if ( command == SetDirectionCmd ){
    	const char* paramString=newValues;
        G4double azimuth,zenith;
	G4String  coord_sys,angle_unit;
	std::istringstream is((char*)paramString);
        is >> coord_sys >> zenith >> azimuth >> angle_unit ;
	zenith *=G4UnitDefinition::GetValueOf(angle_unit);
	azimuth *=G4UnitDefinition::GetValueOf(angle_unit);
	G4bool test;
	test= myGenerator->SetDirection(coord_sys,zenith,azimuth);
  }
  else if ( command == SetDirectionFromPitchAngleCmd){
  	const char* paramString=newValues;
        G4double pitch_angle,phi;
	G4String angle_unit;
	std::istringstream is((char*)paramString);
	is >> pitch_angle >> phi >> angle_unit;
	pitch_angle*=G4UnitDefinition::GetValueOf(angle_unit);
        phi*=G4UnitDefinition::GetValueOf(angle_unit);
	myGenerator->DefineDirectionByPitchAngle(pitch_angle,phi);
  
  
  }     
  
  else if( command == SetRigidityCmd)
           myGenerator->SetRigidity(SetRigidityCmd->GetNewDoubleValue(newValues));
  
  
  else if ( command == AddValuesToRigidityVectorCmd ){
     	const char* paramString=newValues;
        G4double val0,dval;
	G4int   nvalues;
        std::istringstream is((char*)paramString);
        is >> val0 >> dval >> nvalues;
        myGenerator->AddValuesToRigidityVector(nvalues,val0,dval);
  }
   
  else if( command == SetDefaultRigidityVectorCmd) 
                       myGenerator->SetDefaultRigidityVector();    
       
  else if( command == ResetRigidityVectorCmd) 
                       myGenerator->ResetRigidityVector();       
  else if( command == SetVerbosityCmd){
        G4int verbosity =SetVerbosityCmd->GetNewIntValue(newValues);
        myGenerator->SetVerbosity(verbosity); 
	SpaceCoordinateConvertor::getInstance()->SetVerbosity(verbosity);
  }
  
  else if( command == PrintBfieldAtPrimaryCmd)
   	myGenerator->PrintBfieldAtPrimary();
}
