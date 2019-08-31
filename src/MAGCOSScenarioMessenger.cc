// 9/22/2014: Hexc & Olesya - Added CLHEP namespace
//
#include"MAGCOSApplicationScenario.hh"
#include"MAGCOSScenarioMessenger.hh"
#include"MAGCOSPrimaryGeneratorAction.hh"
#include"G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//units
#include "G4UnitsTable.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

////////////////////////////////////////////////////////////////////////////////
//
MAGCOSScenarioMessenger::MAGCOSScenarioMessenger
            (MAGCOSApplicationScenario* ApplicationScenario )
{ G4String guidance;
  
  //directory 
  
  theApplicationScenario = ApplicationScenario;
  ScenarioDir= new G4UIdirectory("/MAGCOS/SCENARIO/");
  ScenarioDir->SetGuidance("Start the different applications");
  
  //parameters
  //----------

 
  G4UIparameter* coord_sys_param = 
               new G4UIparameter("Coordinate system ",'s',false);
  coord_sys_param->SetParameterCandidates("GEO GEOID GSM SM MAG GEI GSE");
  
  G4UIparameter* direction_coord_sys_param = 
               new G4UIparameter("Coordinate system ",'s',false);
  direction_coord_sys_param->SetParameterCandidates("GEO GEOID GSM SM MAG GEI GSE");
  
  G4UIparameter* mag_shell_coord_sys_param = 
               new G4UIparameter("Magnetic coordinate system ",'s',false);
  mag_shell_coord_sys_param->SetParameterCandidates("SM MAG"); 
  
  
  
  G4UIparameter* length_unit_param = new G4UIparameter("Length Unit",'s',false);
  length_unit_param->SetParameterCandidates("Re km m");
  
  G4UIparameter* angle_unit_param = new G4UIparameter("Angle Unit",'s',false);
  angle_unit_param->SetParameterCandidates("degree deg rad radian mrad milliradian");
  
  G4UIparameter* altitude_param = new G4UIparameter("Altitude",'d',false);
  
  G4UIparameter* longitude_param = new G4UIparameter("Longitude",'d',false);
  G4UIparameter* dlong_param = new G4UIparameter("delta longitude",'d',false);
  G4UIparameter* nLong_param = new G4UIparameter("number of longitudes",'i',false);
  
  G4UIparameter* latitude_param = new G4UIparameter("Latitude",'d',false);
  G4UIparameter* dlat_param = new G4UIparameter("delta latitude",'d',false);
  G4UIparameter* nLat_param = new G4UIparameter("number of latitudes",'i',false);
  
  G4UIparameter* zenith_param = new G4UIparameter("Zenith",'d',false);
  G4UIparameter* dZenith_param = new G4UIparameter("delta zenith",'d',false);
  G4UIparameter* nZenith_param = new G4UIparameter("number of zenith angles",'i',false);
  
  G4UIparameter* azimuth_param = new G4UIparameter("Azimuth",'d',false);
  G4UIparameter* dAzimuth_param = new G4UIparameter("delta azimuth",'d',false);
  G4UIparameter* nAzimuth_param = new G4UIparameter("number of azimuth angles",'i',false);
  
  G4UIparameter* L_param = new G4UIparameter("L",'d',false);
  L_param->SetParameterRange(" L > 0.");
  L_param->SetGuidance("L shell parameter in Re");
  
  G4UIparameter* output_file_param = new G4UIparameter("OutputFile name",'s',false);
  G4UIparameter* input_file_param = new G4UIparameter("InputFile name",'s',false);
  
  G4UIparameter* start_time_param = new G4UIparameter("start time",'d',false);
  G4UIparameter* dTime_param = new G4UIparameter("delta time",'d',false);
  G4UIparameter* nTime_param = new G4UIparameter("number of different times",'i',false);
  G4UIparameter*  time_unit = new G4UIparameter("Time Unit",'s',false);
  time_unit->SetParameterCandidates("s second hour  minute day");
  
  // Scenario commands
  //-----------------

  
  //Bline 
  
  BlineCmd = new G4UIcmdWithoutParameter("/MAGCOS/SCENARIO/TraceBline",this);
  guidance = "Compute and trace the magnetic field line passing through ";
  guidance += "the primary position";
  BlineCmd->SetGuidance(guidance);
  BlineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  //forward trajectory
  ParticleTrajectoryCmd = new G4UIcmdWithoutParameter("/MAGCOS/SCENARIO/TraceParticleTrajectory",this);
  ParticleTrajectoryCmd->SetGuidance("Compute and trace a particle trajectory");
  ParticleTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
   
  //bacward trajectory
  ReverseParticleTrajectoryCmd = new G4UIcmdWithoutParameter("/MAGCOS/SCENARIO/ReverseParticleTrajectory",this);
  guidance = "Compute and trace a particle trajectory bacward in time";
  ReverseParticleTrajectoryCmd->SetGuidance(guidance);
  ReverseParticleTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  //asymptotic directions
  ComputeRigidityFilterCmd =
    new G4UIcmdWithAString("/MAGCOS/SCENARIO/ComputeAsymptoticDirections",this);
  guidance ="Compute the asymptotic directions and cutoff rigidities ";
  guidance +="for the primary position and direction that you should have ";
  guidance +="previously defined by using the /MAGCOS/SOURCE commands";
  ComputeRigidityFilterCmd->SetGuidance(guidance);
  ComputeRigidityFilterCmd->SetParameterName("Output_file name",false);
  ComputeRigidityFilterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  /*// asymptotic directions
  //----------------------
  
  ComputeDirectionFilterCmd = new
             G4UIcommand("/MAGCOS/SCENARIO/AsymptoticDirVsDirection",this);
  ComputeDirectionFilterCmd
    ->SetGuidance("Computing of asymptotic direction and filter value for different zenith and azimuth");
  ComputeDirectionFilterCmd
    ->SetGuidance(
   "[usage] /MAGCOS/SCENARIO/AsymptoticDirVsDirection Rigidity unit CoordSys CosZen0 dCosZen nZen  Azim0 dAzim nAzim  OutputFile" ); 
  
  ComputeDirectionFilterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   param = new G4UIparameter("Rigidity",'d',false);
  param->SetDefaultValue("5.");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("unit",'s',false);
  param->SetDefaultValue("GV");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  
  
  param = new G4UIparameter("CoordSys",'s',false);
  param->SetDefaultValue("GEOID");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  
  
  
  param = new G4UIparameter("CosZen0",'d',false);
  param->SetDefaultValue("1.");
  ComputeDirectionFilterCmd->SetParameter(param);
  
 
  param = new G4UIparameter("dCosZen",'d',false);
  param->SetDefaultValue("-.1");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("nZen",'i',false);
  param->SetDefaultValue("10.");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("Azim0",'d',false);
  param->SetDefaultValue("0");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("dAzim",'d',false);
  param->SetDefaultValue("40");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("nAzim",'i',false);
  param->SetDefaultValue("9");
  ComputeDirectionFilterCmd->SetParameter(param);
  
  param = new G4UIparameter("OutputFile",'s',true);
  param->SetDefaultValue("Filter.txt");
  ComputeDirectionFilterCmd->SetParameter(param);*/
  
  
  
  // rigidity cutoff for different position
  
  RCutoffVsPositionCmd = new
             G4UIcommand("/MAGCOS/SCENARIO/RCutoffVsPosition",this);
  RCutoffVsPositionCmd
    ->SetGuidance("Computing cutoff rigidities  at different longitudes and latitudes");
  RCutoffVsPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsPositionCmd->SetParameter(coord_sys_param);
  RCutoffVsPositionCmd->SetParameter(altitude_param);
  RCutoffVsPositionCmd->SetParameter(length_unit_param);
  RCutoffVsPositionCmd->SetParameter(latitude_param);
  RCutoffVsPositionCmd->SetParameter(dlat_param);
  RCutoffVsPositionCmd->SetParameter(nLat_param);
  RCutoffVsPositionCmd->SetParameter(longitude_param);
  RCutoffVsPositionCmd->SetParameter(dlong_param);
  RCutoffVsPositionCmd->SetParameter(nLong_param);
  RCutoffVsPositionCmd->SetParameter(zenith_param);
  RCutoffVsPositionCmd->SetParameter(azimuth_param);
  RCutoffVsPositionCmd->SetParameter(angle_unit_param);
  RCutoffVsPositionCmd->SetParameter(output_file_param);
  
  // rigidity cutoff for different position on dipole-magnetic shell
  
  RCutoffVsPositionOnLShellCmd = new
             G4UIcommand("/MAGCOS/SCENARIO/RCutoffVsPositionOnLShell",this);
  guidance ="Computing of cutoff rigidities at different longitudes and latitudes ";	     
  guidance +="on a given dipole magnetic shell";	     
  RCutoffVsPositionOnLShellCmd->SetGuidance(guidance);
  RCutoffVsPositionOnLShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsPositionOnLShellCmd->SetParameter(mag_shell_coord_sys_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(L_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(latitude_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(dlat_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(nLat_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(longitude_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(dlong_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(nLong_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(zenith_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(azimuth_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(angle_unit_param);
  RCutoffVsPositionOnLShellCmd->SetParameter(output_file_param);
  
  // rigidity cutoff for spenvis trajectory
  
  RCutoffVsSpenvisTrajectoryCmd = new
             G4UIcommand("/MAGCOS/SCENARIO/RCutoffVsSpenvisTrajectory",this);
  RCutoffVsSpenvisTrajectoryCmd
    ->SetGuidance("Computing cutoff rigidities along a trajectory");
  RCutoffVsSpenvisTrajectoryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(input_file_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(output_file_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(zenith_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(azimuth_param);
  RCutoffVsSpenvisTrajectoryCmd->SetParameter(angle_unit_param);
  
  // rigidity cutoff for spenvis position grid
  
  RCutoffVsSpenvisPositionGridCmd = new
             G4UIcommand("/MAGCOS/SCENARIO/RCutoffVsSpenvisPositionGrid",this);
  RCutoffVsSpenvisPositionGridCmd
    ->SetGuidance("Computing cutoff rigidities along a trajectory");
  RCutoffVsSpenvisPositionGridCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(input_file_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(output_file_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(zenith_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(azimuth_param);
  RCutoffVsSpenvisPositionGridCmd->SetParameter(angle_unit_param);
 

  
  // rigidity cutoff for different directions
 
  RCutoffVsDirectionCmd = new
             G4UIcommand("/MAGCOS/SCENARIO/RCutoffVsDirection",this);
  guidance ="Computing of cutoff rigidities for different directions ";
  guidance +="defined by a grid of different zenith and azimuth angles";
  RCutoffVsDirectionCmd->SetGuidance(guidance);
  RCutoffVsDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsDirectionCmd->SetParameter(direction_coord_sys_param);
  RCutoffVsDirectionCmd->SetParameter(zenith_param);
  RCutoffVsDirectionCmd->SetParameter(dZenith_param);
  RCutoffVsDirectionCmd->SetParameter(nZenith_param);
  RCutoffVsDirectionCmd->SetParameter(azimuth_param);
  RCutoffVsDirectionCmd->SetParameter(dAzimuth_param);
  RCutoffVsDirectionCmd->SetParameter(nAzimuth_param);
  RCutoffVsDirectionCmd->SetParameter(output_file_param);
  
  
  
 // Cutoff for different times
 
  RCutoffVsTimeCmd = new
             G4UIcommand("/MAGCOS/SCENARIO/RCutoffVsTime",this);
  RCutoffVsTimeCmd->SetGuidance("Computing of cutoff rigidities  for different times");
  RCutoffVsTimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  RCutoffVsTimeCmd->SetParameter(start_time_param);
  RCutoffVsTimeCmd->SetParameter(dTime_param);
  RCutoffVsTimeCmd->SetParameter(nTime_param);
  RCutoffVsTimeCmd->SetParameter(time_unit);
  RCutoffVsTimeCmd->SetParameter(output_file_param);
  
  
  //AutoDetectionOfPenumbra
  SetAutoDetectionOfPenumbra = new 
      G4UIcmdWithABool("/MAGCOS/SCENARIO/AutomaticDetectionOfPenumbra",this);
  guidance =" If the parameter Auto is true the penumbra is detected ";
  guidance += " automatically for the case where cutoff rigidities ";
  guidance += "are computed in function of position, direction of ";
  guidance += "incidence and time ";
  SetAutoDetectionOfPenumbra->SetGuidance(guidance);  
  SetAutoDetectionOfPenumbra->SetParameterName("Auto",true);
  SetAutoDetectionOfPenumbra->SetDefaultValue(true);    
  SetAutoDetectionOfPenumbra->AvailableForStates(G4State_PreInit,G4State_Idle);
 
 
  //Command for registering results in SpenvisCSV file
  
  SetRegisterResultsInSpenvisCSVFileCmd = new 
      G4UIcmdWithABool("/MAGCOS/SCENARIO/RegisterResultsInSpenvisCSVFile",this);
  guidance =" If true (fase)l the results of cutoff rigidity ";
  guidance += "and asymptotic direction computation are (not) registered in a SpenvisCSV file";
  SetRegisterResultsInSpenvisCSVFileCmd->SetGuidance(guidance);  
  SetRegisterResultsInSpenvisCSVFileCmd->SetParameterName("Register",false);
  SetRegisterResultsInSpenvisCSVFileCmd->SetDefaultValue(true);    
  SetRegisterResultsInSpenvisCSVFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  SetSpenvisCSVFileNameCmd =
    new G4UIcmdWithAString("/MAGCOS/SCENARIO/SetSpenvisCSVFileName",this);
  guidance ="Set the the name of the Spenvis csv file";
  SetSpenvisCSVFileNameCmd->SetGuidance(guidance);
  SetSpenvisCSVFileNameCmd->SetParameterName("SpenvisCSV file  name",false);
  SetSpenvisCSVFileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  
  
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSScenarioMessenger::~MAGCOSScenarioMessenger()
{
}		  
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSScenarioMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if (command == BlineCmd) theApplicationScenario->Bline();
  else if (command == ParticleTrajectoryCmd) 
               theApplicationScenario->ParticleTrajectory();
  else if (command == ReverseParticleTrajectoryCmd) 
               theApplicationScenario->ReverseParticleTrajectory();	       
  else if ( command == ComputeRigidityFilterCmd) 
                  theApplicationScenario->ComputeRigidityFilter(newValues);
  else if ( command == ComputeDirectionFilterCmd ){
  	const char* paramString=newValues;
        G4double  cos_zen0,dcos_zen,azim0,dazim,rigidity;
	G4int   nzen,nazim;
	G4String OutputFile,CoordSys,unit;
        std::istringstream is((char*)paramString);
        is >> rigidity >> unit >> CoordSys>>cos_zen0 >> dcos_zen >> nzen >>
	azim0 >> dazim >> nazim>>OutputFile;
	rigidity *=G4UnitDefinition::GetValueOf(unit);
        theApplicationScenario->ComputeDirectionFilter
	          		(CoordSys,rigidity,cos_zen0,dcos_zen,nzen,
		   		 azim0*degree,dazim*degree,nazim,OutputFile);
  }		  			      
  else if ( command == RCutoffVsPositionCmd ){
  	const char* paramString=newValues;
        G4double  Alt,lat0,dlat,long0,dlong,zenith,azimuth;
	G4int   nlat,nlong;
	G4String OutputFile,CoordSys,LengthUnit,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> CoordSys>> Alt >> LengthUnit>>
	lat0 >> dlat >> nlat>>long0 >> dlong >> nlong>>
	zenith >> azimuth >> AngleUnit >> OutputFile;
	Alt *=G4UnitDefinition::GetValueOf(LengthUnit);
	lat0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlat*=G4UnitDefinition::GetValueOf(AngleUnit);
	long0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlong*=G4UnitDefinition::GetValueOf(AngleUnit);
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theApplicationScenario->RCutoffVsPosition
	          	(CoordSys,Alt,lat0,dlat,nlat,long0,dlong,nlong,zenith,
		   	 azimuth,OutputFile);
  }
  else if ( command == RCutoffVsPositionOnLShellCmd ){
  	const char* paramString=newValues;
       	G4double  L,lat0,dlat,long0,dlong,zenith,azimuth;
	G4int   nlat,nlong;
	G4String OutputFile,CoordSys,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> CoordSys>> L >>
	lat0 >> dlat >> nlat>>long0 >> dlong >> nlong>>
	zenith >> azimuth >> AngleUnit >> OutputFile;
	lat0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlat*=G4UnitDefinition::GetValueOf(AngleUnit);
	long0*=G4UnitDefinition::GetValueOf(AngleUnit);
	dlong*=G4UnitDefinition::GetValueOf(AngleUnit);
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theApplicationScenario->RCutoffVsPositionOnDipoleMagneticShell
	                 (CoordSys,L,lat0,dlat,nlat,long0,dlong,nlong,zenith,
		         azimuth,OutputFile);
  }
  else if ( command == RCutoffVsSpenvisTrajectoryCmd ){
  	const char* paramString=newValues;
        G4double  zenith,azimuth;
	G4String OutputFile,InputFile,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> InputFile>> OutputFile >>zenith >> azimuth >> AngleUnit;
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theApplicationScenario->RCutoffVsSpenvisTrajectory
	          	(InputFile, OutputFile, zenith, azimuth);
  }
  else if ( command == RCutoffVsSpenvisPositionGridCmd ){
  	const char* paramString=newValues;
        G4double  zenith,azimuth;
	G4String OutputFile,InputFile,AngleUnit;
        std::istringstream is((char*)paramString);
        is >> InputFile>> OutputFile >>zenith >> azimuth >> AngleUnit;
	zenith*=G4UnitDefinition::GetValueOf(AngleUnit);
	azimuth*=G4UnitDefinition::GetValueOf(AngleUnit);
        theApplicationScenario->RCutoffVsSpenvisPositionGrid
	          	(InputFile, OutputFile, zenith, azimuth);
  }
  
  else if ( command == RCutoffVsDirectionCmd ){
  	const char* paramString=newValues;
        G4double  zen0,dzen,azim0,dazim;
	G4int   nzen,nazim;
	G4String OutputFile,CoordSys;
        std::istringstream is((char*)paramString);
        is >> CoordSys>>zen0 >> dzen >> nzen>>azim0 >> dazim >> nazim>>OutputFile;
        theApplicationScenario->RCutoffVsDirection
	          (CoordSys,zen0*degree,dzen*degree,nzen,azim0*degree,dazim*degree,nazim,OutputFile);
  }		  
  
  else if ( command == RCutoffVsTimeCmd ){
  	const char* paramString=newValues;
       	G4double  time0,dtime;
	G4int   ntime;
	G4String TimeUnit,OutputFile;
        std::istringstream is((char*)paramString);
        is >> time0 >> dtime >> ntime;
	is >>TimeUnit;
	is>>OutputFile;
	G4double t_unit = G4UnitDefinition::GetValueOf(TimeUnit);
        time0 *=t_unit;
	dtime *=t_unit;
        theApplicationScenario->RCutoffVsTime(time0,dtime,ntime,OutputFile);
  }
  else if (command == SetAutoDetectionOfPenumbra){
  	G4bool aBool = SetAutoDetectionOfPenumbra->GetNewBoolValue(newValues);
	theApplicationScenario->SetAutomaticDetectionPenumbra(aBool);
  }
  else if (command == SetRegisterResultsInSpenvisCSVFileCmd){
  	G4bool aBool = SetRegisterResultsInSpenvisCSVFileCmd->GetNewBoolValue(newValues);
	theApplicationScenario->SetRegisterResultsInSpenvisFile(aBool);
  }
  else if ( command == SetSpenvisCSVFileNameCmd) 
                  theApplicationScenario->SetSpenvisFileName(newValues);
  
}












