#include"MAGCOSApplicationScenario.hh"
#include"MAGCOSScenarioMessenger.hh"
#include"MAGCOSPrimaryGeneratorAction.hh"
#include"MAGCOSEventAction.hh"
#include"MAGCOSTrackingAction.hh"
#include"MAGCOSEquationOfMotion.hh"
#include"MAGCOSMagneticField.hh"
#include"MAGCOSAnalysisManager.hh"
#include"DurationManager.hh"
#include"G4TransportationManager.hh"
#include"G4FieldManager.hh"
#include"G4UImanager.hh"
#include"G4ParticleDefinition.hh"
#include"G4Circle.hh"
#include"G4Event.hh"
#include"G4UnitsTable.hh"
#include"G4ios.hh"
#include"G4RunManager.hh"
#include"G4GeneralParticleSource.hh"
#include"G4Proton.hh"
#include"G4ChargedGeantino.hh"
#include"MAGCOSUnits.hh"
#include"SpaceCoordinateConvertor.hh"
#include"MAGCOSSpenvisManager.hh"
#include <time.h>
 
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSApplicationScenario::MAGCOSApplicationScenario()
{ theMessenger = new MAGCOSScenarioMessenger(this);
  RegisterAsymptoticDirection=false;
  AutomaticDetectionOfThePenumbra=true;
  theSpenvisManager = MAGCOSSpenvisManager::getInstance();
  RegisterResultsInSpenvisFile = false;
  SpenvisFileName = "spenvis_results.cvs";
  TotalDurationReached = false;
  
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSApplicationScenario::~MAGCOSApplicationScenario()
{ delete theMessenger;
}		  
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::BeginOfRunAction(const G4Run* )
{ LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
  
  DurationManager* theDurationManager = DurationManager::getInstance();
  TotalDurationReached = !theDurationManager->CheckDurationAtBeginOfRun();
  if (TotalDurationReached) {
  	G4RunManager::GetRunManager()->AbortRun(true);
	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  }
  
}  
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::EndOfRunAction(const G4Run* )
{
}
////////////////////////////////////////////////////////////////
void MAGCOSApplicationScenario::Bline()
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }	
 
  //get a pointer to the particle source 
  G4GeneralParticleSource* theParticleSource =
       ((MAGCOSPrimaryGeneratorAction*)
                  G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())->
	                                               GetParticleSource() ;
  
  G4ParticleDefinition* theOldParticleDefinition
                   =theParticleSource->GetParticleDefinition();
  theParticleSource->SetParticleDefinition(G4ChargedGeantino::ChargedGeantino());
 
 // get a pointer to the magnetic field  
  MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
     G4TransportationManager::GetTransportationManager()->
                              GetFieldManager()->GetDetectorField();
 
 //get a pointer to the equation of motion
  MAGCOSEquationOfMotion* theEquation =  theField->GetEquationOfMotion();
 
 
 
 //select the bline equation
  theEquation->SetEquationType("BFIELD_LINE");
 
 
 //spenvis csv file
  theSpenvisManager->InitialiseTrajectoryCSVBlock("Magnetic field line");
  
 
 
 //southward integration
  theEquation->SetReverseTimeMode(true);
  G4RunManager::GetRunManager()->BeamOn(1);
  
 //northward integration
  theEquation->SetReverseTimeMode(false);
  if (!TotalDurationReached) G4RunManager::GetRunManager()->BeamOn(1); 
  
  //Back to Lorentz Equation
  theEquation->SetEquationType("LORENTZ_MOTION");
  
  
  //Back to old particle definition
   
  if (theOldParticleDefinition)  
      theParticleSource->SetParticleDefinition(theOldParticleDefinition);
     
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::ParticleTrajectory()
{ 
 if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
 // get a pointer to the magnetic field  
  MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
     G4TransportationManager::GetTransportationManager()->
                              GetFieldManager()->GetDetectorField();
 
 //get a pointer to the equation of motion
  MAGCOSEquationOfMotion* theEquation =  theField->GetEquationOfMotion();
 
 
 
 //select the bline equation
  theEquation->SetEquationType("LORENTZ_MOTION"); 
  
 
 //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  thePrimaryGenerator->SetInitialiseTrajectoryCSVBlockForForwardCase(true);

 //integration
  theEquation->SetReverseTimeMode(false);
  G4RunManager::GetRunManager()->BeamOn(1);
  
  thePrimaryGenerator->SetInitialiseTrajectoryCSVBlockForForwardCase(false);
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::ReverseParticleTrajectory()
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  // get a pointer to the magnetic field  
  MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
     G4TransportationManager::GetTransportationManager()->
                              GetFieldManager()->GetDetectorField();
 
  //get a pointer to the equation of motion
  MAGCOSEquationOfMotion* theEquation =  theField->GetEquationOfMotion();
 
 
 
  //select the bline equation
  theEquation->SetEquationType("LORENTZ_MOTION");
  
  
  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  thePrimaryGenerator->SetInitialiseTrajectoryCSVBlockForBackwardCase(true);
 
 
  //integration
  theEquation->SetReverseTimeMode(true);
  G4RunManager::GetRunManager()->BeamOn(1);
  theEquation->SetReverseTimeMode(false);
  thePrimaryGenerator->SetInitialiseTrajectoryCSVBlockForBackwardCase(false);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::ComputeRigidityFilter(G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  RegisterAsymptoticDirection=true;
  
   //csv file
  if (RegisterResultsInSpenvisFile)
  	theSpenvisManager->InitialiseAsymptoticDirectionCSVFile(); 
 
  MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  theAnalysisManager->OpenAsymptoticDirectionFile(outputfile_name);
  ComputeRigidityCutoff();
  theAnalysisManager->CloseAsymptoticDirectionFile(Rc,Rm,Rs); 
			  
  if (RegisterResultsInSpenvisFile)
  		theSpenvisManager->SaveCSVFile(SpenvisFileName);
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
  RegisterAsymptoticDirection=false; 
 
}
///////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::ComputeDirectionFilter
                          (G4String CoordSys,G4double rigidity,G4double cos_zen0, 
			   G4double delta_cos_zen, G4int nzen,
                           G4double azim0, G4double delta_azim, G4int nazim,
			   G4String outputfile_name)
{ 
  if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  //Open the ascii file
  MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  theAnalysisManager->OpenAsymptoticDirVsDirFile(CoordSys,rigidity,outputfile_name);
  
 
  // get a pointer to the magnetic field  
  MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
               G4TransportationManager::GetTransportationManager()->
                                                 GetFieldManager()->GetDetectorField();
 
  //get a pointer to the equation of motion
  MAGCOSEquationOfMotion* theEquation =  theField->GetEquationOfMotion();
 
 
 
  //select the bline equation
  theEquation->SetEquationType("LORENTZ_MOTION");
  theEquation->SetReverseTimeMode(true);
 
 
  //pointer on MAGCOSTrackingAction
  MAGCOSTrackingAction* theTrackingAction =
    (MAGCOSTrackingAction*)
          G4RunManager::GetRunManager()->GetUserTrackingAction();
  theTrackingAction->SetRegisterLastPoint(true);
 
  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

  thePrimaryGenerator->SetRigidity(rigidity);
  
  for (int i=0; i<nzen;i++){
  	G4double zenith=std::acos( cos_zen0+double(i)*delta_cos_zen);
   	for (int j=0; j<nazim;j++){
      		G4double azimuth=azim0+double(j)*delta_azim;
       		if (thePrimaryGenerator->SetDirection(CoordSys,zenith,azimuth)){ 
          		G4RunManager::GetRunManager()->BeamOn(1);
          		theAnalysisManager
	     			->RegisterAsymptoticDirVsDir
	                            (zenith,azimuth,LastPositions[0],LastMomentaOnCharge[0],
				                    FilterValues[0]);
						    
	  		LastPositions.clear();
	  		LastMomentaOnCharge.clear();
	  		FilterValues.clear();   
	  	}
		else { 
	 		j=nazim;
	  		i=nzen;
	 	}          
       
    	}
  }
  theEquation->SetReverseTimeMode(false);
  theTrackingAction->SetRegisterLastPoint(false); 
  theAnalysisManager->CloseAsciiFile(); 
} 
///////////////////////////////////////////////////////////////////////////////
// 
void MAGCOSApplicationScenario::RCutoffVsPosition
                          (G4String CoordSys, G4double Altitude,
                           G4double lat0,G4double delta_lat, G4int nlat,
                           G4double long0, G4double delta_long, G4int nlong,
			   G4double zenith,G4double azimuth,
			   G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  //Open the ascii file
  MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  theAnalysisManager->OpenCutoffVsPositionFile(outputfile_name,CoordSys,
	                                       Altitude, zenith, azimuth);
  //test cvs
  if (RegisterResultsInSpenvisFile)
  	theSpenvisManager->InitialiseCutoffVsPosCSVFile(CoordSys,zenith, azimuth);
 
   
  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
  	(MAGCOSPrimaryGeneratorAction*)
          	G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	  

  //get a pointer to the user event  action and set DrawTrajewctory to false
  MAGCOSEventAction* theEventAction = 
  	(MAGCOSEventAction*)
          	G4RunManager::GetRunManager()->GetUserEventAction();	  
  G4bool draw_trajectory = theEventAction->GetDrawTrajectory();
  theEventAction->SetDrawTrajectory(false);	  
  
  MAGCOSMagneticField* theField =
		const_cast<MAGCOSMagneticField*>(
			dynamic_cast<const MAGCOSMagneticField*> 
				(G4TransportationManager::GetTransportationManager()
					->GetFieldManager()->GetDetectorField()));
	 
	
  for (int i=0; i<nlat;i++){
  	G4double latitude=lat0+double(i)*delta_lat;
   	for (int j=0; j<nlong;j++){
      		G4double longitude=long0+double(j)*delta_long;
       		G4double test = 
          		thePrimaryGenerator
	    			->SetPositionAndDirection(CoordSys,Altitude,longitude,latitude, 
	                              zenith,azimuth);
		G4ThreeVector GEOPosition = 
				thePrimaryGenerator->GetGEOPosition();
		
		std::vector<double> ILat;		
		ILat.push_back(SpaceCoordinateConvertor::getInstance()
					->ComputeGEODipoleInvariantLatitude("GEO",GEOPosition ));		      
		if (getenv("INVARIANT_LATITUDE")){
			std::vector<double> LMc = theField->ComputeMcIlwainLParameter(GEOPosition,89.999*degree);
			for (int i=0;i<2;i++){
				
				if (LMc[i] >= 1. ) ILat.push_back(std::acos(std::sqrt(1/LMc[i])));
				else ILat.push_back(-999999.*degree);
				
			}	
		}	
	
					      
			
				      
       		if (test){ 
         		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	  		else ComputeRigidityCutoff();
			if (!TotalDurationReached){
          			theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude,longitude,ILat);
	  			if (RegisterResultsInSpenvisFile){
					std::vector<double > values;
					values.clear();
					values.push_back(theField->GetReferenceDate().SpenvisModifiedJulianDay());
					values.push_back(Altitude/km);
					values.push_back(latitude/degree);
					values.push_back(longitude/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
				}
			}		
		}
       		else{
         		j=nlong;
	  		i=nlat;
	 	}
       }
  }
  //G4cout<<"SC1"<<std::endl;
  if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);
  //G4cout<<"SC2"<<std::endl;
  theAnalysisManager->CloseAsciiFile();   
  //G4cout<<"SC3"<<std::endl;

  theEventAction->SetDrawTrajectory(draw_trajectory);
  
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl;
		    
}


///////////////////////////////////////////////////////////////////////////////////////////
// 
void MAGCOSApplicationScenario::RCutoffVsSpenvisPositionGrid(G4String SpenvisCSVFileName,G4String outputfile_name,
  				    			     G4double zenith , 
				    			     G4double azimuth)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  
  //Read the position grid
  
  std::vector<double> altitude;
  std::vector<double>  latitude;
  std::vector<double>  longitude;
  
  
  G4bool test = theSpenvisManager
  			->ReadSpenvisPositionGridCSVFile(SpenvisCSVFileName,
							 altitude,
							 latitude,
							 longitude,
							 1);
   
  if (test) {
  	
	
	
	//Open the ascii file
  	MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  	theAnalysisManager->OpenCutoffVsSpenvisPositionGridFile(outputfile_name,
								"GEOID",
	                                        		zenith, 
								azimuth);
 	//test cvs
  	if (RegisterResultsInSpenvisFile)
  		theSpenvisManager->InitialiseCutoffVsPosCSVFile("GEOID",zenith, azimuth);
 
   
  	//get a pointer to the primary generator action
  	MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
  		(MAGCOSPrimaryGeneratorAction*)
          		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	  

  	//get a pointer to the user event  action and set DrawTrajewctory to false
 	MAGCOSEventAction* theEventAction = 
  		(MAGCOSEventAction*)
          		G4RunManager::GetRunManager()->GetUserEventAction();	  
  	 
	G4bool draw_trajectory = theEventAction->GetDrawTrajectory();
  	theEventAction->SetDrawTrajectory(false);	  
  
  	MAGCOSMagneticField* theField =
		const_cast<MAGCOSMagneticField*>(
			dynamic_cast<const MAGCOSMagneticField*> 
				(G4TransportationManager::GetTransportationManager()
					->GetFieldManager()->GetDetectorField()));
	 
	for (unsigned int i=0; i<altitude.size();i++){
	 	G4double test1 = thePrimaryGenerator
					->SetPositionAndDirection("GEOID",
							     	  altitude[i],
							          longitude[i],
							          latitude[i], 
	                              			          zenith,
							          azimuth);
		G4ThreeVector GEOPosition = 
				thePrimaryGenerator->GetGEOPosition();
	        std::vector<double> ILat;		
		ILat.push_back(SpaceCoordinateConvertor::getInstance()
					->ComputeGEODipoleInvariantLatitude("GEO",GEOPosition ));		      
		if (getenv("INVARIANT_LATITUDE")){
			std::vector<double> LMc = theField->ComputeMcIlwainLParameter(GEOPosition,89.999*degree);
			for (int i=0;i<2;i++){
				
				if (LMc[i] >= 1. ) ILat.push_back(std::acos(std::sqrt(1/LMc[i])));
				else ILat.push_back(-999999.*degree);
				
			}	
		}
		if (test1){ 
         		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	  		else ComputeRigidityCutoff();
          		if (!TotalDurationReached){
				theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude[i],longitude[i],ILat);
	  			if (RegisterResultsInSpenvisFile){
					std::vector<double > values;
					values.clear();
					values.push_back(theField->GetReferenceDate().SpenvisModifiedJulianDay());
					values.push_back(altitude[i]/km);
					values.push_back(latitude[i]/degree);
					values.push_back(longitude[i]/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
				}
			}		
		}
		
	}
	if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);						     
  	theAnalysisManager->CloseAsciiFile();
	theEventAction->SetDrawTrajectory(draw_trajectory);
	clock2=clock();
    	double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  	G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl; 
	
  }
  else {
  	G4cout<<"A problem occured when reading the CSV file."<<std::endl; 
	G4cout<<"The rigidity cutoffs will not be computed."<<std::endl;
	
  }
 
  
		    
}
///////////////////////////////////////////////////////////////////////////////////////////
// 
void MAGCOSApplicationScenario::RCutoffVsSpenvisTrajectory(G4String SpenvisCSVFileName,G4String outputfile_name,
  				    			     G4double zenith, 
				    			     G4double azimuth)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  
  //Read the position grid
  
  std::vector<double> altitude;
  std::vector<double>  latitude;
  std::vector<double>  longitude;
  std::vector<double>  smjd;
  
  
  G4bool test = theSpenvisManager
  			->ReadSpenvisTrajectoryCSVFile(SpenvisCSVFileName,
							 altitude,
							 latitude,
							 longitude,
							 smjd,
							 1);
   
  if (test) {
  	
	
	DateAndTime theDate = DateAndTime(2000,1,1,0,0,0);
	DateAndTime StartDate = DateAndTime(2000,1,1,0,0,0);
	StartDate.ConvertToSpenvisModifiedJulianDay(smjd[0]);
	

	//Open the ascii file
  	MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  	theAnalysisManager->OpenCutoffVsSpenvisPositionGridFile(outputfile_name,
								"GEOID",
	                                        		zenith, 
								azimuth);
 	//cvs file
  	if (RegisterResultsInSpenvisFile)
  		theSpenvisManager->InitialiseCutoffVsPosCSVFile("GEOID",
								zenith, 
								azimuth);
 
   
  	//get a pointer to the primary generator action
  	MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
  		(MAGCOSPrimaryGeneratorAction*)
          		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	  

  	//get a pointer to the user event  action and set DrawTrajewctory to false
 	MAGCOSEventAction* theEventAction = 
  		(MAGCOSEventAction*)
          		G4RunManager::GetRunManager()->GetUserEventAction();	  
  	 
	G4bool draw_trajectory = theEventAction->GetDrawTrajectory();
  	theEventAction->SetDrawTrajectory(false);	  
  
  	MAGCOSMagneticField* theField =
		const_cast<MAGCOSMagneticField*>(
			dynamic_cast<const MAGCOSMagneticField*> 
				(G4TransportationManager::GetTransportationManager()
					->GetFieldManager()->GetDetectorField()));
	 
	theField->SetStartDate(StartDate);
	G4cout<<altitude.size()<<std::endl;
	for (unsigned int i=0; i<altitude.size();i++){
		theDate.ConvertToSpenvisModifiedJulianDay(smjd[i]);
		G4double t =theDate.DifferenceInSeconds(StartDate);
		theField->SetTimeOfB(t);
	 	G4double test1 = thePrimaryGenerator
					->SetPositionAndDirection("GEOID",
							     	  altitude[i],
							          longitude[i],
							          latitude[i], 
	                              			          zenith,
							          azimuth);
		G4ThreeVector GEOPosition = 
				thePrimaryGenerator->GetGEOPosition();
	        std::vector<double> ILat;		
		ILat.push_back(SpaceCoordinateConvertor::getInstance()
					->ComputeGEODipoleInvariantLatitude("GEO",GEOPosition ));		      
		if (getenv("INVARIANT_LATITUDE")){
			std::vector<double> LMc = theField->ComputeMcIlwainLParameter(GEOPosition,89.999*degree);
			for (int i=0;i<2;i++){
				
				if (LMc[i] >= 1. ) ILat.push_back(std::acos(std::sqrt(1/LMc[i])));
				else ILat.push_back(-999999.*degree);
				
			}	
		}
		if (test1){ 
         		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	  		else ComputeRigidityCutoff();
			if (!TotalDurationReached){
          			theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude[i],longitude[i],ILat);
	  			if (RegisterResultsInSpenvisFile){
					std::vector<double > values;
					values.clear();
					values.push_back(smjd[i]);
					values.push_back(altitude[i]/km);
					values.push_back(latitude[i]/degree);
					values.push_back(longitude[i]/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
				}
			}		
		}
		
	}
	if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);						     
  	theAnalysisManager->CloseAsciiFile();
	theEventAction->SetDrawTrajectory(draw_trajectory);
	clock2=clock();
    	double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  	G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl; 
	
  }
  else {
  	G4cout<<"A problem occured when reading the CSV file."<<std::endl; 
	G4cout<<"The rigidity cutoffs will not be computed."<<std::endl;
	
  }
 
}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::RCutoffVsPositionOnDipoleMagneticShell
                          (G4String CoordSys, G4double L,
                           G4double lat0,G4double delta_lat, G4int nlat,
                           G4double long0, G4double delta_long, G4int nlong,
			   G4double zenith,G4double azimuth,
			   G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  
  clock1=clock();
  //Open the ascii file
  MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  theAnalysisManager->OpenCutoffVsPositionOnLShellFile
                       (outputfile_name,CoordSys, L, zenith, azimuth);



  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	  

  //get a pointer to the user event  action and set DrawTrajewctory to false
  MAGCOSEventAction* theEventAction = 
  	(MAGCOSEventAction*)
          	G4RunManager::GetRunManager()->GetUserEventAction();	  
  G4bool draw_trajectory = theEventAction->GetDrawTrajectory();
  theEventAction->SetDrawTrajectory(false);
  
  
  MAGCOSMagneticField* theField =
		const_cast<MAGCOSMagneticField*>(
			dynamic_cast<const MAGCOSMagneticField*> 
				(G4TransportationManager::GetTransportationManager()
					->GetFieldManager()->GetDetectorField()));
	 
		  

  for (int i=0; i<nlat;i++){
  	G4double latitude=lat0+double(i)*delta_lat;
   	for (int j=0; j<nlong;j++){
      		G4double longitude=long0+double(j)*delta_long;
       		G4double test = thePrimaryGenerator
	    				->SetPositionOnDipoleMagneticShell
	                        		(CoordSys,L,latitude,longitude);
		G4ThreeVector GEOPosition = 
				thePrimaryGenerator->GetGEOPosition();
		std::vector<double> ILat;		
		ILat.push_back(SpaceCoordinateConvertor::getInstance()
					->ComputeGEODipoleInvariantLatitude("GEO",GEOPosition ));		      
		if (getenv("INVARIANT_LATITUDE")){
			std::vector<double> LMc = theField->ComputeMcIlwainLParameter(GEOPosition,89.999*degree);
			for (int i=0;i<2;i++){
				
				if (LMc[i] >= 1. ) ILat.push_back(std::acos(std::sqrt(1/LMc[i])));
				else ILat.push_back(-999999.*degree);
				
			}	
		}			      
       		if (test){ 
         		thePrimaryGenerator
	                 	->SetDirection(CoordSys, zenith, azimuth);
	  		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	  		else ComputeRigidityCutoff();
			if (!TotalDurationReached){
          			theAnalysisManager
	     				->RegisterCutoffVsPosition(Rc,Rm,Rs,latitude,
						longitude,ILat);
			}		
	  	}
       		else {
         		j=nlong;
	  		i=nlat;
	 	}
       }
  }
  theAnalysisManager->CloseAsciiFile();   
  
  theEventAction->SetDrawTrajectory(draw_trajectory);
  
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl;
}
/////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::RCutoffVsDirection
                          (G4String CoordSys,G4double zen0, G4double delta_zen, G4int nzen,
                           G4double azim0, G4double delta_azim, G4int nazim,
			   G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  clock1=clock();
  
  //Open the ascii file
  MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  theAnalysisManager->OpenCutoffVsDirectionFile(outputfile_name,CoordSys);

 
  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
  	(MAGCOSPrimaryGeneratorAction*)
          	G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
 
  //initialise the cvs file
  
  if (RegisterResultsInSpenvisFile) {
  	G4ThreeVector theGEOPosition = thePrimaryGenerator->GetGEOPosition();
  	SpaceCoordinateConvertor* theCoordinateConvertor  = 
  					SpaceCoordinateConvertor::getInstance();  
  	G4double altitude,longitude,latitude;
  	if (CoordSys == "GEOID" )
  		theCoordinateConvertor
			->ComputeGEOIDCoordinatesFromGEOPosition(theGEOPosition,
                      						 altitude, 
								 longitude,
								 latitude);
  	else {
  		G4ThreeVector thePosition =
			theCoordinateConvertor->Transform(theGEOPosition,
                  			               	  "GEO",
                                		 	   CoordSys);
		altitude = thePosition.mag()-Re;
		latitude = 90.*degree -thePosition.theta();
		longitude = thePosition.phi();
		if (longitude > 180. *degree) longitude = longitude -360.*degree;  					   
  	}
  	theSpenvisManager->InitialiseCutoffVsDirCSVFile(CoordSys, 
   				     altitude,
				     latitude,
				     longitude);
  }  
  
  
  //get a pointer to the user event  action and set DrawTrajewctory to false
  MAGCOSEventAction* theEventAction = 
  	(MAGCOSEventAction*)
          	G4RunManager::GetRunManager()->GetUserEventAction();	  
  G4bool draw_trajectory = theEventAction->GetDrawTrajectory();
  theEventAction->SetDrawTrajectory(false);

  for (int i=0; i<nzen;i++){
  	G4double zenith=zen0+double(i)*delta_zen;
   	for (int j=0; j<nazim;j++){
      		G4double azimuth=azim0+double(j)*delta_azim;
		
       		G4double test = 
              		thePrimaryGenerator->SetDirection(CoordSys,zenith,azimuth);
       		
		if (test){ 
            		if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
	    		else ComputeRigidityCutoff();
			if (!TotalDurationReached){
	  			theAnalysisManager
	     				->RegisterCutoffVsDirection(Rc,Rm,Rs,zenith,azimuth);
	  			
				if (RegisterResultsInSpenvisFile){
					std::vector<double > values;
					values.clear();
					values.push_back(zenith/degree);
					values.push_back(azimuth/degree);
					values.push_back(Rm);
					values.push_back(Rc);
					values.push_back(Rs);
					theSpenvisManager->WriteData(values);
				}
			}	
		}
		else{ 
	 		j=nazim;
	  		i=nzen;
	 	}          
       
       }
  }
  theAnalysisManager->CloseAsciiFile();
  if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);
  
  theEventAction->SetDrawTrajectory(draw_trajectory);
  
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl; 
}  
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::RCutoffVsTime
                         (G4double time0, G4double delta_t, G4int ntime, 
			  G4String outputfile_name)
{ if (TotalDurationReached) {
  	G4cout<< "The maximum allowed total duration is reached no run will be executed anymore"<<std::endl;
  	return;
  }
  clock_t clock1,clock2;
  clock1=clock();
  
  //Open the ascii file
  MAGCOSAnalysisManager* theAnalysisManager= MAGCOSAnalysisManager::getInstance();
  theAnalysisManager->OpenCutoffVsTimeFile(outputfile_name);

   
   //initialise the cvs file
  
  if (RegisterResultsInSpenvisFile) { 
  
  //get a pointer to the primary generator action
  	MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
  		(MAGCOSPrimaryGeneratorAction*)
          		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
 
  	G4ThreeVector theGEOPosition = thePrimaryGenerator->GetGEOPosition();
	G4ThreeVector theGEODirection =  thePrimaryGenerator->GetGEODirection();
  	G4double altitude,longitude,latitude, zenith,azimuth;
	altitude = theGEOPosition.mag()-Re;
	latitude = 90.*degree - theGEOPosition.theta();
	longitude = theGEOPosition.phi();
	if (longitude > 180.*degree) longitude += -360.*degree; 
	theGEODirection.rotateZ(-theGEOPosition.phi());
	theGEODirection.rotateY(-theGEOPosition.theta());
	theGEODirection=-theGEODirection;
	zenith = theGEODirection.theta();
	azimuth = theGEODirection.phi();
	
	theSpenvisManager->InitialiseCutoffVsTimeCSVFile("GEO", 
   				     			altitude,
				     			latitude,
				     			longitude,
							zenith,
							azimuth);
  }  
   
  // get a pointer to the magnetic field  
  MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
  G4TransportationManager::GetTransportationManager()->
                              GetFieldManager()->GetDetectorField();
 
  //get a pointer to the user event  action and set DrawTrajewctory to false
  MAGCOSEventAction* theEventAction = 
  	(MAGCOSEventAction*)
          	G4RunManager::GetRunManager()->GetUserEventAction();	  
  G4bool draw_trajectory = theEventAction->GetDrawTrajectory();
  theEventAction->SetDrawTrajectory(false);
  
  for (int i=0; i<ntime;i++){
  	G4double time=time0+double(i)*delta_t;
   	theField->SetTimeOfB(time/s);
   	//theField->PrintStormParameter();
   	if (AutomaticDetectionOfThePenumbra) ComputeCutoff();
   	else ComputeRigidityCutoff();
	if (!TotalDurationReached){
   		theAnalysisManager
	     		->RegisterCutoffVsTime(Rc,Rm,Rs,time);
		if (RegisterResultsInSpenvisFile){
				std::vector<double > values;
				values.clear();
				values.push_back(theField->GetReferenceDate().SpenvisModifiedJulianDay());
				values.push_back(Rm);
				values.push_back(Rc);
				values.push_back(Rs);
				theSpenvisManager->WriteData(values);
		}
	}	
  }

  theAnalysisManager->CloseAsciiFile();
  if (RegisterResultsInSpenvisFile) theSpenvisManager->SaveCSVFile(SpenvisFileName);
  theEventAction->SetDrawTrajectory(draw_trajectory);
  
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl;   

}
///////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::RegisterTrackLastPoint(G4ThreeVector position, 
                             G4ThreeVector momentum,G4double filter_value)
{ LastPositions.push_back(position);
  LastMomentaOnCharge.push_back(momentum);
  FilterValues.push_back(G4int(filter_value));
 
  if (RegisterAsymptoticDirection){
    	MAGCOSAnalysisManager::getInstance()
                ->RegisterAsymptoticDirection(position,momentum,
		                                  G4int(filter_value));
  	if (RegisterResultsInSpenvisFile) 
	    theSpenvisManager->RegisterAsymptoticDirection(position,momentum,
		                                           G4int(filter_value));
  
  } 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::ComputeRigidityCutoff()
{ if (TotalDurationReached) return;
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();

  // get a pointer to the magnetic field  
  MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
  	G4TransportationManager::GetTransportationManager()->
                              GetFieldManager()->GetDetectorField();
 
  //Set the  equation of motion
  MAGCOSEquationOfMotion* theEquation =  theField->GetEquationOfMotion();
  theEquation->SetEquationType("LORENTZ_MOTION");
  theEquation->SetReverseTimeMode(true);
 
  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
  	(MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  thePrimaryGenerator->SelectTypeOfPrimaries("RigidityFilter");
  thePrimaryGenerator->ResetRigidityIndex();
 
  //pointer on MAGCOSTrackingAction
  MAGCOSTrackingAction* theTrackingAction =
   	(MAGCOSTrackingAction*)
          	G4RunManager::GetRunManager()->GetUserTrackingAction();
  theTrackingAction->SetRegisterLastPoint(true);
 
  //Compute trajectory and rigidity cutoff

  G4int n_particles = thePrimaryGenerator->GetNumberOfRigidity();  
  if (n_particles> 0){
    	G4RunManager::GetRunManager()->BeamOn(n_particles);
    
    	// RM RS and RC computation
    	G4double Rss=0,Rmm=0;
    	G4double Rcc=0.;
    	G4double last_rigidity=100.,fac=0.,last_filter_value=1.;
    	for (int j=0; j<n_particles;j++){
      		G4double rigidity = LastMomentaOnCharge[j].mag()/1000.;
       		G4int filter_value = FilterValues[j];
		if (filter_value <0) filter_value =0;
       		if ( (fac==0.) && filter_value==0){
        		fac=1.;
         		Rmm=last_rigidity;
         	}
       		Rcc += fac*filter_value*(rigidity-last_rigidity);
       		last_rigidity = rigidity;
      		if ( filter_value==1) Rss=rigidity;
      		last_filter_value=filter_value;
      	}
     
     	Rcc+=Rmm;
     	Rc=Rcc;
     	Rs=Rss;
     	Rm=Rmm;
      
  }
  else G4cout<<"The rigidity vector is empty"<<G4endl;
  if (Rm==Rc) Rs=Rc;
  
 
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
 
 
  //go back to normal conditions
 
  theEquation->SetReverseTimeMode(false);
  thePrimaryGenerator->SelectTypeOfPrimaries("Standard");
  theTrackingAction->SetRegisterLastPoint(false);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSApplicationScenario::ComputeCutoff()
{ 
 //In this method the cutoff rigidity is computed in the following way:
 // a) A start rigidity Sc is computed from the Störmer formula by considering the 
 //   geomagnetic  dipole
 // b) It is looked in the following rigidity series [Sc,Sc+1.GV,Sc+2.GV,....]
 // 	for the first rigidity with an allowed trajectory. This rigidity is called R1st .
 // c)It is looked in the serie [R1st-0.01GV, R1st-0.02GV,...,R1st-n*0.01GV,.... ]
 //    for the first rigidity with a forbidden trajectory. We call this rigidity R1stForb.
 // d)We compute trajectories for [R1st +0.01GV, R1st+0.02GV,...,R1st+n*0.01GV,.... ]
 //   till we have find a 1GV  band of rigidity with only allowed trajectories. Note that
 // the band [R1stForb+0.01GV, R1st] is considered as allowed and if a serie of rigidity
 // is found directly above [R1stForb+0.01GV, R1st] such that together they form a 1GV
 // band of allowed trajectories the 1GV band is considered as being found.
 //  We define Rm as the lowest rigidity in the 1GV allowed band.
 // e)It is looked in the serie [R1stForb, R1stForb-0.01GV,...,R1stForb-n*0.01GV,.... ]
  //   to the first 1GV band with all trajectories forbidden. Rs is defined by the 
  // highest rigidity of this band +0.01 GV.
  // f) Nallowed the number of allowed trajectory betwenn
  //   Rm and Rs-0.01 GV is computed during a) to e). 
  //   Rc =Rm- Nallowed*0.01*GV
   

  if (TotalDurationReached) return;
  
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();

  // get a pointer to the magnetic field  
  MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
     G4TransportationManager::GetTransportationManager()->
                              GetFieldManager()->GetDetectorField();
 
  //Set the  equation of motion
  MAGCOSEquationOfMotion* theEquation =  theField->GetEquationOfMotion();
  theEquation->SetEquationType("LORENTZ_MOTION");
  theEquation->SetReverseTimeMode(true);
 
  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
  (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	  
	  
  //GEOPosition GEOdirection
  G4ThreeVector GEOPosition = thePrimaryGenerator->GetGEOPosition();
  G4ThreeVector GEODirection = thePrimaryGenerator->GetGEODirection();

  //GEODipole
  G4double DipoleB0=theField->GetDipoleB0();
  G4double DipolePhi=theField->GetDipolePhi();
  G4double DipoleTheta=theField->GetDipoleTheta();
  G4ThreeVector DipoleSchift=theField->GetDipoleSchift();

  // compute the Störmer cutoff in a geodipole field 
  //------------------------------------------------	  
  G4ThreeVector PositionInDipole=GEOPosition-DipoleSchift;
  PositionInDipole.rotateZ(-DipolePhi).rotateY(-DipoleTheta);
  G4double phi=PositionInDipole.phi();
  G4double r=PositionInDipole.r()/re;
  G4double cos_geomagnetic_lat=sin(PositionInDipole.theta());
  G4double cos3=cos_geomagnetic_lat*cos_geomagnetic_lat*cos_geomagnetic_lat;
  G4double cos4=cos3*cos_geomagnetic_lat;
  G4ThreeVector EastDirection=G4ThreeVector(-sin(phi),cos(phi),0); 
  G4ThreeVector DirectionInDipole=GEODirection;
  DirectionInDipole.rotateZ(-DipolePhi).rotateY(-DipoleTheta);
  G4double cos_gam=-EastDirection.dot(DirectionInDipole);
 
  // G4cout<<cos_gam<<" cos gam"<<endl;
  G4double c=2.99792458e8*m/s;
  G4double term=1.+std::sqrt(1.-cos3*cos_gam);
  G4double SCutoff=DipoleB0*re*c*cos4/(r*r*term*term);
  SCutoff=G4int(SCutoff/(0.01*GV))*(0.01*GV);
  SCutoff = std::max(SCutoff, 0.01*GV);
  
  //pointer on MAGCOSTrackingAction
  MAGCOSTrackingAction* theTrackingAction =
       (MAGCOSTrackingAction*)
               G4RunManager::GetRunManager()->GetUserTrackingAction();
  theTrackingAction->SetRegisterLastPoint(true);
  
  
  G4int nmax=100;
  G4int n=0;
  G4double rigidity=SCutoff;
  G4double first_rigidity=SCutoff;
  G4double rig_min = .01*GV;
  G4double rig_max = std::max(100. * GV, SCutoff+10.*GV);
  G4int n_forbidden=0;
  G4int n_allowed=0;
  //G4cout<<SCutoff/GV<<" GV SCutoff"<<std::endl;
  
 //find a first allowed trajectory
 //---------------------------
  first_rigidity-=1.*GV;
  G4int filter_value= 0;
  while (filter_value ==0){
  	first_rigidity+= 1.*GV;
	thePrimaryGenerator->SetRigidity(first_rigidity);
      	G4RunManager::GetRunManager()->BeamOn(1);
	if (TotalDurationReached) return;
      	filter_value = FilterValues[0];
	if (filter_value <0) filter_value =0;
  	
  }
 // G4cout<<"first "<<first_rigidity/GV<<std::endl;
  Rm=first_rigidity;
  
  
 //Find the first forbidden trajectory below first_rigidity 
  G4double first_forbid_rigidity=first_rigidity;
  while (filter_value ==1 &&  first_forbid_rigidity>= rig_min+0.01*GV){
  	first_forbid_rigidity-=0.01*GV;
	thePrimaryGenerator->SetRigidity(first_forbid_rigidity);
	G4RunManager::GetRunManager()->BeamOn(1);
	if (TotalDurationReached) return;
      	filter_value = FilterValues[0];
	if (filter_value <0) filter_value =0;
	
  }
  if (filter_value ==0) {
  	Rm=first_forbid_rigidity+0.01*GV;
  }
  else {
  	Rm=first_forbid_rigidity;
  	Rm/=GV;
  	Rc=Rm;
  	Rs=Rm;
  	LastPositions.clear();
 	 LastMomentaOnCharge.clear();
  	FilterValues.clear();
 
 
  	//go back to normal 
  	theEquation->SetReverseTimeMode(false);
  	thePrimaryGenerator->SelectTypeOfPrimaries("Standard");
  	theTrackingAction->SetRegisterLastPoint(false);
	
	if (Rm<0.015 *GV){
		Rm=0.;
		Rc=0.;
		Rs=0.;
  	}
	return;
  }
  Rs=Rm;
  
  
  // find first 1GV band of allowed trajectories above first_forbid_rigidity 
  //------------------------------------------------------------------------ 
  
  n=int((first_rigidity-Rm)/(0.01*GV))+1;
  // G4cout<<"Rm "<<Rm/GV<<std::endl;
  // G4cout<<"n "<<n<<std::endl;
  rigidity= first_rigidity;
  
  while (n< nmax && rigidity <= rig_max-.01*GV) {
     	rigidity += .01*GV;
     	// G4cout<<rigidity<<endl;
      	thePrimaryGenerator->SetRigidity(rigidity);
      	G4RunManager::GetRunManager()->BeamOn(1);
	if (TotalDurationReached) return;
      	filter_value = FilterValues[0];
	if (filter_value <0) filter_value =0;
      	n++;
     
    	// G4cout<<filter_value<<" filter"<<std::endl;
      
      	if (filter_value == 0 ){ 
        	n_forbidden+=1;
		Rm=rigidity + 0.01*GV;
	    	n=0;
	}      		 	 
  }
  
  // G4cout<<"Rm "<<Rm/GV<<std::endl;
  Rs=first_forbid_rigidity+0.01*GV;
  
  
  //find first 1GV band of forbidden trajectories below first_forbid_rigidity
  //---------------------------------
  
  n_allowed=int ((Rm- first_forbid_rigidity)/(0.01*GV));
   //we add 1 at n_forbidden because the lowest computed rigidity is forbidden
  n_forbidden+=1;
  n_allowed-=n_forbidden;
  //G4cout<<"n_allowed "<<n_allowed<<std::endl;
 
  n=1;
  rigidity=first_forbid_rigidity;
  while (n< 100 && rigidity >= rig_min + 0.01*GV ){ 
     	rigidity -= .01*GV;
      	thePrimaryGenerator->SetRigidity(rigidity);
      	G4RunManager::GetRunManager()->BeamOn(1);
	if (TotalDurationReached) return;
      	G4int filter_value = FilterValues[0];
	if (filter_value <0) filter_value =0;    
      	n++;
      	if (filter_value == 1 ){ 
        	n=0;
	  	n_allowed+=1; 
          	Rs=rigidity;
		//G4cout<<"new Rs :"<<Rs<<std::endl;
	}	 	 
  }
  //G4cout<<"n_allowed "<<n_allowed<<std::endl;
  //G4cout<<Rm/GV<<std::endl;
  if (Rm<0.015 *GV){
	Rm=0.;
	Rc=0.;
	Rs=0.;
  }
  
  else {
  	Rc=Rm-n_allowed*.01*GV; 
  }
  Rc/=GV;
  Rs/=GV;
  Rm/=GV;
  
  if (Rm==Rs) Rc=Rs;
  
  LastPositions.clear();
  LastMomentaOnCharge.clear();
  FilterValues.clear();
 
 
  //go back to normal 
  theEquation->SetReverseTimeMode(false);
  thePrimaryGenerator->SelectTypeOfPrimaries("Standard");
  theTrackingAction->SetRegisterLastPoint(false);
  
  
  
  
}



