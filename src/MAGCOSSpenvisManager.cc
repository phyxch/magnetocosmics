#include "MAGCOSSpenvisManager.hh"
#include "MAGCOSPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "SpaceCoordinateConvertor.hh"
#include "G4Step.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include "MAGCOSUnits.hh"
#include "MAGCOSMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4GeneralParticleSource.hh"

#ifdef  USE_UNILIB
#include"unilib_c_pub.h"
#endif

MAGCOSSpenvisManager* MAGCOSSpenvisManager::instance = 0;
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSSpenvisManager::MAGCOSSpenvisManager()  
{ theSpenvisCSVCollection = 0;
  theSpenvisCSVFiles.clear();
  RegisterParticleTrajectory = false;
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSSpenvisManager::~MAGCOSSpenvisManager() 
{ if (theSpenvisCSVCollection) delete  theSpenvisCSVCollection;
  
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSSpenvisManager* MAGCOSSpenvisManager::getInstance()
{
  if (instance == 0) instance = new MAGCOSSpenvisManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::WriteData(std::vector< double> values,SpenvisCSV* aBlock)
{ //G4cout<<values.size()<<std::endl;
  if (!aBlock) theSpenvisCSVFiles[0].AddDataRow(values); 
  else aBlock->AddDataRow(values);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::WriteBfieldInformation(SpenvisCSV* aBlock)
{ SpenvisCSV* theBlock;
  if (!aBlock) theBlock = &(theSpenvisCSVFiles[0]);
  else theBlock =aBlock;
  const MAGCOSMagneticField* theField =
		dynamic_cast<const MAGCOSMagneticField*> 
			(G4TransportationManager::GetTransportationManager()
					->GetFieldManager()->GetDetectorField());

 
#ifdef  USE_UNILIB
  if ( theField->GetUseUnilibModel() && theField->Getul_ifail() >= 0){
 		theBlock->AddMetaVariableStr("BFIELD_MODEL", "UNILIB");
	
 
  }
  else {
#endif
	G4String Bintern = theField->GetIntFieldModelName();
		
	theBlock->AddMetaVariableStr("BFIELD_MODEL", "MAGCOS");
	//Date of reference for the magnetic field
	DateAndTime theDate = theField->GetReferenceDate();
	AddMetaVariableToCSVFile(G4String("YEAR"),
				 G4String(""),
				 theDate.year,
				 G4String("%4.0f"),theBlock);
	AddMetaVariableToCSVFile(G4String("MONTH"),
				 G4String(""),
				 theDate.month,
				 G4String("%2.0f"),theBlock);
	AddMetaVariableToCSVFile(G4String("DAY"),
				 G4String(""),
				 theDate.day,
				 G4String("%2.0f"),
				 theBlock);
	AddMetaVariableToCSVFile(G4String("HOUR"),
				 G4String(""),
				 theDate.hour,
				 G4String("%2.0f"),
				 theBlock);
	AddMetaVariableToCSVFile(G4String("MINUTE"),
			 	 G4String(""),
				 theDate.min,
				 G4String("%2.0f"),
				 theBlock);
	AddMetaVariableToCSVFile(G4String("SECOND"),
					 G4String(""),
					 theDate.sec,
				         G4String("%2.0f"),
				         theBlock);
	AddMetaVariableToCSVFile(G4String("MSECOND"),
			         G4String(""),
			         theDate.msec,
				 G4String("%4.0f"),
				         theBlock);
		
		
	//Internal field			 		 			 			 				 
	theBlock->AddMetaVariableStr("BFIELD_INT", Bintern);
	if (Bintern == "IGRF"){
		AddMetaVariableToCSVFile(G4String("BFIELD_NMIGRF"),
					 G4String(""),
					 double(theField->Getnm_igrf()),
					 G4String("%2.0f"),
					 theBlock);
		}
	if (Bintern == "DIPOLE") {
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPB0"),
					 G4String("nT"),
					 theField->GetDipoleB0()/nanotesla,
					 G4String("%5.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPTH"),
					 G4String("degree"),
					 theField->GetDipoleTheta()/degree,
					 G4String("%4.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPPHI"),
					 G4String("degree"),
					 theField->GetDipolePhi()/degree,
					 G4String("%4.2f"),
					 theBlock);
		G4ThreeVector DipSchift =theField->GetDipoleSchift();
		std::vector<double> dip_schift;
		dip_schift.push_back(DipSchift.x()/km);
		dip_schift.push_back(DipSchift.y()/km);
		dip_schift.push_back(DipSchift.z()/km);
		theBlock->AddMetaVariable(G4String("BFIELD_DIPSCHIFT"),
 					  dip_schift, 
					  "%6.2f", G4String("km"));				 
			
	}
	G4String Bextern = theField->GetExtFieldModelName();
	theBlock->AddMetaVariableStr("BFIELD_EXT", Bextern);
		
	//extern field parameter
	if (Bextern != "NOFIELD") {
		AddMetaVariableToCSVFile(G4String("BFIELD_DIPTILT"),
					 G4String("degree"),
					 theField->GetDipolePS()/degree,
					 G4String("%4.2f"),
					 theBlock);
	}
	if (Bextern == "TSY89") {
		AddMetaVariableToCSVFile(G4String("BFIELD_KPIOPT"),
					 G4String(""),
					 theField->Getiopt(),
					 G4String("%1.0f"),
					 theBlock);
	}
	if (Bextern == "TSY96" || Bextern == "TSY2001" ) {
		AddMetaVariableToCSVFile(G4String("BFIELD_PDYN"),
					 G4String("Pa"),
					 theField->GetPdyn(),
					 G4String("%4.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_DST"),
					 G4String("nT"),
					 theField->GetDst()/nanotesla,
					 G4String("%6.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_IMFY"),
					 G4String("nT"),
					 theField->GetImfy()/nanotesla,
					 G4String("%5.2f"),
					 theBlock);
		AddMetaVariableToCSVFile(G4String("BFIELD_IMFZ"),
					 G4String("nT"),
					 theField->GetImfz()/nanotesla,
				 	 G4String("%5.2f"),
					 theBlock);
		if  (Bextern == "TSY2001") {
			AddMetaVariableToCSVFile(G4String("BFIELD_G1"),
					       	 G4String(""),
					       	 theField->GetG1(),
						 G4String("%6.2f"),
						 theBlock);
			AddMetaVariableToCSVFile(G4String("BFIELD_G2"),
					       	 G4String(""),
					       	 theField->GetG2(),
						 G4String("%6.2f"),
						 theBlock);		 
			
		}
						 
						 			 
	}
		
		
		
#ifdef  USE_UNILIB
 }
#endif

}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::SaveCSVFile(G4String name)
{ G4cout<<"CSV1"<<std::endl;
  theSpenvisCSVCollection = new SpenvisCSVCollection();
  for (unsigned int i=0; i<theSpenvisCSVFiles.size(); i++){
  	theSpenvisCSVCollection->AddCSVBlock(block_names[i],
						theSpenvisCSVFiles[i]);
	
  }
  G4cout<<"CSV2"<<std::endl;
  theSpenvisCSVCollection->OutputCollection(name);
  G4cout<<"CSV3"<<std::endl;
  for (unsigned int i=0; i<theSpenvisCSVFiles.size(); i++){
  	theSpenvisCSVCollection->DelCSVBlock(block_names[i]);
	
  }
  G4cout<<"CSV4"<<std::endl;
/*  theSpenvisCSVFiles.clear();
  block_names.clear(); */
  delete theSpenvisCSVCollection;
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::ClearCSVFiles()
{ theSpenvisCSVFiles.clear();
  block_names.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::InitialiseAsymptoticDirectionCSVFile()
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  
  //get a pointer to the primary generator action
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  G4double GEOIDalt, GEOIDlat, GEOIDlong, GEOIDzenith, GEOIDazimuth; 
  thePrimaryGenerator->GetGEOIDPosition( GEOIDalt,
  					 GEOIDlat,
					 GEOIDlong);
  					
  thePrimaryGenerator->GetGEOIDDirection( GEOIDzenith, GEOIDazimuth);  
  G4ThreeVector GEOPosition = thePrimaryGenerator->GetGEOPosition();
  std::vector<double> pos ;
  pos.push_back(GEOPosition.x()/km);
  pos.push_back(GEOPosition.y()/km);
  pos.push_back(GEOPosition.z()/km); 
  theSpenvisCSVFiles[0].AddMetaVariable("GEO_POSITION",pos,"%6.2f","km");
  pos.clear();
  
  G4ThreeVector GEODirection = thePrimaryGenerator->GetGEODirection();
  std::vector<double> dir ;
  dir.push_back(GEODirection.x());
  dir.push_back(GEODirection.y());
  dir.push_back(GEODirection.z());
  theSpenvisCSVFiles[0].AddMetaVariable("GEO_DIRECTION",dir,"%6.2f","");
  dir.clear();
  
  AddMetaVariableToCSVFile(G4String("GEODETIC_ALT"),
			   G4String("km"),
			   GEOIDalt/km,
			   G4String("%6.2f"),
			   0);
  AddMetaVariableToCSVFile(G4String("GEODETIC_LAT"),
			   G4String("degree"),
			   GEOIDlat/degree,
			   G4String("%6.2f"),
			   0);			   
  AddMetaVariableToCSVFile(G4String("GEODETIC_LON"),
			   G4String("degree"),
			   GEOIDlong/degree,
			   G4String("%6.2f"),
			   0);
  AddMetaVariableToCSVFile(G4String("GEODETIC_ZEN"),
			   G4String("degree"),
			   GEOIDzenith/degree,
			   G4String("%6.2f"),
			   0);			   
  AddMetaVariableToCSVFile(G4String("GEODETIC_AZIM"),
			   G4String("degree"),
			   GEOIDazimuth/degree,
			   G4String("%6.2f"),
			   0);
  
  //particle type
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
  //definition of variable
  
  theSpenvisCSVFiles[0].AddVariable("rigidity","GV",1,"rigidity","%6.3f");
  theSpenvisCSVFiles[0].AddVariable("trajectory label","",1,"trajectory label","%1.0f");
  theSpenvisCSVFiles[0].AddVariable("asymptotic latitude","degree",1,"asymptotic latitude","%4.2f");
  theSpenvisCSVFiles[0].AddVariable("asymptotic longitude","degree",1,"asymptotic longitude","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("last position in GEO coordinate ","Re",
  					3,"last position in GEO coordinate","%6.3f");
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::RegisterAsymptoticDirection(G4ThreeVector LastPosition,
                                    G4ThreeVector LastMomentaOnCharge,
			            G4int  FilterValue)
{ //Rigidity 
    
    G4double rigidity =LastMomentaOnCharge.mag()/1000.;
    
    //Asymptotic direction
    
    G4ThreeVector AsympDir= -1.*LastMomentaOnCharge;
    G4double Aslat=90.-(AsympDir.theta()/degree);
    G4double Aslong=AsympDir.phi()/degree;
    
    
    //LastPosition
    
    G4ThreeVector LastGEOPosition;
    LastGEOPosition = LastPosition/re;
    
    //give data
    std::vector<double > values;
    values.clear();
    values.push_back(rigidity);
    values.push_back(double(FilterValue));
    values.push_back(Aslat);
    values.push_back(Aslong);
    values.push_back(LastGEOPosition.x());
    values.push_back(LastGEOPosition.y());
    values.push_back(LastGEOPosition.z());
    WriteData(values);
    
   					    
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::InitialiseCutoffVsPosCSVFile(G4String coor_sys, 
   				     			G4double zenith,
				     			G4double azimuth)
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  //additional meta variable
  theSpenvisCSVFiles[0].AddMetaVariableStr("COORDINATE_SYSTEM",coor_sys);
  AddMetaVariableToCSVFile(G4String("DIRECTION_ZEN"),
			   G4String("degree"),
			   zenith/degree,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("DIRECTION_AZIM"),
			   G4String("degree"),
			   azimuth/degree,
			   G4String("%6.2f"),
			   0);
  			   
  
  
  //get a pointer to the primary generator action
  
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  //particle type
  
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
  
  
  //variables
  
  theSpenvisCSVFiles[0].AddVariable("SMJD","day",1,"Spenvis Modified julian date","%16.8f");
  theSpenvisCSVFiles[0].AddVariable("Altitude","km",1,"Altitude","%9.3f");
  theSpenvisCSVFiles[0].AddVariable("Latitude","degree",1,"Latitude","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Longitude","degree",1,"Longitude","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Ru","GV",1,"Upper cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rc","GV",1,"Effective cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rl","GV",1,"Lower cutoff rigidity","%.3e");
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::InitialiseCutoffVsDirCSVFile(G4String coor_sys,
							G4double altitude,
				     			G4double latitude,
				     			G4double longitude 
   				     			)
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  //additional meta variable
  theSpenvisCSVFiles[0].AddMetaVariableStr("COORDINATE_SYSTEM",coor_sys);
  AddMetaVariableToCSVFile(G4String("ALTITUDE"),
			   G4String("km"),
			   altitude/km,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("LATITUDE"),
			   G4String("degree"),
			   latitude/degree,
			   G4String("%6.2f"),
			   0);
   AddMetaVariableToCSVFile(G4String("LONGITUDE"),
			   G4String("degree"),
			   longitude/degree,
			   G4String("%6.2f"),
			   0);
  
  //get a pointer to the primary generator action
  
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  //particle type
  
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
   			   			   
  			   
  
  //variables
  
  theSpenvisCSVFiles[0].AddVariable("Zenith","degree",1,"Zenith angle of incoming direction","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Azimuth","degree",1,"Azimuth angle of incoming direction","%5.2f");
  theSpenvisCSVFiles[0].AddVariable("Ru","GV",1,"Upper cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rc","GV",1,"Effective cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rl","GV",1,"Lower cutoff rigidity","%.3e");
  
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::InitialiseCutoffVsTimeCSVFile(G4String coor_sys,
							G4double altitude,
				     			G4double latitude,
				     			G4double longitude,
							G4double zenith,
							G4double azimuth 
   							)
{ block_names.clear();
  block_names.push_back(G4String("00001"));
  theSpenvisCSVFiles.clear();
  theSpenvisCSVFiles.push_back(SpenvisCSV());
  WriteBfieldInformation();
  
  //additional meta variable
  theSpenvisCSVFiles[0].AddMetaVariableStr("COORDINATE_SYSTEM",coor_sys);
  AddMetaVariableToCSVFile(G4String("ALTITUDE"),
			   G4String("km"),
			   altitude/km,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("LATITUDE"),
			   G4String("degree"),
			   latitude/degree,
			   G4String("%6.2f"),
			   0);
  AddMetaVariableToCSVFile(G4String("LONGITUDE"),
			   G4String("degree"),
			   longitude/degree,
			   G4String("%6.2f"),
			   0);	
  
  AddMetaVariableToCSVFile(G4String("DIRECTION_ZEN"),
			   G4String("degree"),
			   zenith/degree,
			   G4String("%6.2f"),
			   0);
  
  AddMetaVariableToCSVFile(G4String("DIRECTION_AZIM"),
			   G4String("degree"),
			   azimuth/degree,
			   G4String("%6.2f"),
			   0);	
  
  			   		   			   
  			   		   
  //get a pointer to the primary generator action
  
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  //particle type
  
  G4String particle_name= thePrimaryGenerator
   				->GetParticleSource()
					->GetParticleDefinition()
						->GetParticleName();
  theSpenvisCSVFiles[0].AddMetaVariableStr("PARTICLE",particle_name);											                  
   			   			   
  			   
  
  //variables
  
  theSpenvisCSVFiles[0].AddVariable("SMJD","day",1,"Spenvis Modified julian date","%16.8f");
  theSpenvisCSVFiles[0].AddVariable("Ru","GV",1,"Upper cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rc","GV",1,"Effective cutoff rigidity","%.3e");
  theSpenvisCSVFiles[0].AddVariable("Rl","GV",1,"Lower cutoff rigidity","%.3e");
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::AddMetaVariableToCSVFile(G4String variable_name,
    				   G4String variable_unit,
				   G4double value, G4String format,
				    SpenvisCSV* aBlock)
{SpenvisCSV* theBlock;
 if (!aBlock) theBlock = &(theSpenvisCSVFiles[0]);
 else theBlock =aBlock;
 std::vector<double> values;
 values.push_back(value);
 theBlock->AddMetaVariable(variable_name,
 			   values, 
			   format, variable_unit);
 
}
/////////////////////////////////////////////////////////////////////////////
//
bool MAGCOSSpenvisManager::ReadSpenvisPositionGridCSVFile(G4String input_file,
   				   std::vector<double>& altitude,
   				   std::vector<double>& latitude,
				   std::vector<double>& longitude,
				   G4int nblock)
{ theSpenvisCSVCollection = new SpenvisCSVCollection();
  theSpenvisCSVCollection->ReadCollection(input_file);
  G4String block_name;
  if (nblock <= 0) block_name="00001";
  else {
  	std::stringstream astream;
	astream <<nblock;
	astream >>block_name;
	if (nblock  <10) block_name="0000"+block_name;
	else if  (nblock  <100) block_name="000"+block_name;
	else if  (nblock  <1000) block_name="00"+block_name;
	else if  (nblock  <10000) block_name="0"+block_name;
	
  }
  SpenvisCSV aSpenvisCSV = theSpenvisCSVCollection->GetCSVBlock(block_name);
  delete theSpenvisCSVCollection;
  theSpenvisCSVCollection = 0;
  
  G4int nb_var = aSpenvisCSV.GetNumVariables();
  G4int nb_lines = aSpenvisCSV.GetNumDataLines();
  altitude.clear();
  latitude.clear();
  longitude.clear();
  
  if (nb_var && nb_lines >0){
  	for (int i=0; i< nb_lines;i++){
		G4double alt,lon,lat;
		alt = aSpenvisCSV.GetDataRecord("Altitude",i)[0]*km;
		lat = aSpenvisCSV.GetDataRecord("Latitude",i)[0]*degree;
		lon = aSpenvisCSV.GetDataRecord("Longitude",i)[0]*degree;
		altitude.push_back(alt);
		latitude.push_back(lat);
		longitude.push_back(lon);		
	}
  	
  		
	return true;
		
  }
  return false;
  
   
  
   
 
}				   
////////////////////////////////////////////////////////////////////////////////
//
bool MAGCOSSpenvisManager::ReadSpenvisTrajectoryCSVFile(G4String input_file,
   				   std::vector<double>& altitude,
   				   std::vector<double>& latitude,
				   std::vector<double>& longitude,
				   std::vector<double>& smjd,
				   G4int nblock)
{ theSpenvisCSVCollection = new SpenvisCSVCollection();
  theSpenvisCSVCollection->ReadCollection(input_file);
  G4String block_name;
  if (nblock <= 0) block_name="00001";
  else {
  	std::stringstream astream;
	astream <<nblock;
	astream >>block_name;
	if (nblock  <10) block_name="0000"+block_name;
	else if  (nblock  <100) block_name="000"+block_name;
	else if  (nblock  <1000) block_name="00"+block_name;
	else if  (nblock  <10000) block_name="0"+block_name;
	
  }
  SpenvisCSV aSpenvisCSV = theSpenvisCSVCollection->GetCSVBlock(block_name);
  delete theSpenvisCSVCollection;
  theSpenvisCSVCollection = 0;
  
  G4int nb_var = aSpenvisCSV.GetNumVariables();
  G4int nb_lines = aSpenvisCSV.GetNumDataLines();
  altitude.clear();
  latitude.clear();
  longitude.clear(),
  smjd.clear();
  if (nb_var && nb_lines >0){
  	for (int i=0; i< nb_lines;i++){
		G4double alt,lon,lat,mjd;
		alt = aSpenvisCSV.GetDataRecord("Altitude",i)[0]*km;
		lat = aSpenvisCSV.GetDataRecord("Latitude",i)[0]*degree;
		lon = aSpenvisCSV.GetDataRecord("Longitude",i)[0]*degree;
		mjd = aSpenvisCSV.GetDataRecord("MJD",i)[0];
		altitude.push_back(alt);
		latitude.push_back(lat);
		longitude.push_back(lon);
		smjd.push_back(mjd);
		
	}
  	
  		
	return true;
		
  }
  return false;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::InitialiseTrajectoryCSVBlock(G4String Type)
{ if (RegisterParticleTrajectory){
	G4bool BlineCase = false;
	if (Type.contains("Magnetic")) BlineCase = true;
	theParticleTrajectoryCSVBlocks.push_back(SpenvisCSV());
	SpenvisCSV* theBlock = &(theParticleTrajectoryCSVBlocks.back());
	WriteBfieldInformation(theBlock);
	theBlock->AddMetaVariableStr("TRAJECTORY_TYPE",Type);
	
	//Start position and direction  
	//get a pointer to the primary generator action
  	
	MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
   		(MAGCOSPrimaryGeneratorAction*)
          		G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  	G4double GEOIDalt, GEOIDlat, GEOIDlong, GEOIDzenith, GEOIDazimuth; 
  	thePrimaryGenerator->GetGEOIDPosition( GEOIDalt,
  					       GEOIDlat,
					       GEOIDlong);
  					
  	thePrimaryGenerator->GetGEOIDDirection( GEOIDzenith, GEOIDazimuth);  
  	G4ThreeVector GEOPosition = thePrimaryGenerator->GetGEOPosition();
  	std::vector<double> pos ;
  	pos.push_back(GEOPosition.x()/re);
  	pos.push_back(GEOPosition.y()/re);
  	pos.push_back(GEOPosition.z()/re); 
  	theBlock->AddMetaVariable("GEO_STARTPOSITION",pos,"%6.2f","re");
  	pos.clear();
  
 	G4ThreeVector GEODirection = thePrimaryGenerator->GetGEODirection();
	std::vector<double> dir ;
  	if (!BlineCase) {
  		dir.push_back(GEODirection.x());
  		dir.push_back(GEODirection.y());
  		dir.push_back(GEODirection.z());
  		theBlock->AddMetaVariable("GEO_STARTDIRECTION",dir,"%6.2f","");
  		dir.clear();
	}	
	
	SpaceCoordinateConvertor* theConvertor = 
			SpaceCoordinateConvertor::getInstance();
	G4ThreeVector Position = theConvertor->Transform(GEOPosition,"GEO","GSM");	
	G4ThreeVector Direction = theConvertor->Transform(GEODirection,"GEO","GSM");
	pos.push_back(Position.x()/re);
  	pos.push_back(Position.y()/re);
  	pos.push_back(Position.z()/re); 
  	theBlock->AddMetaVariable("GSM_STARTPOSITION",pos,"%6.2f","re");
  	pos.clear();
	if (!BlineCase) {
  		dir.push_back(Direction.x());
  		dir.push_back(Direction.y());
  		dir.push_back(Direction.z());
  		theBlock->AddMetaVariable("GSM_STARTDIRECTION",dir,"%6.2f","");
  		dir.clear();
	}	
	
	Position = theConvertor->Transform(GEOPosition,"GEO","GSE");	
	Direction = theConvertor->Transform(GEODirection,"GEO","GSE");
	pos.push_back(Position.x()/re);
  	pos.push_back(Position.y()/re);
  	pos.push_back(Position.z()/re); 
  	theBlock->AddMetaVariable("GSE_STARTPOSITION",pos,"%6.2f","km");
  	pos.clear();
	if (!BlineCase) {
  		dir.push_back(Direction.x());
  		dir.push_back(Direction.y());
  		dir.push_back(Direction.z());
		theBlock->AddMetaVariable("GSE_STARTDIRECTION",dir,"%6.2f","");
		dir.clear();
	}	
  	
	
	Position = theConvertor->Transform(GEOPosition,"GEO","MAG");	
	Direction = theConvertor->Transform(GEODirection,"GEO","MAG");
	pos.push_back(Position.x()/re);
  	pos.push_back(Position.y()/re);
  	pos.push_back(Position.z()/re); 
  	theBlock->AddMetaVariable("MAG_STARTPOSITION",pos,"%6.2f","km");
  	pos.clear();
	if (!BlineCase) {
  		dir.push_back(Direction.x());
  		dir.push_back(Direction.y());
  		dir.push_back(Direction.z());
		theBlock->AddMetaVariable("MAG_STARTDIRECTION",dir,"%6.2f","");
  		dir.clear();
	}	
  	
	
	
			
  
  	AddMetaVariableToCSVFile(G4String("GEODETIC_ALT"),
			         G4String("km"),
			         GEOIDalt/km,
			         G4String("%6.2f"),
			         theBlock);
        AddMetaVariableToCSVFile(G4String("GEODETIC_LAT"),
			         G4String("degree"),
			         GEOIDlat/degree,
			         G4String("%6.2f"),
			         theBlock);			   
  	AddMetaVariableToCSVFile(G4String("GEODETIC_LON"),
			         G4String("degree"),
			         GEOIDlong/degree,
			         G4String("%6.2f"),
			         theBlock);
    	AddMetaVariableToCSVFile(G4String("GEODETIC_ZEN"),
			         G4String("degree"),
			         GEOIDzenith/degree,
			         G4String("%6.2f"),
			         theBlock);			   
  	AddMetaVariableToCSVFile(G4String("GEODETIC_AZIM"),
			         G4String("degree"),
			         GEOIDazimuth/degree,
			         G4String("%6.2f"),
			         theBlock);
	
	// particle type, energy, rigidity
  
  	if (!BlineCase) {
		G4ParticleDefinition* theParticleDefinition = 
					thePrimaryGenerator
   						->GetParticleSource()
							->GetParticleDefinition();
		G4String particle_name= theParticleDefinition->GetParticleName();
  		theBlock->AddMetaVariableStr("PARTICLE",particle_name);
		
		G4double energy = thePrimaryGenerator
   						->GetParticleSource()
							->GetParticleEnergy();
		G4cout<<energy<<std::endl;					
		AddMetaVariableToCSVFile(G4String("ENERGY"),
			         G4String("GeV"),
			         energy/GeV,
			         G4String("%.3e"),
			         theBlock);
		
		G4double m0 = theParticleDefinition->GetPDGMass();
		G4double E=energy+ m0;
		G4double rigidity =std::abs( std::sqrt( E*E - m0*m0)
				 		/theParticleDefinition->GetPDGCharge());		 
		
		AddMetaVariableToCSVFile(G4String("RIGIDITY"),
			         G4String("GV"),
			         rigidity/GV,
			         G4String("%.3e"),
			         theBlock);		 
				 
	}
	
	//definition of variable
  
  	theBlock->AddVariable("GEO_POS","Re",3,"position in GEO cartesian coordinate","%6.3f");
  	theBlock->AddVariable("MAG_POS","Re",3,"position in MAG cartesian coordinate","%6.3f");
  	theBlock->AddVariable("GSM_POS","Re",3,"position in GSM cartesian coordinate","%6.3f");
  	theBlock->AddVariable("GSE_POS","Re",3,"position in GSE cartesian coordinate","%6.3f");
       
  
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::RegisterTrajectoryPoint(G4ThreeVector GEOPosition)
{ if (RegisterParticleTrajectory){
	SpenvisCSV* theBlock = &(theParticleTrajectoryCSVBlocks.back());
	std::vector<double> pos;
	pos.clear();
	SpaceCoordinateConvertor* theConvertor = 
			SpaceCoordinateConvertor::getInstance();
	G4ThreeVector Position =GEOPosition/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
	
	Position = theConvertor->Transform(GEOPosition,"GEO","MAG")/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
	
	Position = theConvertor->Transform(GEOPosition,"GEO","GSM")/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
		
	Position = theConvertor->Transform(GEOPosition,"GEO","GSE")/re;
	pos.push_back(Position.x());
	pos.push_back(Position.y());
	pos.push_back(Position.z());
	
	WriteData(pos,theBlock);
	
		
  	
  }	
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSpenvisManager::SaveTrajectoryCSVBlocks(G4String name)
{ theSpenvisCSVCollection = new SpenvisCSVCollection();
  for (unsigned int i=0; i<theParticleTrajectoryCSVBlocks.size(); i++){
  		G4String block_name;
		std::stringstream astream;
		astream<<i+1;
		astream >>block_name;
		if (i+1  <10) block_name="0000"+block_name;
		else if  (i+1  <100) block_name="000"+block_name;
		else if  (i+1  <1000) block_name="00"+block_name;
		else if  (i+1  <10000) block_name="0"+block_name;
		theSpenvisCSVCollection->AddCSVBlock(block_name,
						  theParticleTrajectoryCSVBlocks[i]);
	
  }
  theSpenvisCSVCollection->OutputCollection(name);
  delete theSpenvisCSVCollection;
}
