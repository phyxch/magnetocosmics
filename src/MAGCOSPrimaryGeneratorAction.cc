#include "MAGCOSPrimaryGeneratorAction.hh"
#include "MAGCOSPrimaryMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "magneto_fsubroutine_def.hh"
#include "SpaceCoordinateConvertor.hh"
#include "MAGCOSUnits.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "MAGCOSMagneticField.hh"
#include "G4Proton.hh"
#include "MAGCOSSpenvisManager.hh"

MAGCOSPrimaryGeneratorAction::
      		MAGCOSPrimaryGeneratorAction()
{ myParticleSource= new G4GeneralParticleSource();
  myMessenger = new MAGCOSPrimaryMessenger(this);
  pGeneratePrimaries= 
         &MAGCOSPrimaryGeneratorAction::GenerateStandardPrimaries;
  SetDefaultRigidityVector();
  
  // first position and direction 
  G4bool test;  
  test = SetPositionAndDirection("GEOID", 20.*km, 7.98*degree, 46.55*degree, 0.,0.);
#ifdef USE_OLD_GPS  
  myParticleSource->SetParticleDefinition(G4Proton::Proton());
#else
  myParticleSource->SetParticleDefinition(G4Proton::Proton());
#endif  
  verbosity =0;
  InitialiseTrajectoryCSVBlockForForwardCase =false;
  InitialiseTrajectoryCSVBlockForBackwardCase =false;	
 
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSPrimaryGeneratorAction::~MAGCOSPrimaryGeneratorAction()
{ delete myParticleSource;}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ (this->*pGeneratePrimaries)(anEvent);
   if (InitialiseTrajectoryCSVBlockForForwardCase)
   		MAGCOSSpenvisManager::getInstance()->InitialiseTrajectoryCSVBlock("Forward particle trajectory");
   else if (InitialiseTrajectoryCSVBlockForBackwardCase)
   		MAGCOSSpenvisManager::getInstance()->InitialiseTrajectoryCSVBlock("Backward particle trajectory");		
   
  
   if (verbosity >0) 
          PrintPrimaryVertexInformation(anEvent->GetPrimaryVertex());
    	  
    	
}
///////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::GenerateStandardPrimaries(G4Event* anEvent)
{ // myParticleSource->Set("Mono");
  
  
  myParticleSource->GeneratePrimaryVertex(anEvent);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::
      GeneratePrimariesForComputingCutOffRigidity(G4Event* anEvent)
{ 
  Rigidity=rigidity_values[rigidity_index];
  G4double P = Rigidity* myParticleSource
                           ->GetParticleDefinition()->GetPDGCharge();
  G4double E0 = myParticleSource->GetParticleDefinition()->GetPDGMass()/GeV;
  G4double Etot=sqrt(E0*E0+P*P)*GeV;
#ifdef USE_OLD_GPS    
  myParticleSource->SetEnergyDisType("Mono");
  myParticleSource->SetMonoEnergy(Etot-E0*GeV);
#else
  myParticleSource->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
  myParticleSource->GetCurrentSource()->GetEneDist()->SetMonoEnergy(Etot-E0*GeV);	
#endif  
  rigidity_index++;
  myParticleSource->GeneratePrimaryVertex(anEvent);
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool MAGCOSPrimaryGeneratorAction
       ::SetPositionAndDirection(const G4String CoordSys, 
                                 const G4double anAlt,const G4double aLong,
				 const G4double aLat,const G4double aZenith,
				 const G4double  anAzimuth) 
{ G4bool test = SetPosition(CoordSys,anAlt,aLong,aLat);
  if (test)  test = SetDirection(CoordSys,aZenith,anAzimuth);
  return test;  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction
       ::SetPositionAndDirection(const G4String CoordSys, 
                                 const G4ThreeVector position,
				 const G4ThreeVector direction) 
{ G4bool test = SetPosition(CoordSys,position);
  if (test) test = SetDirection(CoordSys,direction);  
}
/////////////////////////////////////////////////////////////////////////////////
//
G4bool MAGCOSPrimaryGeneratorAction
       ::SetPosition(const G4String CoordSys, 
                     const G4double anAlt,const G4double aLong,
                     const G4double aLat)
{ if (CoordSys != "GSM" &&  CoordSys != "GEO" && CoordSys != "SM" &&  CoordSys != "MAG"
     && CoordSys != "GEOID" && CoordSys != "GEI" && CoordSys != "GSE"){
         G4cout<<CoordSys<<" is not a good system of coordinate"<<G4endl;
	 return false;
  }  
 
  SpaceCoordinateConvertor* theConvertor=
                            SpaceCoordinateConvertor::getInstance();
 
  if (CoordSys == "GEOID"){
  	GEOIDaltitude=anAlt;
   	GEOIDlongitude=aLong;
   	GEOIDlatitude=aLat;
   	theConvertor->ComputeGEOPositionFromGEOID(anAlt,aLat,aLong,GEOPosition);
  }
  else {
 	G4ThreeVector position=G4ThreeVector(1.,0.,0.);
  	position.setPhi(aLong);
  	position.setTheta(90.*degree-aLat);
  	position.setMag(re+anAlt);
  	GEOPosition=theConvertor->Transform(position,CoordSys,"GEO");
  	theConvertor->ComputeGEOIDCoordinatesFromGEOPosition
                     (GEOPosition,GEOIDaltitude, GEOIDlongitude, GEOIDlatitude);
  }
#ifdef USE_OLD_GPS  
  myParticleSource->SetCentreCoords(GEOPosition);
#else
  myParticleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(GEOPosition);	
#endif
  
  return true; 		     
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool MAGCOSPrimaryGeneratorAction::SetPosition(const G4String CoordSys, 
                                               const G4ThreeVector aPosition)
{ if (CoordSys != "GSM" &&  CoordSys != "GEO" && CoordSys != "SM" && 
           CoordSys != "GEI" && CoordSys != "GSE" && CoordSys != "MAG"){
   	G4cout<<CoordSys<<" is not a good system of coordinate"<<G4endl;
    	return false;
  }  
 
  SpaceCoordinateConvertor* theConvertor=
                            SpaceCoordinateConvertor::getInstance();
  GEOPosition=theConvertor->Transform(aPosition,CoordSys,"GEO");
  theConvertor->ComputeGEOIDCoordinatesFromGEOPositionAndDirection
                     (GEOPosition,GEODirection,
		      GEOIDaltitude, GEOIDlongitude, GEOIDlatitude,
		      GEOIDzenith,GEOIDazimuth);
#ifdef USE_OLD_GPS  
  myParticleSource->SetCentreCoords(GEOPosition);
#else
  myParticleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(GEOPosition);	
#endif  
  
  return true; 		     
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool MAGCOSPrimaryGeneratorAction::SetPositionOnDipoleMagneticShell
                     (const G4String Reference, const G4double L, const G4double
		     latitude, const G4double  longitude)
{
  if (Reference != "MAG" && Reference != "SM"){
         	G4cout<<Reference<<" is not a good reference"<<G4endl;
	  	return false;
  }
 
  G4String system;
  system=Reference;
 
  SpaceCoordinateConvertor* theConvertor=
                            SpaceCoordinateConvertor::getInstance();

  G4double coslat=std::cos(latitude);
  G4double r=L*re*coslat*coslat;
 
  G4ThreeVector Position=r*G4ThreeVector(1.,0.,0.);
  Position.setPhi(longitude); 
  Position.setTheta(90.*degree -latitude);
  

  GEOPosition = theConvertor->Transform(Position,system,"GEO");
			  

 // dipole Schift 

  G4ThreeVector GEODipoleSchift = ((MAGCOSMagneticField*)
  G4TransportationManager::GetTransportationManager()
      ->GetFieldManager()->GetDetectorField())->GetDipoleSchift();
 
  GEOPosition += GEODipoleSchift;
  
  theConvertor->ComputeGEOIDCoordinatesFromGEOPositionAndDirection
                     (GEOPosition,GEODirection,
		      GEOIDaltitude, GEOIDlongitude, GEOIDlatitude,
		      GEOIDzenith,GEOIDazimuth);

#ifdef USE_OLD_GPS  
  myParticleSource->SetCentreCoords(GEOPosition);
#else
  myParticleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(GEOPosition);	
#endif  
  return true; 
 
}
////////////////////////////////////////////////////////////////////////////////
//		     
G4bool MAGCOSPrimaryGeneratorAction
       ::SetDirection(const G4String CoordSys, 
                      const G4double aZenith,const G4double anAzimuth)
{if (CoordSys != "GSM" &&  CoordSys != "GEO" && CoordSys != "SM" &&  CoordSys != "MAG"
     && CoordSys != "GEOID" && CoordSys != "GEI" && CoordSys != "GSE")
         {G4cout<<CoordSys<<" is not a good system of coordinate"<<G4endl;
	  return false;}  
 
 SpaceCoordinateConvertor* theConvertor=
                            SpaceCoordinateConvertor::getInstance();
 
 if (CoordSys == "GEOID")
  {
   theConvertor->ComputeGEODirectionAndPositionFromGEOID
	            (GEOIDaltitude, GEOIDlatitude, GEOIDlongitude,
		     aZenith, anAzimuth,
		     GEOPosition, GEODirection);
  }
 else
 {
  G4ThreeVector position=theConvertor->Transform(GEOPosition,"GEO",CoordSys);
  G4ThreeVector direction =G4ThreeVector(0.,0.,1.);
  direction.setTheta(aZenith);
  direction.setPhi(180.*degree-anAzimuth);
  direction=-direction.rotateY(position.theta()).rotateZ(position.phi());
  GEODirection=theConvertor->Transform(direction,CoordSys,"GEO");
 }
 
 
 theConvertor->ComputeGEOIDCoordinatesFromGEOPositionAndDirection
                     (GEOPosition,GEODirection,
		      GEOIDaltitude, GEOIDlongitude, GEOIDlatitude,
		      GEOIDzenith,GEOIDazimuth);
		      
#ifdef USE_OLD_GPS  
  myParticleSource->SetParticleMomentumDirection(GEODirection);
#else
  myParticleSource->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(GEODirection);	
#endif
 
 return true; 		     
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool MAGCOSPrimaryGeneratorAction::SetDirection(const G4String CoordSys, 
                                               const G4ThreeVector aDirection)
{ if (CoordSys != "GSM" &&  CoordSys != "GEO" && CoordSys != "SM" && 
           CoordSys != "GEI" && CoordSys != "GSE" && CoordSys != "MAG"){
   	G4cout<<CoordSys<<" is not a good system of coordinate"<<G4endl;
    	return false;
  }  
 
  SpaceCoordinateConvertor* theConvertor=
                            SpaceCoordinateConvertor::getInstance();
  GEODirection=theConvertor->Transform(aDirection,CoordSys,"GEO");
  
  theConvertor->ComputeGEOIDCoordinatesFromGEOPositionAndDirection
                     (GEOPosition,GEODirection,
		      GEOIDaltitude, GEOIDlongitude, GEOIDlatitude,
		      GEOIDzenith,GEOIDazimuth);
#ifdef USE_OLD_GPS  
  myParticleSource->SetParticleMomentumDirection(GEODirection);
#else
  myParticleSource->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(GEODirection);	
#endif
  return true;
}
/////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::DefineDirectionByPitchAngle(G4double pitch_angle, G4double phi)
{
#ifdef USE_OLD_GPS  
  G4ThreeVector geo_pos = myParticleSource->GetCentreCoords();
#else
  G4ThreeVector geo_pos = myParticleSource->GetCurrentSource()
  					  ->GetPosDist()
  					  ->GetCentreCoords();	
#endif

 const MAGCOSMagneticField* theField =
                    reinterpret_cast<const MAGCOSMagneticField*>
	                 (G4TransportationManager::GetTransportationManager()->
 		                                GetFieldManager()->GetDetectorField()); 
  G4ThreeVector Bfield =theField->GetFieldValue(geo_pos);
  
  // construction of an orthogonal system with 
  // z axis parallel to Bfield
  // y axis in the direction  of vector zaxis X position
  // x axis completes the system
 
  G4ThreeVector z_axis = Bfield/Bfield.mag();
  G4ThreeVector y_axis = z_axis.cross(geo_pos);
  y_axis = y_axis/y_axis.mag();
  G4ThreeVector x_axis=y_axis.cross(z_axis);
  G4double cos_t =std::cos(pitch_angle);
  G4double sin_t =std::sin(pitch_angle);
 
 
  GEODirection = sin_t*std::cos(phi)*x_axis;
  GEODirection += sin_t*std::sin(phi)*y_axis;
  GEODirection += cos_t*z_axis;
  
  SpaceCoordinateConvertor* theConvertor=
                            SpaceCoordinateConvertor::getInstance();
  theConvertor->ComputeGEOIDCoordinatesFromGEOPositionAndDirection
                     (GEOPosition,GEODirection,
		      GEOIDaltitude, GEOIDlongitude, GEOIDlatitude,
		      GEOIDzenith,GEOIDazimuth);
 
#ifdef USE_OLD_GPS  
  myParticleSource->SetParticleMomentumDirection(GEODirection);
#else
  myParticleSource->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(GEODirection);	
#endif
  
}	
/////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::SetRigidity(G4double aRigidity)
{ G4double E0 = myParticleSource->GetParticleDefinition()->GetPDGMass();
  G4double P = aRigidity* myParticleSource
                           ->GetParticleDefinition()->GetPDGCharge();
  G4double Etot=sqrt(E0*E0+P*P);
#ifdef USE_OLD_GPS    
  myParticleSource->SetEnergyDisType("Mono");
  myParticleSource->SetMonoEnergy(Etot-E0);
#else
  myParticleSource->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
  myParticleSource->GetCurrentSource()->GetEneDist()->SetMonoEnergy(Etot-E0);	
#endif
 		      
}
/////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::AddValuesToRigidityVector
                                     (G4int nvalues,G4double val1, G4double step)
{ G4double val2=val1+G4double(nvalues)*step;
  if (val1 >=0 && val2 >=0){ 
   	G4double val_max= std::max(val1,val2);
    	G4double val_min= std::min(val1,val2);
    	step=std::abs(step);
    	if (rigidity_values.empty()){
     		for (G4int i=0; i<nvalues;i++)
               		rigidity_values.push_back(val_max - step * i);
     	}
    	else if (val_min > rigidity_values[0]){
     		for (G4int i=0; i<nvalues;i++)
               		rigidity_values.insert(rigidity_values.begin(),val_max - step * i);
     	}
    	else if  (val_max < rigidity_values.back()){
     		for (G4int i=0; i<nvalues;i++)
               		rigidity_values.push_back(val_max - step * i);
     	}
   	else { 
     		G4cout<<"error when adding values to rigidity vector"<<G4endl; 
     	}
  }
  else { 
  	G4cout<<"error when adding values to rigidity vector"<<G4endl;
  }   
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::SetDefaultRigidityVector()
{ rigidity_values.clear();
  AddValuesToRigidityVector(100,20,-0.1);
  AddValuesToRigidityVector(900,10,-0.01);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::SelectTypeOfPrimaries(G4String aString)
{ if (aString == "RigidityFilter") 
      pGeneratePrimaries= 
           &MAGCOSPrimaryGeneratorAction::GeneratePrimariesForComputingCutOffRigidity;
  else pGeneratePrimaries= 
           &MAGCOSPrimaryGeneratorAction::GenerateStandardPrimaries;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction
                 ::PrintPrimaryVertexInformation(const G4PrimaryVertex* aPrimaryVertex)
{G4ThreeVector geo_pos = aPrimaryVertex->GetPosition()/re;
 const G4PrimaryParticle* aPrimary = aPrimaryVertex->GetPrimary();
 
 G4ThreeVector geo_dir = aPrimary->GetMomentum();
 G4double P = geo_dir.mag();
 geo_dir/=P;
 G4double rigidity = std::abs(P/aPrimary->GetCharge()/GV);
 G4String particle_name =aPrimary->GetG4code()->GetParticleName();
 G4double E0=aPrimary->GetG4code()->GetPDGMass();
 G4double E= std::sqrt(E0*E0 + P*P);
 G4double Ekin= (E-E0)/GeV;
 
 G4cout.precision(5);
 //G4cout<<setiosflags(0x1000);
 G4cout<<"New primary "<<std::endl;
 G4cout<<"Particle, Energy, Rigidity : "<<particle_name<<", "
                                        <<Ekin<<", "   
                                        <<rigidity<<std::endl;
					
 G4cout<<"GEO position : X "<<geo_pos.x()<<", Y "
                            <<geo_pos.y()<<", Z "   
                            <<geo_pos.z()<<", theta "
			    <<geo_pos.theta()/degree<<", phi "
			    <<geo_pos.phi()/degree<<std::endl; 
 G4cout<<"GEO direction : X "<<geo_dir.x()<<", Y "
                             <<geo_dir.y()<<", Z "   
                             <<geo_dir.z()<<", theta "
			     <<geo_dir.theta()/degree<<", phi "
			     <<geo_dir.phi()/degree<<std::endl; 
 if (verbosity>1){
 	SpaceCoordinateConvertor* theConvertor=
                            SpaceCoordinateConvertor::getInstance();
  	G4ThreeVector gsm_pos = theConvertor->Transform(geo_pos,"GEO","GSM");
  	G4ThreeVector mag_pos = theConvertor->Transform(geo_pos,"GEO","MAG");
  	G4ThreeVector sm_pos = theConvertor->Transform(geo_pos,"GEO","SM");
	G4ThreeVector gei_pos = theConvertor->Transform(geo_pos,"GEO","GEI");
  	G4ThreeVector gse_pos = theConvertor->Transform(geo_pos,"GEO","GSE");
  
  	G4ThreeVector gsm_dir = theConvertor->Transform(geo_dir,"GEO","GSM");
  	G4ThreeVector mag_dir = theConvertor->Transform(geo_dir,"GEO","MAG");
  	G4ThreeVector sm_dir = theConvertor->Transform(geo_dir,"GEO","SM");
	G4ThreeVector gei_dir = theConvertor->Transform(geo_dir,"GEO","GEI");
  	G4ThreeVector gse_dir = theConvertor->Transform(geo_dir,"GEO","GSE");
  	
	G4double altitude,longitude,latitude;
  	theConvertor->ComputeGEOIDCoordinatesFromGEOPosition(geo_pos*re,
                                                       		altitude, 
						       		longitude, 
						       		latitude); 
						       
  	altitude/=km;
	
	G4cout<<"GSE position : X "<<gse_pos.x()<<", Y "
                            <<gse_pos.y()<<", Z "   
                            <<gse_pos.z()<<", theta "
			    <<gse_pos.theta()/degree<<", phi "
			    <<gse_pos.phi()/degree<<std::endl; 
        G4cout<<"GSE direction : X "<<gse_dir.x()<<", Y "
                             <<gse_dir.y()<<", Z "   
                             <<gse_dir.z()<<", theta "
			     <<gse_dir.theta()/degree<<", phi "
			     <<gse_dir.phi()/degree<<std::endl; 	  
        G4cout<<"GEI position : X "<<gei_pos.x()<<", Y "
                            <<gei_pos.y()<<", Z "   
                            <<gei_pos.z()<<", theta "
			    <<gei_pos.theta()/degree<<", phi "
			    <<gei_pos.phi()/degree<<std::endl; 
        G4cout<<"GEI direction : X "<<gei_dir.x()<<", Y "
                             <<gei_dir.y()<<", Z "   
                             <<gei_dir.z()<<", theta "
			     <<gei_dir.theta()/degree<<", phi "
			     <<gei_dir.phi()/degree<<std::endl; 
 	G4cout<<"GSM position : X "<<gsm_pos.x()<<", Y "
                            <<gsm_pos.y()<<", Z "   
                            <<gsm_pos.z()<<", theta "
			    <<gsm_pos.theta()/degree<<", phi "
			    <<gsm_pos.phi()/degree<<std::endl; 
        G4cout<<"GSM direction : X "<<gsm_dir.x()<<", Y "
                             <<gsm_dir.y()<<", Z "   
                             <<gsm_dir.z()<<", theta "
			     <<gsm_dir.theta()/degree<<", phi "
			     <<gsm_dir.phi()/degree<<std::endl; 
          
	G4cout<<"MAG position : X "<<mag_pos.x()<<", Y "
                            <<mag_pos.y()<<", Z "   
                            <<mag_pos.z()<<", theta "
			    <<mag_pos.theta()/degree<<", phi "
			    <<mag_pos.phi()/degree<<std::endl; 
        G4cout<<"MAG direction : X "<<mag_dir.x()<<", Y "
                             <<mag_dir.y()<<", Z "   
                             <<mag_dir.z()<<", theta "
			     <<mag_dir.theta()/degree<<", phi "
			     <<mag_dir.phi()/degree<<std::endl; 
	
	G4cout<<"SM position : X "<<sm_pos.x()<<", Y "
                            <<sm_pos.y()<<", Z "   
                            <<sm_pos.z()<<", theta "
			    <<sm_pos.theta()/degree<<", phi "
			    <<sm_pos.phi()/degree<<std::endl; 
        G4cout<<"SM direction : X "<<sm_dir.x()<<", Y "
                             <<sm_dir.y()<<", Z "   
                             <<sm_dir.z()<<", theta "
			     <<sm_dir.theta()/degree<<", phi "
			     <<sm_dir.phi()/degree<<std::endl; 
	
	G4cout<<"GEOID position alt, lat, long : "<<altitude<<", "
                           <<latitude/degree<<", "   
                           <<longitude/degree<<std::endl;		   		   
			   		   
	
			   
	}		   
       // G4cout<<setiosflags(0x800);                     		    
 
 }
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSPrimaryGeneratorAction::PrintBfieldAtPrimary()
{ 
#ifdef USE_OLD_GPS  
  G4ThreeVector geo_pos = myParticleSource->GetCentreCoords();
#else
  G4ThreeVector geo_pos = myParticleSource->GetCurrentSource()
  					  ->GetPosDist()
  					  ->GetCentreCoords();	
#endif
  
  const MAGCOSMagneticField* theField =
                    reinterpret_cast<const MAGCOSMagneticField*>
	                 (G4TransportationManager::GetTransportationManager()->
		                                GetFieldManager()->GetDetectorField()); 
  theField->PrintBfield(geo_pos); 	
}
 
   			   			 


