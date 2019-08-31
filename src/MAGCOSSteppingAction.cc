#include "MAGCOSSteppingAction.hh"
#include "MAGCOSEventAction.hh"
#include "G4SteppingManager.hh"
#include "MAGCOSMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UserLimits.hh"
#include "SpaceCoordinateConvertor.hh"
#include "DurationManager.hh"
#include "G4EventManager.hh"
MAGCOSSteppingAction::MAGCOSSteppingAction()
{ stop_altitude = 19.99*km;}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  DurationManager* theDurationManager = DurationManager::getInstance();
  if (!theDurationManager->CheckDuration()) {
  	G4cout<< "The event will be aborted"<<std::endl;
	G4EventManager::GetEventManager()->AbortCurrentEvent();
  }
  
  
  
  
  const G4VPhysicalVolume* currentVolume = aStep->GetPostStepPoint()
                ->GetPhysicalVolume();
  G4Track* aTrack=aStep->GetTrack();
  G4double StepLength = aTrack->GetStepLength(); 	
  G4int nStep =aTrack->GetCurrentStepNumber();	
  if (nStep > 10 && StepLength == 0.) aTrack->SetTrackStatus(fStopAndKill); 
  if (currentVolume){		
   	const G4String name = currentVolume->GetName();
    	if (name =="Earth") { 
     		G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
      		G4double altitude,longitude,latitude;
    		//  G4cout<<"ok"<<endl;
      		SpaceCoordinateConvertor::getInstance()
              		->ComputeGEOIDCoordinatesFromGEOPosition
	              		(position,altitude,longitude,latitude);
    		// G4cout<<altitude/km<<std::endl;		      
    		//  G4cout<<"ok1"<<endl;
      		if (altitude <= stop_altitude )
		                       aTrack->SetTrackStatus(fStopAndKill);
     	}
    	else {  
     		G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
      
      		//check if the particle is outside the magnetosphere
      		MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
      		G4TransportationManager::GetTransportationManager()->
                              		GetFieldManager()->GetDetectorField();
      		if (theField->OutsideMagnetosphere(position)){
           		aTrack->SetTrackStatus(fStopAndKill);
	    		
			return;
		}
	}	
		
	G4UserLimits* theUserLimits = 
                  	currentVolume->GetLogicalVolume()->GetUserLimits();			 
      
      	//stop the particle if the track length is greater than the maximum
      	// allowed track length
      
      	G4double max_track_length =
                      		theUserLimits->GetUserMaxTrackLength(*aTrack);
      
      	if (aTrack->GetTrackLength() >= max_track_length){
           		aTrack->SetTrackStatus(fStopAndKill);
			
	  		return;
	} 
			 
			 
     	//stop the particle if the proper time of the particle is greater than the maximum allowed
     	// proper time of a trajectory
      
      	G4double max_local_time =
                     		theUserLimits->GetUserMaxTime(*aTrack);
      	if (aTrack->GetLocalTime() >= max_local_time){
           	aTrack->SetTrackStatus(fStopAndKill);
	   	G4cout<<"Stop  at local time ="<<aTrack->GetLocalTime()/s<< "second"<<std::endl;
	    	return;
	}   			 					    			 
     	 
  }
}
