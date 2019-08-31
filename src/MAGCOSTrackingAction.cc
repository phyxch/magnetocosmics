#include "MAGCOSTrackingAction.hh"
#include "G4ThreeVector.hh"
#include "MAGCOSApplicationScenario.hh"
#include "MAGCOSMagneticField.hh"
#include"G4TransportationManager.hh"
#include"G4FieldManager.hh"
#include"G4RunManager.hh"
#include"G4UserLimits.hh"
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSTrackingAction::MAGCOSTrackingAction ()
{RegisterLastPoint=false;}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSTrackingAction::~MAGCOSTrackingAction ()
{;}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSTrackingAction::PreUserTrackingAction (const G4Track* )
{
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSTrackingAction::PostUserTrackingAction (const G4Track* aTrack)
{if (RegisterLastPoint)
  {
   G4ThreeVector LastPosition = aTrack->GetPosition();
   G4double aCharge =std::abs(aTrack->GetDynamicParticle()
                                          ->GetDefinition()->GetPDGCharge());
   G4ThreeVector LastMomentumOnCharge = aTrack->GetMomentum() / aCharge;
   
   //check if particle come from outside the magnetosphere
   
   MAGCOSMagneticField* theField = (MAGCOSMagneticField*)
     G4TransportationManager::GetTransportationManager()->
                              GetFieldManager()->GetDetectorField();
   G4double FilterValue =0;
   if (theField->OutsideMagnetosphere(LastPosition)) FilterValue=1;
   else {
   	G4UserLimits* theUserLimits = 
                  	aTrack->GetLogicalVolumeAtVertex()->GetUserLimits();
   	G4double max_track_length =
                      		theUserLimits->GetUserMaxTrackLength(*aTrack);
   	if (aTrack->GetTrackLength()>=max_track_length) FilterValue=-1;		
   }
   
  // send information to MAGCOSApplicationScenario
   
   MAGCOSApplicationScenario* theApplicationScenario =
      (MAGCOSApplicationScenario*) G4RunManager::GetRunManager()->GetUserRunAction();
   theApplicationScenario
             ->RegisterTrackLastPoint(LastPosition,LastMomentumOnCharge,FilterValue);
  } 
}
////////////////////////////////////////////////////////////////////////////////
