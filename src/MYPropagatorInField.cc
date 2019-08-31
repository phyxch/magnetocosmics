// Updated on 9/16/2014, Hexc, Olesya: added CLHEP namespace for g4.10 version

#include "MYPropagatorInField.hh"
#include "G4ios.hh"
#include "iomanip"

using namespace CLHEP;

const G4double  MYPropagatorInField::fEpsilonMinDefault = 5.0e-7;  
const G4double  MYPropagatorInField::fEpsilonMaxDefault = 0.05;

MYPropagatorInField::MYPropagatorInField(G4Navigator    *theNavigator, 
		                         G4FieldManager *detectorFieldMgr)
  : fDetectorFieldMgr(detectorFieldMgr), 
    fCurrentFieldMgr(detectorFieldMgr), 
    fNavigator(theNavigator),
    End_PointAndTangent(G4ThreeVector(0.,0.,0.),
			G4ThreeVector(0.,0.,0.),0.0,0.0,0.0,0.0,0.0),
    fVerboseLevel(0),
    fEpsilonMin(fEpsilonMinDefault),
    fEpsilonMax(fEpsilonMaxDefault),  
    fmax_loop_count(10000),
    fNoZeroStep(0)
{
  // this->fChordFinder = new G4ChordFinder( (G4MagneticField*)0, 1e-6 );
  
  fActionThreshold_NoZeroSteps= 2; 
  fSevereActionThreshold_NoZeroSteps= 10; 
  fAbandonThreshold_NoZeroSteps= 50; 
  // fMidPoint_CurveLen_of_LastAttempt= -1;
  fFull_CurveLen_of_LastAttempt= -1; 
  fLast_ProposedStepLength= -1; 
  
  fLargestAcceptableStep= 1000.0 * meter;
  
  //Motion in G4
  fBSEquationInG4=new BSEquationInG4();
  fBSEquationInG4->SetNavigator(theNavigator);
     
}

///////////////////////////////////////////////////////////////////////////
//
// Compute the next geometric Step

G4double MYPropagatorInField::ComputeStep( G4FieldTrack& pFieldTrack,
					   G4double  CurrentProposedStepLength,
					   G4double&  currentSafety,                // IN/OUT
					   G4VPhysicalVolume* )
{ 
  // Set The Equation of motion for MotionInG4
  
  // should be change too complicated
  //------------------
  const G4MagIntegratorStepper* aStepper =fCurrentFieldMgr->GetChordFinder()
    ->GetIntegrationDriver()->GetStepper();
  
  G4MagIntegratorStepper* copyStepper = 
    const_cast<G4MagIntegratorStepper*> (aStepper);   
  G4EquationOfMotion* theEquation  = copyStepper->GetEquationOfMotion(); 
  fBSEquationInG4->SetEquationOfMotion(theEquation);
  fBSEquationInG4->SetCrossingDelta(
				    fCurrentFieldMgr->GetDeltaIntersection());	      
  
  //Compute the integration
  //---------------
  G4double StepLength= 
    fBSEquationInG4->ComputeTrajectory(pFieldTrack,
				       CurrentProposedStepLength,
				       currentSafety);
  
  return StepLength; 
  
}

//////////////////////////////////////////////////////////////
void MYPropagatorInField::SetChargeMomentumMass(G4ChargeState ChargeState,            // in e+ units
						G4double Momentum,          // in GeV/c 
						G4double Mass)              // in ? units
{ 
  if  (!fBSEquationInG4->GetEquationOfMotion())
    {const G4MagIntegratorStepper* aStepper =fCurrentFieldMgr->GetChordFinder()
	->GetIntegrationDriver()->GetStepper();
      
      G4MagIntegratorStepper* copyStepper = 
	const_cast<G4MagIntegratorStepper*> (aStepper);	  
      G4EquationOfMotion* theEquation  = copyStepper->GetEquationOfMotion(); 
      fBSEquationInG4->SetEquationOfMotion(theEquation);
    }
  fBSEquationInG4->GetEquationOfMotion()->SetChargeMomentumMass(ChargeState, Momentum, Mass); 
}
