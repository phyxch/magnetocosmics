// 9/18/2014: Hexc, Olesya - Add G4ChargeState.hh header in order to use SetChargeMomentumMass
//    in G4EquationOfMotion

#ifndef MYPropagatorInField_hh 
#define MYPropagatorInField_hh  1

#include "globals.hh"
#include "G4FieldTrack.hh"
// #include "G4VPhysicalVolume.hh"
// class  G4VPhysicalVolume;

#include "G4Navigator.hh"
#include "G4ChordFinder.hh"

#include "G4FieldManager.hh"
#include "BSEquationInG4.hh"

#include "G4ChargeState.hh"

// #include "G4MagIntegratorDriver.hh"

class MYPropagatorInField
{
  
public:  // with description
  
  MYPropagatorInField( G4Navigator    *theNavigator, 
		       G4FieldManager *detectorFieldMgr);
  
  MYPropagatorInField( G4Navigator   *theNavigator );
  ~MYPropagatorInField(){};
  
  G4double ComputeStep(G4FieldTrack  &pFieldTrack,
		       G4double pCurrentProposedStepLength,
		       G4double       &pNewSafety, 
		       G4VPhysicalVolume *pPhysVol=0 );
  
  void SetChargeMomentumMass( G4ChargeState ChargeState,  // in e+ units
			      G4double Momentum,          // in GeV/c 
			      G4double Mass); 
  
  inline G4ThreeVector  EndPosition() const;       
  inline G4ThreeVector  EndMomentumDir() const;
  inline G4bool         IsParticleLooping() const;

  // Return the state after the Step
  
  
  inline G4int  SetVerboseLevel( G4int Verbose );
  inline G4int  Verbose() const;
  
  inline G4double  GetDeltaIntersection() const;
  // Accuracy for boundary intersection.
  inline G4double  GetDeltaOneStep() const;
  // Accuracy for one tracking/physics step.
  
  inline void    SetAccuraciesWithDeltaOneStep(G4double deltaOneStep);  
  // Sets both accuracies for the Global (Detector) field, 
  // maintaining a particular ratio for accuracties 
  // of volume Intersection and Integration (in One Step).
  
  inline void    SetDeltaIntersection(G4double deltaIntersection);
  // Set accuracy of  intersection of a volume.  (only)
  inline void    SetDeltaOneStep(G4double deltaOneStep);  
  // Set accuracy for integration of one step.   (only)
  
  inline G4int   GetMaxLoopCount() const;
  inline void    SetMaxLoopCount(G4int new_max);
  // A maximum for the number of steps that a (looping) particle can take.
  
  void printStatus( 
		   const G4FieldTrack&  StartFT,
		   const G4FieldTrack&  CurrentFT, 
		   G4double             requestStep, 
		   G4double             safety,
		   G4int                Step, 
		   G4VPhysicalVolume*   startVolume);
  // Print Method - useful mostly for debugging.
  
  inline G4FieldTrack GetEndState() const;
  
  // Minimum for Relative accuracy of any Step 
  inline G4double  GetMinimumEpsilonStep() const;
  inline void      SetMinimumEpsilonStep(G4double newEpsMin);
  
  inline void      SetLargestAcceptableStep(G4double newBigDist);
  inline G4double  GetLargestAcceptableStep();
  
public:  // without description
  
  // void  SetGlobalFieldMgr( G4FieldManager *detectorFieldMgr );
  // The Field Manager of the Detector.
  
  inline G4FieldManager*  GetCurrentFieldManager();
  
public:  // no description
  
  inline void SetThresholdNoZeroStep(G4int noAct, G4int noHarsh, G4int noAbandon);
  inline G4int GetThresholdNoZeroSteps(G4int i); 
  
protected:
  
  // ----------------------------------------------------------------------
  //  DATA Members
  // ----------------------------------------------------------------------
  
private:
  
  G4FieldManager *fDetectorFieldMgr; 
  // The  Field Manager of the whole Detector.  (default)
  
  G4FieldManager *fCurrentFieldMgr;
  // The  Field Manager of the current volume (may be the one above.)
  
  G4Navigator   *fNavigator;
  
  //  STATE information
  // ------------------
  
  G4double    fEpsilonStep;
  // Relative accuracy for current Step (Calc.)
  
  G4FieldTrack    End_PointAndTangent;
  // End point storage
  
  G4bool      fParticleIsLooping;
  
  G4int  fVerboseLevel;
  // For debuging purposes
  
  //  Values for the small possible relative accuracy of a step
  //       (corresponding to the greatest possible integration accuracy)
  
  // Limits for the Relative accuracy of any Step 
  G4double  fEpsilonMin; 
  G4double  fEpsilonMax;
  static const G4double  fEpsilonMinDefault;         // Can be 1.0e-5 to 1.0e-10 ...
  static const G4double  fEpsilonMaxDefault;         // Can be 1.0e-1 to 1.0e-3 ...
  
  G4int  fmax_loop_count;
  
  //  Variables to keep track of "abnormal" case - which causes loop
  //
  G4int     fNoZeroStep;                //  Counter of zeroStep
  G4int     fActionThreshold_NoZeroSteps;       //  Threshold: above this - act
  G4int     fSevereActionThreshold_NoZeroSteps; //  Threshold to act harshly
  G4int     fAbandonThreshold_NoZeroSteps;      //  Threshold to abandon
  // G4double  fMidPoint_CurveLen_of_LastAttempt= -1;
  G4double  fFull_CurveLen_of_LastAttempt; 
  G4double  fLast_ProposedStepLength; 
  G4double  fLargestAcceptableStep;
  
  //  in G4
  BSEquationInG4*   fBSEquationInG4;
  
  
public:
    
    inline G4double GetEpsilonStep() const
   { 
      return fEpsilonStep; 
   }

  inline void SetEpsilonStep(G4double newEps)
   {
      fEpsilonStep=newEps;
   }   
   
};

//  Defines the constructor.
//
#include "MYPropagatorInField.icc"

#endif 


