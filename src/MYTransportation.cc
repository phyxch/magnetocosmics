// Updated on 8/6/2014: Hexc, Olesya
// Modified for using Geant4.9.6

//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: MYTransportation.cc,v 1.3 2014/09/18 14:15:23 hexc Exp $
// GEANT4 tag $Name:  $
// 
// ------------------------------------------------------------
//	GEANT 4  include file implementation
//
// ------------------------------------------------------------
//
// This class is a process responsible for the transportation of 
// a particle, ie the geometrical propagation that encounters the 
// geometrical sub-volumes of the detectors.
//
// It is also tasked with part of updating the "safety".
//
// =======================================================================
// Modified:   
//            29 June 2001, J. Apostolakis, D. Cote-Ahern, P. Gumplinger: 
//                          correction for spin tracking   
//            20 Febr 2001, J. Apostolakis:  update for new FieldTrack
//            22 Sept 2000, V. Grichine:     update of Kinetic Energy
//             9 June 1999, J. Apostolakis & S.Giani: protect full relocation
//                               used in DEBUG for track that started on surface
//                               and went step < tolerance
//    				Also forced fast relocation in all DEBUG cases
//    				 & changed #if to use DEBUG instead of VERBOSE
// Created:  19 March 1997, J. Apostolakis
// =======================================================================

#include "MYTransportation.hh"

///////////////////////////////////////////////////////////////////////////////
//
// Constructor

MYTransportation::MYTransportation() : G4VProcess(G4String("MYTransportation") )
{
  G4TransportationManager* transportMgr ; 
  
  transportMgr = G4TransportationManager::GetTransportationManager() ; 
  
  fLinearNavigator =   transportMgr->GetNavigatorForTracking() ; 
  fFieldPropagator =   0 ; 
  
  // fFieldExists= false ; 
  
  fParticleIsLooping = false ; 
  
  // fGlobalFieldMgr=    transportMgr->GetFieldManager() ;  
  // fFieldPropagator=   transportMgr->GetPropagatorInField() ;
  fFieldPropagator=new
    MYPropagatorInField(fLinearNavigator,transportMgr->GetFieldManager());
  // Find out if an electromagnetic field exists
  // 
  // fFieldExists= transportMgr->GetFieldManager()->DoesFieldExist() ;
  // 
  //   The above code is problematic, because it only works if
  // the field manager has informed about the detector's field 
  // before this transportation process is constructed.
  // I cannot foresee how the transportation can be informed later. JA 
  //   The current answer is to ignore this data member and use 
  // the member function DoesGlobalFieldExist() in its place ...
  //    John Apostolakis, July 7, 1997
  
  fCurrentTouchableHandle = new G4TouchableHistory();
  
  // Initial value for safety and point-of-origin of safety
  
  fPreviousSafety    = 0.0 ; 
  fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;
  
  fEndGlobalTimeComputed= false;
  fCandidateEndGlobalTime= 0;
  // We added this line based G4EnclosingCylinder.cc
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

/////////////////////////////////////////////////////////////////////////////

MYTransportation::~MYTransportation()
{
}

//////////////////////////////////////////////////////////////////////////////
//
// Responsibilities:
//    Find whether the geometry limits the Step, and to what length
//    Calculate the new value of the safety and return it.
//    Store the final time, position and momentum.

G4double MYTransportation::
AlongStepGetPhysicalInteractionLength(  const G4Track&  track,
		                              G4double  ,
		                              G4double  currentMinimumStep,
		                              G4double& currentSafety,
		                              G4GPILSelection* selection  )
{
  G4double geometryStepLength, newSafety ; 
  fParticleIsLooping = false ;

  if(track.GetCurrentStepNumber()==1) {
     // reset safety value 
     fPreviousSafety    = 0.0 ; 
     fPreviousSftOrigin = G4ThreeVector(0.,0.,0.) ;
     // ChordFinder reset internal state
     //removed by Laurent
    /* if ( this->DoesGlobalFieldExist() )
        fFieldPropagator->GetChordFinder()->ResetStepEstimate();*/
     // We need to update the current transportation's touchable handle
     // to the track's one 
     fCurrentTouchableHandle = track.GetTouchableHandle();
  }

  // GPILSelection is set to defaule value of CandidateForSelection
  // It is a return value

  *selection = CandidateForSelection ;

  // Get initial Energy/Momentum of the track
  
  const G4DynamicParticle*  pParticle  = track.GetDynamicParticle() ;
  G4ThreeVector startMomentumDir       = pParticle->GetMomentumDirection() ;
  G4ThreeVector startPosition          = track.GetPosition() ;

  // G4double   theTime        = track.GetGlobalTime() ;

  // The Step Point safety is now generalised to mean the limit of assumption
  // of all processes, so it is not the previous Step's geometrical safety.
  // We calculate the starting point's safety here.

  G4ThreeVector OriginShift = startPosition - fPreviousSftOrigin ;
  G4double      MagSqShift  = OriginShift.mag2() ;
  if( MagSqShift >= fPreviousSafety*fPreviousSafety )
  {
     currentSafety = 0.0 ;
  }
  else
  {
     currentSafety = fPreviousSafety - sqrt(MagSqShift);
  }
  // Is the particle charged ?

  G4ParticleDefinition* pParticleDef   = pParticle->GetDefinition();
  G4double              particleCharge = pParticleDef->GetPDGCharge(); 
  G4ChargeState         particleChargeState(particleCharge,
					    pParticleDef->GetPDGMagneticMoment(),
					    pParticleDef->GetPDGSpin(),
					    0, 0);

  G4bool   fieldExertsForce = false ;
  fGeometryLimitedStep = false ;
  // fEndGlobalTimeComputed     = false ;

  // There is no need to locate the current volume. It is Done elsewhere:
  //   On track construction 
  //   By the tracking, after all AlongStepDoIts, in "Relocation"
  //  Does the particle have an (EM) field force exerting upon it?
  
  if( (particleCharge != 0.0) )
  {
     
     fieldExertsForce= this->DoesGlobalFieldExist() ;

     // Future: will/can also check whether current volume's field is Zero or
     //  set by the user (in the logical volume) to be zero.
  }
  //  Choose the calculation of the transportation: Field or not 
 
  if( !fieldExertsForce ) 
  {
     G4double linearStepLength ;
     
     if( currentMinimumStep <= currentSafety )
     {
       // The Step is guaranteed to be taken

       geometryStepLength   = currentMinimumStep ;
       fGeometryLimitedStep = false ;
     }
     else
     {
       //  Find whether the straight path intersects a volume

       linearStepLength = fLinearNavigator->ComputeStep( startPosition, 
                                                         startMomentumDir,
                                                         currentMinimumStep, 
                                                         newSafety) ;
       // Remember last safety origin & value.

       fPreviousSftOrigin = startPosition ;
       fPreviousSafety    = newSafety ; 

       // The safety at the initial point has been re-calculated:

       currentSafety = newSafety ;
			    
       if( linearStepLength <= currentMinimumStep)
       {
	       // The geometry limits the Step size (an intersection was found.)

         geometryStepLength   = linearStepLength ;
         fGeometryLimitedStep = true ;
       }
       else
       {
	       // The full Step is taken.

         geometryStepLength   = currentMinimumStep ;
         fGeometryLimitedStep = false ;
       }
     }
     endpointDistance = geometryStepLength ;

     // Calculate final position

     fTransportEndPosition = startPosition + geometryStepLength*startMomentumDir ;

     // Momentum direction, energy and polarisation are unchanged by transport

     fTransportEndMomentumDir   = startMomentumDir ; 
     fTransportEndKineticEnergy = track.GetKineticEnergy() ;
     fTransportEndSpin          = track.GetPolarization();
     fParticleIsLooping         = false ;
     fMomentumChanged           = false ; 
     fEndGlobalTimeComputed     = false ;
  }
  else
  {
     G4double       momentumMagnitude = pParticle->GetTotalMomentum() ;
     G4ThreeVector  EndUnitMomentum ;
     G4double       lengthAlongCurve ;
     G4double       restMass = pParticleDef->GetPDGMass() ;
 
     fFieldPropagator->SetChargeMomentumMass( particleChargeState,    // charge in e+ units
                                              momentumMagnitude, // Momentum in Mev/c 
                                              restMass           ) ;  
     
     G4ThreeVector spin           = track.GetPolarization() ;
     G4FieldTrack  aFieldTrack =
       G4FieldTrack( startPosition, 
		     track.GetMomentumDirection(),
		     0.0, 
		     track.GetKineticEnergy(),
		     restMass,
		     track.GetVelocity(),
		     track.GetGlobalTime(),   // Laboratory fr.
		     track.GetProperTime(),   // Particle rest fr.
		     &spin                   ) ;
     
     if( currentMinimumStep > 0 ) 
       {
	 //  Do the Transport in the field (non recti-linear)
	 
	 lengthAlongCurve = fFieldPropagator->ComputeStep( aFieldTrack,
							   currentMinimumStep, 
							   currentSafety,
							   track.GetVolume()) ;
	 if( lengthAlongCurve < currentMinimumStep)
	   {
	         geometryStepLength   = lengthAlongCurve ;
	         fGeometryLimitedStep = true ;
        }
        else
        {
  	       geometryStepLength   = currentMinimumStep ;
	         fGeometryLimitedStep = false ;
	      }
     }
     else
     {
        geometryStepLength   = lengthAlongCurve= 0.0 ;
	      fGeometryLimitedStep = false ;
     }
     // Remember last safety origin & value.

     fPreviousSftOrigin = startPosition ;
     fPreviousSafety    = currentSafety ; 		    
        
     // Get the End-Position and End-Momentum (Dir-ection)

     fTransportEndPosition = aFieldTrack.GetPosition() ;

     // Momentum:  Magnitude and direction can be changed too now ...

     fMomentumChanged         = true ; 
     fTransportEndMomentumDir = aFieldTrack.GetMomentumDir() ;

     fTransportEndKineticEnergy  = aFieldTrack.GetKineticEnergy() ; 

     // if( (track.GetKineticEnergy() - fTransportEndKineticEnergy) 
     //      > perMillion * fTransportEndKineticEnergy             ){

     if( fFieldPropagator->GetCurrentFieldManager()->DoesFieldChangeEnergy() ){

        // If the field can changed energy, then the time must be integrated
        //    - so this should have been updated
        fCandidateEndGlobalTime   = aFieldTrack.GetLabTimeOfFlight();
	fEndGlobalTimeComputed    = true;
	  // was ( fCandidateEndGlobalTime != track.GetGlobalTime() );
     // a cleaner way is to have FieldTrack knowing whether time is updated.
     }else{
        //  The energy is unchanged by field transport,
        //    so the time changed will be calculated elsewhere
        fEndGlobalTimeComputed    = false;
     }

     fTransportEndSpin = aFieldTrack.GetSpin();

     fParticleIsLooping = fFieldPropagator->IsParticleLooping() ;
     endpointDistance   = (fTransportEndPosition - startPosition).mag() ;
  }
  // If we are asked to go a step length of 0, and we are on a boundary
  //  then a boundary will also limit the step -> we must flag this.

  if (currentMinimumStep == 0.0 ) 
  {
     if( currentSafety == 0.0 )	fGeometryLimitedStep = true ;
  }

  // Update the safety starting from the end-point, if it will become 
  //  negative at the end-point.
  
  if( currentSafety < endpointDistance ) 
  {
      G4double endSafety = fLinearNavigator->ComputeSafety( fTransportEndPosition) ;
      currentSafety      = endSafety ;
      fPreviousSftOrigin = fTransportEndPosition ;
      fPreviousSafety    = currentSafety ; 

      // Because the Stepping Manager assumes it is from the start point, 
      //  add the StepLength

      currentSafety     += endpointDistance ;

#ifdef G4DEBUG_TRANSPORT 
      G4cout.precision(16) ;
      G4cout << "***Transportation::AlongStepGPIL ** " << G4endl  ;
      G4cout << "  Called Navigator->ComputeSafety at " << fTransportEndPosition
	     << "    and it returned safety= " << endSafety << G4endl ; 
      G4cout << "  Adding endpoint distance " << endpointDistance 
	     << "   we obtain pseudo-safety= " << currentSafety << G4endl ; 
#endif
  }				    
#ifdef BEFORE_V7 
  fParticleChange.SetTrueStepLength(geometryStepLength)  ;
#else
  fParticleChange.ProposeTrueStepLength(geometryStepLength);	
#endif


  return geometryStepLength ;
}

/////////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)

G4VParticleChange* MYTransportation::AlongStepDoIt( const G4Track& track,
                                                    const G4Step&  stepData )
{
  fParticleChange.Initialize(track) ;

  //  Code for specific process 
#ifdef BEFORE_V7   
  fParticleChange.SetPositionChange(fTransportEndPosition) ;
  fParticleChange.SetMomentumChange(fTransportEndMomentumDir) ;
  fParticleChange.SetEnergyChange(fTransportEndKineticEnergy) ;
  fParticleChange.SetMomentumChanged(fMomentumChanged) ;
  fParticleChange.SetPolarizationChange(fTransportEndSpin);
#else
  fParticleChange.ProposePosition(fTransportEndPosition) ;
  fParticleChange.ProposeMomentumDirection(fTransportEndMomentumDir) ;
  fParticleChange.ProposeEnergy(fTransportEndKineticEnergy) ;
  fParticleChange.SetMomentumChanged(fMomentumChanged) ;
  fParticleChange.ProposePolarization(fTransportEndSpin);
#endif
  
  G4double deltaTime = 0.0 ;

  // Calculate  Lab Time of Flight (ONLY if field Equations used it!)
     // G4double endTime   = fCandidateEndGlobalTime;
     // G4double delta_time = endTime - startTime;
  G4double startTime = track.GetGlobalTime();
  
  if (!fEndGlobalTimeComputed){
  //  The time was not integrated .. make the best estimate possible
     G4double finalVelocity   = track.GetVelocity() ;
     G4double initialVelocity = stepData.GetPreStepPoint()->GetVelocity();
     G4double stepLength = track.GetStepLength();

     if (finalVelocity > 0.0) { 
        G4double meanInverseVelocity ;
        // deltaTime = stepLength/finalVelocity ;  
        meanInverseVelocity= 0.5 * ( 1.0 / initialVelocity + 1.0 / finalVelocity );
        deltaTime = stepLength * meanInverseVelocity ; 
     }else{
        deltaTime = stepLength/initialVelocity;     
     }
     fCandidateEndGlobalTime   = startTime + deltaTime; 
  }else
     deltaTime = fCandidateEndGlobalTime - startTime;
#ifdef BEFORE_V7     
  fParticleChange. SetTimeChange( fCandidateEndGlobalTime ) ;
#else
  fParticleChange.ProposeGlobalTime( fCandidateEndGlobalTime ) ;
#endif
  // Now Correct by Lorentz factor to get "proper" deltaTime
  
  G4double  restMass       = track.GetDynamicParticle()->GetMass() ;
 
  G4double deltaProperTime = deltaTime*( restMass/track.GetTotalEnergy() ) ;
#ifdef BEFORE_V7
  fParticleChange.SetProperTimeChange( track.GetProperTime() + deltaProperTime ) ;
#else
  fParticleChange.ProposeProperTime( track.GetProperTime() + deltaProperTime ) ;
#endif
  //fParticleChange. SetTrueStepLength( track.GetStepLength() ) ;

  // If the particle is caught looping or is stuck (in very difficult boundaries)
  //   in a magnetic field (doing many steps) 
  //   THEN this kills it ...
  if ( fParticleIsLooping )
  {
      // Kill the looping particle 
#ifdef BEFORE_V7 
      fParticleChange.SetStatusChange( fStopAndKill )  ;
#else
      fParticleChange.ProposeTrackStatus( fStopAndKill )  ;
#endif      
#ifdef G4VERBOSE
      G4cout << " MYTransportation is killing track that is looping or stuck " << G4endl
	     << "   This track has " << track.GetKineticEnergy() << " MeV energy."
	     << G4endl;
#endif
      // ClearNumberOfInteractionLengthLeft() ;
  }
  // Another (sometimes better way) is to use a user-limit maximum Step size
  //  to alleviate this problem .. 

  return &fParticleChange ;

}

////////////////////////////////////////////////////////////////////////////////
//
// This ensures that the PostStep action is always called,
//   so that it can do the relocation if it is needed.
// 

G4double MYTransportation::
PostStepGetPhysicalInteractionLength( const G4Track&,
                                            G4double   ,
                                            G4ForceCondition* pForceCond )
{ 
  *pForceCond = Forced ; 

  return DBL_MAX ;  // was kInfinity ; but convention now is DBL_MAX
}

/////////////////////////////////////////////////////////////////////////////
//

G4VParticleChange* MYTransportation::PostStepDoIt( const G4Track& track,
                                                   const G4Step&  stepData )
{
  G4TouchableHandle retCurrentTouchable ;   // The one to return

  //   Initialize ParticleChange  (by setting all its members equal
  //                               to corresponding members in G4Track)
  // fParticleChange.Initialize(track) ;  // To initialise TouchableChange
#ifdef BEFORE_V7
  fParticleChange.SetStatusChange(track.GetTrackStatus()) ;
#else
  fParticleChange.ProposeTrackStatus(track.GetTrackStatus()) ;
#endif
  // If the Step was determined by the volume boundary,
  // logically relocate the particle
  
  if( fGeometryLimitedStep )
  {  
    // fCurrentTouchable will now become the previous touchable, 
    //  and what was the previous will be freed.
    // (Needed because the preStepPoint can point to the previous touchable)
    //G4cout<<"that is the case"<<std::endl;
    fLinearNavigator->SetGeometricallyLimitedStep() ;
    fLinearNavigator->
    LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(),
                                               track.GetMomentumDirection(),
                                               fCurrentTouchableHandle,
                                               true                      ) ;
    
    // Check whether the particle is out of the world volume 
    //   If so it has exited and must be killed.
    
    if( fCurrentTouchableHandle->GetVolume() == 0 )
    {
#ifdef BEFORE_V7    
       fParticleChange.SetStatusChange( fStopAndKill )  ;
#else
       fParticleChange.ProposeTrackStatus( fStopAndKill )  ;
#endif
    }
    retCurrentTouchable = fCurrentTouchableHandle;
    fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;
  }
  else
  {                    // fGeometryLimitedStep  is false
//#ifdef G4DEBUG
#ifdef G4VERBOSE
    // Although the location is changed, we know that the physical 
    //   volume remains constant. 
    // In order to help in checking the user geometry
    //    we perform a full-relocation and check its result 
    //     *except* if we have made a very small step from a boundary
    //      (ie remaining inside the tolerance

    G4bool  startAtSurface_And_MoveEpsilon ;
    startAtSurface_And_MoveEpsilon =
             (stepData.GetPreStepPoint()->GetSafety() == 0.0) && 
             (stepData.GetStepLength() < kCarTolerance ) ;

    if( startAtSurface_And_MoveEpsilon) 
    {

       fLinearNavigator->
       LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(),
                                                  track.GetMomentumDirection(),
                                                  fCurrentTouchableHandle,
                                                  true                     );
       if( fCurrentTouchableHandle->GetVolume() != track.GetVolume() )
       {
         G4cerr << " ERROR: A relocation within safety has caused a volume change! " << G4endl  ; 
         G4cerr << "   The old volume is called " 
	         << track.GetVolume()->GetName() << G4endl ; 
         G4cerr << "   The new volume is called " ;

         if ( fCurrentTouchableHandle->GetVolume() != 0 )
	       {
	           G4cerr << fCurrentTouchableHandle->GetVolume()->GetName() << G4endl ; 
	       }
               else
	       {
	           G4cerr << "Out of World" << G4endl ; 
	       }
         G4cerr.precision(7) ;
         G4cerr << "   The position is " << track.GetPosition() <<  G4endl ;

          // Let us relocate again, for debuging

         fLinearNavigator->
         LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(),
                                                    track.GetMomentumDirection(),
                                                    fCurrentTouchableHandle,
                                                    true                     ) ;
         G4cerr << "   The newer volume is called "  ;

         if ( fCurrentTouchableHandle->GetVolume() != 0 )
	       {
	           G4cerr << fCurrentTouchableHandle->GetVolume()->GetName() << G4endl ;
	       } 
         else
	       {
	           G4cerr << "Out of World" << G4endl ; 
	       }
       }

       assert( fCurrentTouchableHandle->GetVolume()->GetName() ==
               track.GetVolume()->GetName() ) ;

       retCurrentTouchable = fCurrentTouchableHandle ; 
       fParticleChange.SetTouchableHandle( fCurrentTouchableHandle ) ;
       
    }
    else
    {
       retCurrentTouchable = track.GetTouchableHandle() ;
       fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
    }
    //  This must be done in the above if ( AtSur ) fails
    //  We also do it for if (true) in order to get debug/opt to  
    //  behave as exactly the same way as possible.

    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition()) ;
#else
    // ie #ifndef G4DEBUG does a quick relocation
    
    // The serves only to move the Navigator's location

    fLinearNavigator->LocateGlobalPointWithinVolume( track.GetPosition()) ;

    // The value of the track's current Touchable is retained. 
    //    (and it must be correct because we must use it below to
    //      overwrite the (unset) one in particle change)
    //  Although in general this is fCurrentTouchable, at the start of
    //   a step it could be different ... ??

    fParticleChange.SetTouchableHandle( track.GetTouchableHandle() ) ;
    retCurrentTouchable = track.GetTouchableHandle() ;
#endif

  }                   // endif ( fGeometryLimitedStep ) 

  const G4VPhysicalVolume* pNewVol = retCurrentTouchable->GetVolume() ;
  //G4cout<<pNewVol->GetName()<<std::endl;
  const G4Material* pNewMaterial   = 0 ;
 
  if( pNewVol != 0 ) pNewMaterial= pNewVol->GetLogicalVolume()->GetMaterial() ; 

  // ( <const_cast> pNewMaterial ) ;
#ifdef BEFORE_V7
  fParticleChange.SetMaterialChange( (G4Material *) pNewMaterial ) ;
#else
  fParticleChange.SetMaterialInTouchable( (G4Material *) pNewMaterial ) ;
#endif
  //    temporarily until Get/Set Material of ParticleChange, 
  //    and StepPoint can be made const. 
  // Set the touchable in ParticleChange
  //   this must always be done because the particle change always
  //   uses this value to overwrite the current touchable pointer.
  
  fParticleChange.SetTouchableHandle(retCurrentTouchable) ;

  return &fParticleChange ;
}

