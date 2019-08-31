// Updated on 8/6/2014 Hexc, Olesya
// Fixing a missing variable: kCarTolerence
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
// $Id: MYTransportation.hh,v 1.2 2014/08/06 15:59:44 hexc Exp $
// GEANT4 tag $Name:  $
//
// 
// ------------------------------------------------------------
//	GEANT 4  include file implementation
//
// ------------------------------------------------------------
//
//   This class' object is a process responsible for the transportation of 
// a particle, ie the geometrical propagation that encounters the 
// geometrical sub-volumes of the detectors.
//
//   It is also tasked with part of updating the "safety".
//
// =======================================================================
// Created:  19 March 1997, J. Apostolakis
// =======================================================================
#ifndef MYTransportation_hh
#define MYTransportation_hh 1

#include "G4VProcess.hh"
#include "G4FieldManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
//#include "G4PropagatorInField.hh"
#include "MYPropagatorInField.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleChangeForTransport.hh"

class MYTransportation : public G4VProcess 
{
  // Concrete class that does the geometrical transport 
  
public:
  MYTransportation();
  ~MYTransportation(); 
  
  //  G4double          GetContinuousStepLimit  (
  G4double      AlongStepGetPhysicalInteractionLength(
						      const G4Track& track,
						      G4double  previousStepSize,
						      G4double  currentMinimumStep, 
						      G4double& currentSafety,
						      G4GPILSelection* selection
						      );
  
  G4VParticleChange* AlongStepDoIt(
				   const G4Track& track,
				   const G4Step& stepData
				   );
  
  // This only does the relocation
  // 
  G4VParticleChange* PostStepDoIt(
				  const G4Track& track,
				  const G4Step&  stepData
				  );
  
  // This forces the PostStepDoIt action to be called, 
  //   but does not limit the step.
  //
  G4double PostStepGetPhysicalInteractionLength(
						const G4Track& ,
						G4double   previousStepSize,
						G4ForceCondition* pForceCond
						);
  
  //  Access/set the assistant class that Propagate in a Field
  MYPropagatorInField* GetPropagatorInField();
  void SetPropagatorInField( MYPropagatorInField* pFieldPropagator);
  
  //  no operation in  AtRestDoIt
  G4double AtRestGetPhysicalInteractionLength(
					      const G4Track& ,
					      G4ForceCondition* 
					      ) { return -1.0; };
  
  //  no operation in  AtRestDoIt
  G4VParticleChange* AtRestDoIt(
				const G4Track& ,
				const G4Step&
				) {return NULL;};
protected:
  //  Checks whether a field exists for the "global" field manager.
  G4bool               DoesGlobalFieldExist();
  
private:
  // The Propagators used to transport the particle
  G4Navigator*         fLinearNavigator;
  //G4PropagatorInField* fFieldPropagator;
  MYPropagatorInField* fFieldPropagator;
  // Field Manager for the whole Detector
  // G4FieldManager*      fGlobalFieldMgr;     // Used MagneticField CC
  
  // The particle's state after this Step, Store for DoIt
  G4ThreeVector        fTransportEndPosition;
  G4ThreeVector        fTransportEndMomentumDir;
  G4double             fTransportEndKineticEnergy;
  G4ThreeVector        fTransportEndSpin;
  G4bool               fMomentumChanged;
  G4bool               fEnergyChanged;
  G4bool               fEndGlobalTimeComputed; 
  G4double             fCandidateEndGlobalTime;
  
  G4bool               fParticleIsLooping;
  
  G4TouchableHandle    fCurrentTouchableHandle;
  
  // Whether a magnetic field exists ...
  // G4bool         fFieldExists;
  //   The above data member is problematic: it is useful only if
  // it is initialised.  However the transportation process(es) are not
  // capable of doing this initialisation itself (themselves) and there
  // seems no alternative agent capable of doing it right now.
  // Eg, at construction time it is likely that the field manager has
  // not yet been informed about the detector's field 
  // I cannot foresee how the transportation can be informed later. JA 
  //   The current answer is to ignore this data member and use 
  // the member function DoesGlobalFieldExist() in its place ...
  //    John Apostolakis, July 7, 1997
  
  // Flag to determine whether a boundary was reached.
  G4bool fGeometryLimitedStep;
  
  // Remember last safety origin & value.
  G4ThreeVector  fPreviousSftOrigin;
  G4double       fPreviousSafety; 
  
  // New ParticleChange
  G4ParticleChangeForTransport fParticleChange;
  
  G4double endpointDistance;
  G4double kCarTolerance;
};


#include "MYTransportation.icc"

#endif  

// End of MYTransportation.hh
