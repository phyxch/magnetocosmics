
#ifndef MAGCOSPrimaryGeneratorAction_h
#define MAGCOSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "vector"



class G4ParticleGun;
class G4GeneralParticleSource;
class MAGCOSPrimaryMessenger;
class G4Event;
class G4PrimaryVertex;



class MAGCOSPrimaryGeneratorAction : 
          public G4VUserPrimaryGeneratorAction 
{
  public:
    MAGCOSPrimaryGeneratorAction();    
    ~MAGCOSPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void (MAGCOSPrimaryGeneratorAction::*pGeneratePrimaries)(G4Event* anEvent);
    void GenerateStandardPrimaries(G4Event* anEvent);
    void GeneratePrimariesForComputingCutOffRigidity(G4Event* anEvent);
   
    void SetRigidity(G4double aRigidity);
    
   
    G4bool SetPositionAndDirection(const G4String CoordSys, 
                                 const G4double anAlt,const G4double aLong,
				 const G4double aLat,const G4double aZenith,
				 const G4double  anAzimuth); 
    void SetPositionAndDirection(const G4String CoordSys, 
                                 const G4ThreeVector position,
				 const G4ThreeVector direction); 				  
    G4bool SetPosition(const G4String CoordSys, 
                     const G4double anAlt,const G4double aLong,
                     const G4double aLat);
    G4bool SetPosition(const G4String CoordSys, 
                     const G4ThreeVector aPosition);
    G4bool SetPositionOnDipoleMagneticShell
                     (const G4String Reference, const G4double L, 
		      const G4double latitude, const G4double longitude);  		      		     				 
    G4bool SetDirection(const G4String CoordSys, 
                      const G4double aZenith,const G4double anAzimuth);
    G4bool SetDirection(const G4String CoordSys, 
                      const G4ThreeVector aDirection);
    void DefineDirectionByPitchAngle(G4double pitch_angle, G4double phi); 	
				 
     
    
    void AddValuesToRigidityVector(G4int nvalues,G4double val1, G4double step);
    void SetDefaultRigidityVector();
    
    void SelectTypeOfPrimaries(G4String aString);
    void PrintPrimaryVertexInformation(const G4PrimaryVertex* aPrimary);
    void PrintBfieldAtPrimary();
    
    
    inline void ResetRigidityIndex() {rigidity_index=0;}
    inline G4int GetNumberOfRigidity() const {return rigidity_values.size();}
    inline void ResetRigidityVector(){rigidity_values.clear();}
    inline G4ThreeVector GetGEOPosition(){return GEOPosition;}
    inline G4ThreeVector GetGEODirection(){return GEODirection;}
    inline void GetGEOIDPosition(G4double& GEOIDalt,G4double& GEOIDlat,G4double& GEOIDlong)
                       {GEOIDalt=GEOIDaltitude;
		        GEOIDlat=GEOIDlatitude;
			GEOIDlong=GEOIDlongitude;}
    inline void GetGEOIDDirection(G4double& GEOIDzen,G4double& GEOIDaz)
                       {GEOIDzen=GEOIDzenith;
		        GEOIDaz=GEOIDazimuth;}			
    inline G4GeneralParticleSource* GetParticleSource()
                                  {return myParticleSource;}
    inline void SetVerbosity(G4int n){verbosity=n;}
    inline void SetInitialiseTrajectoryCSVBlockForForwardCase(bool aBool)
    				{ InitialiseTrajectoryCSVBlockForForwardCase=aBool;}
    inline void SetInitialiseTrajectoryCSVBlockForBackwardCase(bool aBool)
    				{ InitialiseTrajectoryCSVBlockForBackwardCase =aBool;}

  private:
    G4GeneralParticleSource* myParticleSource;
    MAGCOSPrimaryMessenger* myMessenger;
    
    
    G4ThreeVector GEOPosition;
    G4ThreeVector GEODirection;
    G4double GEOIDaltitude, GEOIDlongitude, GEOIDlatitude;
    G4double GEOIDzenith, GEOIDazimuth;
    G4double Rigidity;
    
    //member data for defining the binnning shema in rigidity for computing
    // rigidity filter
    //---------------------------------------------------------------
    std::vector<G4double> rigidity_values;
    G4int rigidity_index;
    
    G4int verbosity;
    G4bool InitialiseTrajectoryCSVBlockForForwardCase;
    G4bool InitialiseTrajectoryCSVBlockForBackwardCase;
  
    
   		  
   			  
     
};

#endif


