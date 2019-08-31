#ifndef MAGCOSApplicationScenario_h
#define MAGCOSApplicationScenario_h 1

#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4UserRunAction.hh"
#include"vector"

class MAGCOSScenarioMessenger;
class MAGCOSPrimaryGeneratorAction;
class MAGCOSEquationOfMotion;
class MAGCOSSpenvisManager;

class MAGCOSApplicationScenario : public G4UserRunAction 
{public:
  
  MAGCOSApplicationScenario();
  virtual ~MAGCOSApplicationScenario();
  
 public:
  virtual void BeginOfRunAction(const G4Run* aRun);
  virtual void EndOfRunAction(const G4Run* aRun);
  
  //apllication scenario 
  void Bline();
  void ParticleTrajectory();
  void ReverseParticleTrajectory();
  void ComputeRigidityFilter(G4String outputfile_name);
  void ComputeDirectionFilter(G4String CoordSys, G4double rigidity,
                               G4double cos_zen0, G4double delta_cos_zen, G4int nzen,
                               G4double azim0, G4double delta_azim, G4int nazim,
			       G4String outputfile_name);
  void RCutoffVsPosition(G4String CoordSys, 
                          G4double Altitude,
                          G4double lat0,G4double delta_lat, G4int nlat,
                          G4double long0, G4double delta_long, G4int nlong,
			  G4double zenith, G4double azimuth, 
			  G4String outputfile_name);
  void RCutoffVsSpenvisPositionGrid(G4String SpenvisCSVFileName,G4String outputfile_name,
  				    G4double zenith =0., 
				    G4double azimuth = 0.);			  
  void RCutoffVsSpenvisTrajectory(G4String SpenvisCSVFileName,G4String outputfile_name,
  				  G4double zenith =0., 
				  G4double azimuth = 0.);
  void RCutoffVsPositionOnDipoleMagneticShell(G4String CoordSys, 
                          G4double L,
                          G4double lat0,G4double delta_lat, G4int nlat,
                          G4double long0, G4double delta_long, G4int nlong,
			  G4double zenith, G4double azimuth, 
			  G4String outputfile_name);
  void RCutoffVsDirection(G4String CoordSys, G4double zen0, G4double delta_zen, G4int nzen,
                           G4double azim0, G4double delta_azim, G4int nazim, 
			   G4String outputfile_name);
   
  void RCutoffVsTime(G4double time0, G4double delta_t, G4int ntime, 
		     G4String outputfile_name);	
   
   	     		  
 
    
    
    
  void RegisterTrackLastPoint(G4ThreeVector position, 
                             G4ThreeVector momentum,G4double filter_value);
			     
  
  inline void SetAutomaticDetectionPenumbra(G4bool aBool)
   		 {AutomaticDetectionOfThePenumbra = aBool;}
  
  inline void SetRegisterResultsInSpenvisFile(bool aBool) 
  		{RegisterResultsInSpenvisFile = aBool;}		 
    
  inline void SetSpenvisFileName(G4String aName){SpenvisFileName = aName;}
   
   		       
  
 private:  
    
  void ComputeRigidityCutoff();
  void ComputeCutoff();
  
 private:
 
  MAGCOSScenarioMessenger* theMessenger;
   
  //Last position, momentum and filter value for a particle trajectory 
  std::vector<G4ThreeVector>  LastPositions;
  std::vector<G4ThreeVector>  LastMomentaOnCharge;
  std::vector<G4int>  FilterValues;
   
   
  //Rigidity cutoff
  G4double  Rs;
  G4double  Rc;
  G4double  Rm; 
   
  bool RegisterAsymptoticDirection;
    
 // automatic detection of the penumba
 // if true ComputeCutoff() is used
 // is false ComputeRigidityCutoff() is used
  bool AutomaticDetectionOfThePenumbra;  
  
  MAGCOSSpenvisManager* theSpenvisManager;   
  
  bool RegisterResultsInSpenvisFile;
  G4String SpenvisFileName;
  
  bool TotalDurationReached;
  
  
};












#endif
