#ifndef MAGCOSSpenvisManager_HH
#define MAGCOSSpenvisManage_HH 1
#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include"vector"
#include"globals.hh"
#include"fstream"
#include"G4ThreeVector.hh"
#include "SpenvisCSV.hh"
#include "SpenvisCSVCollection.hh"

 
class G4Step;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MAGCOSSpenvisManager
{
public:

   ~MAGCOSSpenvisManager(); 
   static MAGCOSSpenvisManager* getInstance();
 
   
   //General
   //--------
   void WriteData(std::vector< double>,SpenvisCSV* aBlock = 0 ); 
   void SaveCSVFile(G4String name);
   void ClearCSVFile();
   
   void WriteBfieldInformation(SpenvisCSV* aBlock = 0 );
   
   
   //Asymptotic direction vs rigidity
   //--------------------------------
   void InitialiseAsymptoticDirectionCSVFile();
   void RegisterAsymptoticDirection(G4ThreeVector LastPosition,
                                    G4ThreeVector LastMomentaOnCharge,
			            G4int  FilterValue);
   
  
  
  //Rigidity ForDifferent position
  //-----------------------------
   void InitialiseCutoffVsPosCSVFile(G4String coor_sys, 
   				     G4double zenith,
				     G4double azimuth);
   void InitialiseCutoffVsDirCSVFile(G4String coor_sys, 
   				     G4double altitude,
				     G4double latitude,
				     G4double longitude);
   
   void InitialiseCutoffVsTimeCSVFile(G4String coor_sys, 
   				     G4double altitude,
				     G4double latitude,
				     G4double longitude,
				     G4double zenith,
				     G4double azimuth);				     				     
   
   void AddMetaVariableToCSVFile(G4String variable_name,
    				 G4String variable_unit,
				 G4double value,
				 G4String format,
				 SpenvisCSV* aBlock = 0 );
   void ClearCSVFiles();
   
   
   
   //Reading spenvis input file
   
   bool ReadSpenvisPositionGridCSVFile(G4String input_file,
   				   std::vector<double>& altitude,
   				   std::vector<double>&  latitude,
				   std::vector<double>&  longitude,
				   G4int nblock =0);
   bool ReadSpenvisTrajectoryCSVFile(G4String input_file,
   				   std::vector<double>&  altitude,
   				   std::vector<double>&  latitude,
				   std::vector<double>&  longitude,
				   std::vector<double>&  smjd,
				   G4int nblock =0);
   
// particle trajectory or Bfield line 
   void InitialiseTrajectoryCSVBlock(G4String Type);
   void RegisterTrajectoryPoint(G4ThreeVector GEOPosition);
   void SaveTrajectoryCSVBlocks(G4String name);
   
     
   inline void ResetTheParticleTrajectoryCSVBlocks()
   				{theParticleTrajectoryCSVBlocks.clear();}
   inline void SetRegisterParticleTrajectory(bool aBool){RegisterParticleTrajectory =aBool;};				   
   inline bool GetRegisterParticleTrajectory(){return RegisterParticleTrajectory;}
private:
   static MAGCOSSpenvisManager* instance;
   SpenvisCSVCollection* theSpenvisCSVCollection;
   std::vector<SpenvisCSV> theSpenvisCSVFiles;
   std::vector<G4String> block_names;
   
   std::vector<SpenvisCSV> theParticleTrajectoryCSVBlocks;
   bool RegisterParticleTrajectory;
   MAGCOSSpenvisManager();  
 
};

#endif




