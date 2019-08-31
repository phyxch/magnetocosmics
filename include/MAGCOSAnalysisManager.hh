
#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include"G4ios.hh"
#include"G4strstreambuf.hh"
#include"vector"
#include"globals.hh"
#include"fstream"
#include"G4ThreeVector.hh"


/*namespace AIDA
{class ITree;
 class IHistogramFactory;
 class IHistogram1D;
 class IHistogram2D;
 class IAnalysisFactory;
 class ITupleFactory;
 class ITuple;}

using namespace AIDA;*/ 
 
class G4Step;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MAGCOSAnalysisManager
{
public:
  
  ~MAGCOSAnalysisManager();
   static MAGCOSAnalysisManager* getInstance();
   
   //Asymptotic direction vs rigidity
   void OpenAsymptoticDirectionFile(G4String fileName);
   void RegisterAsymptoticDirection(G4ThreeVector LastPosition,
                                G4ThreeVector LastMomentaOnCharge,
			        G4int  FilterValue);
   void CloseAsymptoticDirectionFile(G4double Rc,G4double Rm,G4double Rs);
  
  //Asymptotic direction vs direction 
   void OpenAsymptoticDirVsDirFile(G4String coord_sys, G4double rigidity, G4String fileName);
   void RegisterAsymptoticDirVsDir(G4double zenith, G4double azimuth,G4ThreeVector LastPosition,
                                G4ThreeVector LastMomentaOnCharge,
			        G4int  FilterValue);
   
  
  //Rigidity ForDiffrent position
   void OpenCutoffVsPositionFile(G4String fileName,G4String CoordSys,
	                             G4double Altitude, G4double zenith,
				     G4double azimuth);
   void OpenCutoffVsSpenvisPositionGridFile(G4String fileName,G4String CoordSys,
	                             G4double zenith,
				     G4double azimuth);
   void OpenCutoffVsSpenvisTrajectoryFile(G4String fileName,G4String CoordSys,
	                             G4double zenith,
				     G4double azimuth);				     
   void RegisterCutoffVsPosition(G4double Rc,G4double Rm,G4double Rs,
                                     G4double latitude, 
				     G4double longitude, std::vector<G4double> Ilat);
   
   //Rigidity For Different position on a magnetic dipole shell
   void OpenCutoffVsPositionOnLShellFile(G4String fileName,G4String CoordSys,
	                             G4double L, G4double zenith,
				     G4double azimuth);
  
    //Rigidity For Different directions
   void OpenCutoffVsDirectionFile(G4String fileName,G4String CoordSys);
   void RegisterCutoffVsDirection(G4double Rc,G4double Rm,G4double Rs,
                                     G4double zenith, G4double azimuth);
   
   
   
    //Rigidity For Different times
   void OpenCutoffVsTimeFile(G4String fileName);
   void RegisterCutoffVsTime(G4double Rc,G4double Rm,G4double Rs,
                                     G4double time);
   
   
  
   void CloseAsciiFile();
private:
  static MAGCOSAnalysisManager* instance;
  std::ofstream theAsciiFile;

private:
  MAGCOSAnalysisManager();
   
 

private:

  
  

 /* IAnalysisFactory  *aFact;
  ITree             *theTree;
  IHistogramFactory *histFact;*/
 
  
  
 
};

#endif




