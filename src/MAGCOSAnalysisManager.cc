#include "MAGCOSAnalysisManager.hh"
#include "MAGCOSPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "SpaceCoordinateConvertor.hh"
#include "G4Step.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include "MAGCOSUnits.hh"


MAGCOSAnalysisManager* MAGCOSAnalysisManager::instance = 0;
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSAnalysisManager::MAGCOSAnalysisManager()  
{
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSAnalysisManager::~MAGCOSAnalysisManager() 
{
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSAnalysisManager* MAGCOSAnalysisManager::getInstance()
{
  if (instance == 0) instance = new MAGCOSAnalysisManager;
  return instance;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::OpenAsymptoticDirectionFile(G4String fileName)
{ // Asymptotic direction file
  theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"Rigidity"<<'\t'<<"Filter"<<'\t'<<"Asympt. Lat."<<'\t'
                    <<"Asympt. Long."<<'\t'
                    <<"Position Xgeo"<<'\t'<<"Ygeo"<<'\t'<<"Zgeo"<<G4endl;
		    
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::RegisterAsymptoticDirection(G4ThreeVector LastPosition,
                                G4ThreeVector LastMomentaOnCharge,
			        G4int  FilterValue)
{   //Rigidity 
    
    G4double rigidity =LastMomentaOnCharge.mag()/1000.;
    
    //Asymptotic direction
    
    G4ThreeVector AsympDir= -1.*LastMomentaOnCharge;
    G4double Aslat=90.-(AsympDir.theta()/degree);
    G4double Aslong=AsympDir.phi()/degree;
    
    
    //LastPosition
    
    G4ThreeVector LastGEOPosition;
    LastGEOPosition = LastPosition/re;
    
    //write on asymptotic direction file 
    
    theAsciiFile.precision(8);
    theAsciiFile.setf(std::ios::fixed);					   
    theAsciiFile<<rigidity<<'\t'<<'\t';
    theAsciiFile<<FilterValue<<'\t';
    theAsciiFile<<Aslat<<'\t'<<'\t'<<Aslong<<'\t'<<'\t'<<'\t';
    theAsciiFile<<LastGEOPosition.x()<<'\t'
                   <<LastGEOPosition.y()<<'\t'
		   <<LastGEOPosition.z()<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::OpenAsymptoticDirVsDirFile(G4String coord_sys,G4double rigidity,G4String fileName)
{ // Asymptotic direction file
  theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
  theAsciiFile<<"The coordinate system for defining arrival position is "<< coord_sys<<G4endl;
  theAsciiFile<<"rigidity: "<<rigidity /GV <<" GV " <<G4endl;
  theAsciiFile<<"zenith"<<'\t'<<"azimuth"<<'\t'<<"filter"<<'\t'<<"Asympt. Lat. GEO"<<'\t'
             <<"Asympt. Long. GEO"<<'\t'
             <<"Position Xgeo"<<'\t'<<"Ygeo"<<'\t'<<"Zgeo"<<G4endl;
		    
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::RegisterAsymptoticDirVsDir(G4double zenith, G4double azimuth,
                                                             G4ThreeVector LastPosition,
                                                             G4ThreeVector LastMomentaOnCharge,
			                                     G4int  FilterValue)
{  
    
    //Asymptotic direction
    
    G4ThreeVector AsympDir= -1.*LastMomentaOnCharge;
    G4double Aslat=90.-(AsympDir.theta()/degree);
    G4double Aslong=AsympDir.phi()/degree;
    
    
    //LastPosition
    
    G4ThreeVector LastGEOPosition;
    LastGEOPosition = LastPosition/re;
    
    //write on asymptotic direction file 
    
    theAsciiFile.precision(2);
    theAsciiFile.setf(std::ios::fixed);					   
    theAsciiFile<<zenith/degree<<'\t'<<azimuth/degree<<'\t'<<'\t';
    theAsciiFile<<FilterValue<<'\t';
    theAsciiFile<<Aslat<<'\t'<<'\t'<<Aslong<<'\t'<<'\t'<<'\t';
    theAsciiFile<<LastGEOPosition.x()<<'\t'
                   <<LastGEOPosition.y()<<'\t'
		   <<LastGEOPosition.z()<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::CloseAsymptoticDirectionFile
                                    (G4double Rc,G4double Rm,G4double Rs)				
{//write rigidity cutoff
 
  theAsciiFile<<"Rl "<<Rs<<'\t';
  theAsciiFile<<"Rc "<<Rc<<'\t';
  theAsciiFile<<"Ru "<<Rm<<G4endl;
  
 //close file
  theAsciiFile.close();
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::
            OpenCutoffVsPositionFile(G4String fileName,G4String CoordSys,
	                             G4double Altitude, G4double zenith,
				     G4double azimuth)
{ theAsciiFile.open(fileName, std::ios::out);
  theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"Altitude"<<'\t'<<Altitude/km<<G4endl;
  theAsciiFile<<"Zenith"<<'\t'<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth"<<'\t'<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Latitude"<<'\t'<<"Longitude"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"Ilat1"<<'\t'<<'\t'<<"Ilat2"<<'\t'<<'\t'<<"Ilat3"<<'\t'<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::OpenCutoffVsSpenvisPositionGridFile(G4String fileName,
							    G4String CoordSys,
	                             			    G4double zenith,
				     			    G4double azimuth)
{ theAsciiFile.open(fileName, std::ios::out);
  theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"Zenith [deg]:"<<'\t'<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth [deg]"<<'\t'<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Altitude [km]"<<'\t'<<"Latitude [deg]"<<'\t'<<"Longitude[deg]"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"Ilat1"<<'\t'<<"Ilat2"<<'\t'<<"Ilat3"<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::OpenCutoffVsSpenvisTrajectoryFile(G4String fileName,G4String CoordSys,
	                             G4double zenith,
				     G4double azimuth)
{ theAsciiFile.open(fileName, std::ios::out);
  theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"Zenith [deg]:"<<'\t'<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth [deg]:"<<'\t'<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Mod Julian Day"<<'\t'<<"Altitude [km]"<<'\t'<<"Latitude [deg]"<<'\t'<<"Longitude[deg]"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"ILat1"<<'\t'<<"ILat2"<<'\t'<<"ILat3"<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::RegisterCutoffVsPosition(G4double Rc,G4double Rm,G4double Rs,
                                  G4double latitude, G4double longitude, std::vector<G4double> ILat)
{ 
 
 theAsciiFile<<latitude/degree<<'\t'<<longitude/degree<<'\t';
 if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<ILat[0]/degree<<'\t'<<ILat[1]/degree<<'\t'<<ILat[2]/degree<<'\t'<<'\t';
 theAsciiFile << Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
				       
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::
            OpenCutoffVsPositionOnLShellFile(G4String fileName,G4String CoordSys,
	                             G4double L, G4double zenith,
				     G4double azimuth)
{ 
  theAsciiFile.open(fileName, std::ios::out);
  theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4); 
  
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
  theAsciiFile<<"L:"<<'\t'<<L<<G4endl;
  theAsciiFile<<"Zenith:"<<'\t'<<zenith/degree<<G4endl;
  theAsciiFile<<"Azimuth:"<<'\t'<<azimuth/degree<<G4endl; 
  theAsciiFile<<"Latitude"<<'\t'<<"Longitude"<<'\t'; 
  if (getenv("INVARIANT_LATITUDE")) theAsciiFile<<"Ilat1"<<'\t'<<"Ilat2"<<'\t'<<"Ilat3"<<'\t';
  theAsciiFile<<"Ru"<<'\t'<<'\t'<<"Rc"<<'\t'<<'\t'<<"Rl"<<G4endl; 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::
            OpenCutoffVsDirectionFile(G4String fileName,G4String CoordSys)
{ theAsciiFile.open(fileName, std::ios::out);
  theAsciiFile<<std::setiosflags(std::ios::scientific);;
  theAsciiFile<<std::setprecision(4);
  theAsciiFile<<"System of coordinate: "<<'\t'<<CoordSys<<G4endl;
 
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
      (MAGCOSPrimaryGeneratorAction*)
          G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
 
  if (CoordSys=="GEOID"){
  	G4double GEOIDalt;
   	G4double GEOIDlat;
   	G4double GEOIDlong;
   	thePrimaryGenerator->GetGEOIDPosition(GEOIDalt,GEOIDlat,GEOIDlong);
   	theAsciiFile<<"Altitude:"<<'\t'<<GEOIDalt/km<<" km"<<G4endl;
   	theAsciiFile<<"Latitude:"<<'\t'<<GEOIDlat/degree<<" deg"<<G4endl;
   	theAsciiFile<<"Longitude:"<<'\t'<<GEOIDlong/degree<<" deg"<<G4endl;  
  }
  else{ 
   	G4ThreeVector GEOposition = thePrimaryGenerator->GetGEOPosition();
    	G4ThreeVector position=SpaceCoordinateConvertor::getInstance()
                                  ->Transform(GEOposition,"GEO",CoordSys)/re;
                   
    	theAsciiFile<<"Position:"<<'\t'<<position.x()<<'\t'
                                   <<position.y()<<'\t'
				   <<position.z()<<G4endl;
    
  }
  theAsciiFile<<"Zenith"<<'\t'<<"Azimuth"<<'\t'<<"Ru"<<'\t'
                    <<"Rc"<<'\t'<<"Rl"<<G4endl; 
 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::RegisterCutoffVsDirection
                                    (G4double Rc,G4double Rm,G4double Rs,
                                     G4double zenith, G4double azimuth)
{ theAsciiFile<<zenith/degree<<'\t'<<'\t'<<azimuth/degree<<'\t'<<
                                            Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
 					       
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::
            OpenCutoffVsTimeFile(G4String fileName)
{ theAsciiFile.open(fileName, (std::ios::binary | std::ios::out));
 
  MAGCOSPrimaryGeneratorAction* thePrimaryGenerator =
 	(MAGCOSPrimaryGeneratorAction*)
          	G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
 
  G4ThreeVector GEOposition = thePrimaryGenerator->GetGEOPosition()/re;
  G4ThreeVector GEOdirection = thePrimaryGenerator->GetGEODirection();
  
  
  // G4ThreeVector position=SpaceCoordinateConvertor::getInstance()
  //                                 ->Transform(GEOposition,"GEO","CoordSys")/re;
                   
  theAsciiFile<<"GEO Position:"<<'\t'<<GEOposition.x()<<'\t'
                                     <<GEOposition.y()<<'\t'
				     <<GEOposition.z()<<G4endl;
  theAsciiFile<<"GEO Direction:"<<'\t'<<GEOdirection.x()<<'\t'
                                     <<GEOdirection.y()<<'\t'
				   <<GEOdirection.z()<<G4endl;
   
  theAsciiFile<<"Time"<<'\t'<<"Ru"<<'\t'
                    <<"Rc"<<'\t'<<"Rl"<<G4endl; 
  theAsciiFile.precision(2);
  theAsciiFile.setf(std::ios::fixed);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::RegisterCutoffVsTime
                                    (G4double Rc,G4double Rm,G4double Rs,
                                     G4double time)
{ theAsciiFile<<time/s<<'\t'<<Rm<<'\t'<<Rc<<'\t'<<Rs<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSAnalysisManager::CloseAsciiFile()
{ theAsciiFile.close();
}
