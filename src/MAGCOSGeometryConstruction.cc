// 8/19/2014: Hexc, Olesya.  Cleaning up the code and understanding the geometry
// 
#include "MAGCOSGeometryConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ParticleTable.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "G4UserLimits.hh"
#include "MAGCOSGeometryMessenger.hh"
#include "MAGCOSMagneticField.hh"
#include "MAGCOSUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "MAGCOSSteppingAction.hh"
MAGCOSGeometryConstruction::MAGCOSGeometryConstruction()
{ 
  detectorMessenger = new MAGCOSGeometryMessenger(this);
  G4cout << " Making magnetic field ... " << G4endl;
  theMagnetosphereMagneticField = new MAGCOSMagneticField();
  G4cout << " Making magnetic field is done! " << G4endl;
  TopOfAtmosphere =19.99*km;
  RemoveEarth =false;
}
////////////////////////////////////////////////////////////////////////////////
//
MAGCOSGeometryConstruction::~MAGCOSGeometryConstruction()
{ delete detectorMessenger;
  delete theMagnetosphereMagneticField;
}
////////////////////////////////////////////////////////////////////////////////
//
G4VPhysicalVolume* MAGCOSGeometryConstruction::Construct()
{ 
  // Clean old geometry, if any
  //
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  G4double a,density;
  G4String name, symbol;
  G4int nel;
  
  //Vacuum
  
  a = 1.*g/mole; 
  density = 1.e-10*g/cm3;
  G4Element* elN = new G4Element( "Nitrogen", "N", 7. , 14.00674*g/mole );
  G4Element* elO = new G4Element( "Oxygen", "O", 8. , 15.9994*g/mole );
  
  G4Material* Vacuum =new G4Material(name="Vacuum",density,
				     nel=2,kStateSolid, 
				     293.0*kelvin, 
				     1.e-6*atmosphere );
  Vacuum->AddElement(elN, .7);
  Vacuum->AddElement(elO, .3);
  
  //Visualisation attributes
  G4VisAttributes * VisAttEarth = new G4VisAttributes(G4Colour(0.,0.,1.));
  VisAttEarth->SetVisibility(true);		  
  G4VisAttributes * VisAttWorld = new G4VisAttributes();
  VisAttWorld->SetVisibility(false);
  
  //------------------------------ 
  // World
  //------------------------------ 
  G4double MagSizeX=200*re;
  G4double MagSizeY=200*re;
  G4double MagSizeZ=200*re;
  solidWorld= new G4Box("World",MagSizeX/2.,MagSizeY/2.,MagSizeZ/2.);
  logicWorld= new G4LogicalVolume( solidWorld, Vacuum,
				   "World", 0, 0, 0);
  maxStepSize=5.*re;
  theWorldUserLimits=new G4UserLimits(maxStepSize);
  theWorldUserLimits->SetUserMaxTrackLength(120.*re);
  logicWorld->SetUserLimits(theWorldUserLimits);
  logicWorld->SetVisAttributes(VisAttWorld);
  physiWorld= new G4PVPlacement(0,               // no rotation
				G4ThreeVector(), // at (0,0,0)
				"WorldPV",       // its name
				logicWorld,      // its logical volume
				0,               // its mother  volume
				false,           // no boolean operations
				0);              // no field specific to volume
  
  //----------------------
  //Earth
  //-----------------------
  
  /* solidEarth=new G4Sphere("Earth",0. ,6378.16*km+top_altitude,0.,360.*deg,0.,180.*deg);*/
  if (!RemoveEarth){		      
    solidEarth=new G4Orb("Earth",6378.137*km+TopOfAtmosphere);
    logicEarth=new G4LogicalVolume(solidEarth,Vacuum,"Earth",0,0,0);
    theEarthUserLimits=new G4UserLimits(1.*km);
    logicEarth->SetUserLimits(theEarthUserLimits);
    theEarthUserLimits->SetUserMaxTrackLength(120.*re);
    logicEarth->SetVisAttributes(VisAttEarth);
    physiEarth=new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(0.,0.,0.), // at (0,0,0)
				 "Earth",       // its name
				 logicEarth,      // its logical volume
				 physiWorld,               // its mother  volume
				 false,           // no boolean operations
				 0);  
    
  }	
  return physiWorld;
}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetStopAltitude(G4double stop_altitude)
{
  if (stop_altitude <0.) {
    G4cout<<"The top altitude below which particle are stopped ";
    G4cout<<"should be >0."<<std::endl; 
    return;
  }
  else {
    TopOfAtmosphere = stop_altitude;
    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
    theSteppingAction->SetStopAltitude(stop_altitude);
    
  }	 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetRemoveEarth(G4bool aBool)
{RemoveEarth = aBool;
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetMaxStepLength(G4double maxStep)
{ 
  theWorldUserLimits->SetMaxAllowedStep(maxStep);
  maxStepSize=maxStep;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetMaxTrackLength
(G4double maxTrackLength)
{ 
  theWorldUserLimits->SetUserMaxTrackLength(maxTrackLength);
  theEarthUserLimits->SetUserMaxTrackLength(maxTrackLength);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSGeometryConstruction::SetMaxTimeOfFlight(G4double maxTime)
{
  theWorldUserLimits->SetUserMaxTime(maxTime);
  theEarthUserLimits->SetUserMaxTime(maxTime);
}

