// Updated on 9/16/2014, Hexc, Olesya: added CLHEP namespace for g4.10 version
#ifndef SPACECOORDINATECONVERTOR_HH
#define SPACECOORDINATECONVERTOR_HH 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              SpaceCoordinateConvertor.hh
//
// Version:		VERSION_NUMBER
// Date:		LAST_DATE
// Author:		L Desorgher
// Organisation:	ORGANISATION_NAME
// Project:		PROJECT_NAME
// Customer:		CUSTOMER_NAME
// Contract:		CONTRACT_NAME
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// This class allows to convert vector component  in one space 
//  coordinate system to another one, at user selected date. It computes also
// the goemagnetic dipole parameters (Bo, Dipole Axis, tilt angle and schift)
// from the igrf coefficient.
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// SpaceCoordinateConvertor()
//    Constructor.
//
// ~SpaceCooordinateConvertor()
//    Destructor.
//
// void SetReferenceDate(G4int year,G4int month ,G4int day,
//	                       G4int hour,G4int minute,G4int second)
//    Set the universal time of reference. The coordinate conversion matrices and
//    the geomagnetic dipole parameters (B0, dipole axis and tilt angle )
//    are recomputed for this date. 
// 
// void SetReferenceDate(struct tm  ref_date)
//    Same as the precedent member function exepect that the format of the 
//    date input is the structure tm. 
//  
//			           
// void SetSystemInAndOut(G4String sys_in, G4String sys_out)
//    Set system_in and system_out 
//
// struct tm GetReferenceDate() 
//    return ReferenceDate. The type of reference date is the structure tm.
// 
// G4double GetTiltAngle()
//    return the tilt angle of the geomagnetic dipole at the time defined by 
//    the attribute ReferenceDate
// 
// G4ThreeVector GetGeoDipoleAxisInGEO()
//    return the goemagnetic dipole axis components in GEO ccordinate at time
//    defined by the attribute ReferenceDate 
//
// G4ThreeVector GetGeoDipoleSchiftInGEO()
//    return the schift of the geomagnetic dipole from the Earth's center 
//    in GEO ccordinate at time defined by the attribute ReferenceDate 
//
// G4double GetGeoDipoleB0()
//    return the geomagnetic dipole momentum B0  
//    at the time defined by the attribute ReferenceDate 
//
// G4ThreeVector Transform(G4ThreeVector vec_in, G4String sys_in, G4String sys_out)
//    It converts the component of the vector vector_in from the coordinate 
//    system  sys_in into the coordinate system sys_out. 
//    It returns the component in the sys_out coordinate system.
//    Possible values for sys_in and sys_out are GEO, GSM and MAG.
//    The attribute system_in and system_out are set to sys_in and sys_out 
//    respectively.  
// 
// G4ThreeVector Transform(G4ThreeVector vec_in)
//    It converts the component of the vector vec_in from the coordinate 
//    system  system_in into the coordinate system system_out. System_in and 
//    system_out are class attributes that can be defined by  SetSystemInAndOut
//    It returns the component in the sys_out coordinate system.
//  
// G4ThreeVector TransformGEOinGSM(G4ThreeVector vec_in)
//    It converts the component of the vector vec_in from the GEO coordinate 
//    system   into the GSM coordinate system.
//    It returns the vector components in the GSM coordinate system.
//
// G4ThreeVector TransformGSMinGEO(G4ThreeVector vec_in)
//    It converts the component of the vector vec_in from the GSM coordinate 
//    system   into the GEO coordinate system.
//    It returns the vector components in the GEO coordinate system.
// 
// void ComputeGEODirectionAndPositionFromGEOID
//	            (G4double altitude, G4double latitude, G4double longitude,
//		     G4double zenith, G4double azimuth,
//		     G4ThreeVector& position, G4ThreeVector& direction)
//     Compute the GEO ccordiante system components of the direction and 
//     position of incidence a particle  given in GEIOD coordiante system by the altitude, latitude
//     lobgitude, zenith and azimuth value. 
//
// void ComputeGEOIDCoordinatesFromGEOPosition
//	            (const G4ThreeVector GEOPosition,
//                   G4double& altitude, G4double& longitude, G4double& latitude)
//     Compute the GEOID altitude, latitude and longitude corrsponding to a
//     geographic position given by GEOPosition. 
//       

  
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "vector"
#include "DateAndTime.hh"

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

// the Earth's ellipsoid is taken from IAU(1964)
/*static const G4double GeoidLargeSemiaxis2=6378.16*6378.16*km*km;
static const G4double GeoidSmallSemiaxis2=6356.775*6356.775*km*km;*/

//New ellipsoid WGS84 for 9th IGRF generation
static const G4double GeoidLargeSemiaxis2=6378.137*6378.137*km*km;
static const G4double GeoidSmallSemiaxis2=6356.752*6356.752*km*km;

class SpaceCoordinateConvertor 
{
public:
		       
        
        ~SpaceCoordinateConvertor() ;
	
	
	 static  SpaceCoordinateConvertor* getInstance(); 
	 
	 //Set methods
	 void SetReferenceDate(G4int year,G4int month ,G4int day,
	                       G4int hour,G4int minute,G4int second);
			       
	 void SetReferenceDate(DateAndTime ref_date);    
       	 void SetSystemInAndOut(G4String sys_in, G4String sys_out);
	 inline void Setp_hh_nm(std::vector<float> *p) {p_hh_nm= p;}
	 inline void Setp_gg_nm(std::vector<float> *p) {p_gg_nm= p;}
	 inline void SetVerbosity(G4int n) {verbose=n;}
	 
	 //Get Methods
	 
	 inline DateAndTime GetReferenceDate() {return ReferenceDate;} 
	 
	 inline G4double GetTiltAngle()
	                     {return tilt_angle;}
	 inline G4ThreeVector GetGeoDipoleAxisInGEO()
	                               {return geodipole_axis_in_GEO;}
	 inline G4ThreeVector GetGeoDipoleSchiftInGEO()
	                               {return geodipole_schift_in_GEO;}
	 inline G4double GetGeoDipoleB0()
	                               {return geodipole_B0;}
				       
	
	 G4ThreeVector Transform(G4ThreeVector vec_in,
	                         G4String sys_in,
			         G4String sys_out);
	 G4ThreeVector Transform(G4ThreeVector vec_in);				   
	 G4ThreeVector TransformGEOinGSM(G4ThreeVector vec_in);
	 G4ThreeVector TransformGSMinGEO(G4ThreeVector vec_in);
	 
	 
	 
	 void ComputeGEOPositionFromGEOID
	            (G4double altitude, G4double latitude, G4double longitude,
		     G4ThreeVector& position);
	 
	 void ComputeGEODirectionAndPositionFromGEOID
	            (G4double altitude, G4double latitude, G4double longitude,
		     G4double zenith, G4double azimuth,
		     G4ThreeVector& position, G4ThreeVector& direction);				    	 
         
	 void ComputeGEOIDCoordinatesFromGEOPosition
	             (const G4ThreeVector GEOposition, 
		      G4double& altitude, G4double& longitude, G4double& latitude);  	 	
	void ComputeGEOIDCoordinatesFromGEOPositionAndDirection
	             (const G4ThreeVector GEOposition,const G4ThreeVector GEOdirection, 
		      G4double& altitude, G4double& longitude, G4double& latitude,
		      G4double& zenith,G4double& azimuth);
	 double ComputeGEODipoleInvariantLatitude(G4String sys, G4ThreeVector position);
protected:


private: 
        static SpaceCoordinateConvertor* instance; 

private:
        SpaceCoordinateConvertor(); 	

private:
        DateAndTime ReferenceDate;
	
	// sun position in GEI
	//---------------------
	
	G4double gst; // greenwich sideral time
	G4double sun_rasn; // sun longitude along ecliptic 
	G4double sun_long; // sun right ascension 
	G4double sun_dec;  // sun declination 
	G4double obliq; //ZGSE in GEI is (0.,-sin(obliq),cos(obliq))
	
	//Geomagnetic dipole
	//------------------
	
	G4double tilt_angle;
	G4ThreeVector geodipole_axis_in_GEO; 
	G4ThreeVector geodipole_schift_in_GEO;
	G4double  geodipole_B0; 
	
	//Rotation matrix
	//----------
	
	G4RotationMatrix GEItoGEO,GEItoGSE, GEOtoMAG,GEOtoGSM,GEOtoSM,GEOtoGSE;
	G4RotationMatrix SelectedMatrix;
	
	//Coordinate system 
	//----------------
	 G4String system_in,system_out;
	 
	// Pointer on IGRF coefficient at reference data
	//------------------------------------
	

	std::vector<float>  *p_hh_nm;
	std::vector<float>  *p_gg_nm; 
	
	G4int verbose;
	 
	
	//compute the sun position in GEI from :
	//          Russell, C., T., Cosmic Electrodynamics, 1971, V.2, 184-196
	// valid from 1901-2099
	void ComputeSunPosition(); 
	
	//Compute the geomagnetic dipole parameters 
	// from   IGRF/DGRF parameters valid for the date of reference
	// valid only for the period 1965-2005
	
	void ComputeGeodipoleParameters(); 
	void ComputeGeodipoleParameters1();
	
	//Compute the transformation coordinate matrices 
        void ComputeTranformationMatrices();
	
	//Compute the selected matrix of transformation according to 
	// system_in and system_out 
        void ComputeSelectedMatrix();
	
		 
	//Recompute everything after changing the date
	void Initialise();

} ;

#endif
