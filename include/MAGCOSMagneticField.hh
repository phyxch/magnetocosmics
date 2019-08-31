#ifndef MAGCOSMAGNETICFIELD_HH
#define MAGCOSMAGNETICFIELD_HH 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              MAGCOSMagneticField.hh
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
// This class defines the magnetic field used in the simulation. 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// MAGCOSMagnetiField()
//    Constructor.
//
// ~MAGCOSMagnetiField()
//    Destructor.
//
// void GetFieldValue(const G4double yTrack[] ,G4double B[]) const
//    Define the total magnetic field vector B at the position defined by yTrack
//
// bool OutsideMagnetosphere(G4ThreeVector pos)
//    return true if the position defined by pos is outside the magnetosphere
//    return false if the position definde by pos is inside the magnetosphere
//
// void SetTimeOfB(G4double val);
//    Define TimeOfB
//    TimeOfB  represents the number of second after the StartDate 
//    at which the magnetic field will be computed   	
// 
// void SetStartDate(G4int year, G4int month,G4int day,
//	                       G4int hour, G4int minute,G4int second );
//    Define StartDate 	
//
// void SetIopt(int iopt)
//    Set the iopt parameter used  in the tsyganenko89 model, 
//    defining the geomagnetic activity level
//    IOPT=  1       2        3        4        5        6      7
//                 corresponds to :
//    KP= 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  > =6-
//    Default value is 1
// 
// void SetInternalField(G4String aString)
//    Defines the geomagnetic field model
//    Possible values for aString are IGRF, and GEODIPOLE
//    
// void SetExternalField(G4String aString)
//    Defines the Earth's external magnetospheric magnetic field model
//    Possible values for aString are NOFIELD, TSY89, TSY96 and TSY2001
//
// void SetMagnetopauseModel(G4String aString)
//    Defines the magnetopause model 
//    Possible values are IGRF, TSY89, TSY96 and TSY2001  
//  
// void Setnm_igrf(G4int n)
//     Defines the maximum  order n of harmonics to be considered in IGRF
//     min value 1
//     max value 10
//     Default is 10
//
// void SetDipoleB0(G4double val)
//     Set the geomagnetic dipole momentum BO
//      
// void SetDipolePS(G4double aVal)
//     Set the geomagnetic dipole tilt angle, the dipole axis is computed
//     according to this value
// 
// void SetDipoleAxis(G4double Theta,G4double Phi)
//     Set the geomagnetic dipole axis, the tilt angle is computed according
//     to this value 
//
// void SetDipoleSchift(G4ThreeVector aVec)
//     The center of the geomagnetic dipole is schifted from the Earth's
//     centrum by aVec in GEO coordinate 
// 
// void SetTiltedDipoleParameterFromIGRF()
//     The geomagnetic dipole parameter: Bo, axis and tilt angle are deduced  
//     from the 1st order IGRF cooeficientd for the date define by StarDate and
//     timeOfB. The Geomagnetic dipole is not schifted from the Earth's center
//
// void SetEccentricDipoleParameterFromIGRF()
//     The geomagnetic dipole parameter: Bo, dipoleAxis and tilt_angle 
//     are dededuced from the 1st order IGRF cooeficients for the date define by StarDate and
//     timeOfB. The Geomagnetic dipole is translated from the Earth's center
//     by a vector eddeduced from the 1st and 2nd order IGRF cooeficients.
//
// void SetPdyn(G4double aVal)
//     Set the solar wind dynamic pressure Pdyn (nV^2) parameters used in the tsyganenko
//     96 and 2001 model.
//
// void SetDst(G4double aVal)
//     Set the disturbance storm index (Dst) prameter used in the Tsyganenko 96
//     and 2001 model
// 
// void SetImfy(G4double aVal)
//     Set the GSM y component of the interplanetary magnetic field. This
//     parameter is used in the Tsyganenko96 and 2001 model.
//
// void SetImfz(G4double aVal)
//     Set the GSM z component of the interplanetary magnetic field. This
//     parameter is used in the Tsyganenko96 and 2001 model.
// 
// void SetG1(G4double aVal)
//     Set the G1 parameter used in the Tsyganenko 2001 model
// 
// void SetG2(G4double aVal)
//     Set the G2 parameter used in the Tsyganenko 2001 model
//
// void SetEpsilon(G4double aVal)
//     Set the maximum relative error for integration.
//
// void SetDeltaIntersection(G4double aVal)
//     Set the precison for crossing boundary when integrating the equation
//     of motion.
//
// void SetDeltaChord(G4double aVal)
//     Set the maximum infinitesimal step allowed by the ChordFinder.
//
// void ResetIntegrationParameters()
//     Reset the integration parameters to their default values
//  
// void PrintStormParameter()
//     Print the value of the parameters iopt, Pdyn, Dst, Imfy, Imfz, G1,G2 
//     used in the Tsyganenko 96 and 2001 model 
//
// void ReadTSYParameter(G4string filename)
//     Read the  start date and series of the parameters iopt, Pdyn, Dst, Imfy, Imfz, 
//     G1 and G2 in function of time t. The time t is defined in second 
//     form the start date
//
// void ReadTSYParameter(G4string filename)
//     Read the  start date and series of the parameters iopt, Pdyn, Dst, Imfy, Imfz, 
//     G1 and G2 in function of time t. The time t is defined in second 
//     form the start date 
    

  
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
#include"G4ios.hh"

#include "G4MagneticField.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"
#include "vector"
#include"G4strstreambuf.hh"
#include "DateAndTime.hh"

#ifdef  USE_UNILIB
#include"unilib_c_pub.h"
#endif

class MAGCOSEquationOfMotion;
class MAGCOSFieldMessenger;
class G4ChordFinder;
class G4MagIntegratorStepper;
class SecondInvariant;

class MAGCOSMagneticField : public G4MagneticField
{
public:
	 //constructor destructor	       
         MAGCOSMagneticField() ;
	 
	 ~MAGCOSMagneticField() ;
	              
	 // Gives the magnetic field B at a given position defined by yTrack
	 void GetFieldValue(const G4double yTrack[],G4double B[]) const;
	 G4ThreeVector GetFieldValue(const G4ThreeVector GEOposition) const; 
	 // magnetopause
	 bool OutsideMagnetosphere(G4ThreeVector pos) const;
	 G4ThreeVector FindPositionOnMagnetopause(G4double theta,
	                                          G4double xgsm,
						  G4double precision) const;
	 G4ThreeVector FindStandOffPosition(G4double precision) const;
	 
	 //Set methods 
	 void SetTimeOfB(G4double val);	
	 void SetStartDate(G4int year, G4int month,G4int day,
	                   G4int hour, G4int minute,G4int second );
	 inline void SetStartDate(DateAndTime aDate) { StartDate = aDate;}
	 void SetIopt(G4int val); 
	 void SetInternalField(G4String aString);
	 inline G4String GetIntFieldModelName() const{return IntFieldModelName;}
	 void SetExternalField(G4String aString);
	 inline G4String GetExtFieldModelName() const{return ExtFieldModelName;}
	 void SetMagnetopauseModel(G4String aString);
	 inline void Setnm_igrf(G4int aVal) {nm_igrf=aVal;}
	 void SetDipoleB0(G4double aVal);
	 inline void SetDipolePS(G4double aVal){DipolePS=aVal;}
	 inline void SetDipoleAxis(G4double Theta,G4double Phi){DipoleTheta=Theta;
	                              				DipolePhi=Phi;}
	 inline void SetDipoleSchift(G4ThreeVector aVec){DipoleSchift=aVec;}
	 void SetStepper(G4String aString);
	 void ComputeBfieldAtDifferentPosition(G4String cosys_pos, 
				G4double alt0, G4double dAlt, G4int nAlt,
				G4double lat0, G4double dlat, G4int nlat,
				G4double long0, G4double dlong, G4int nlong,
				G4String file_output) const;
	 
	 // Define the geomagnetic dipole parameter according to IGRF
	 // coefficients
	 void SetTiltedDipoleParameterFromIGRF();
	 void SetEccentricDipoleParameterFromIGRF();		 
	 
	 //Method for geomagnetic and sw parameters used in the Tsyganenko
	 //models 
	 inline void SetPdyn(G4double aVal){Pdyn = aVal;}
	 inline void SetDst(G4double aVal){Dst = aVal;}
	 inline void SetImfx(G4double aVal){Imfy = aVal;}
	 inline void SetImfy(G4double aVal){Imfy = aVal;}
	 inline void SetImfz(G4double aVal){Imfz = aVal;}
	 inline void SetStdOff(G4double aVal){StdOff = aVal;}
	 inline void SetG1(G4double aVal){G1 = aVal;}
	 inline void SetG2(G4double aVal){G2 = aVal;}
	 inline void SetNsw(G4double aVal){Nsw = aVal;}
	 inline void SetVsw(G4double aVal){Vsw = aVal;}
#ifdef  USE_TSY04
	 inline void SetW1(G4double aVal){W1 = aVal;}
	 inline void SetW2(G4double aVal){W2 = aVal;}
	 inline void SetW3(G4double aVal){W3 = aVal;}
	 inline void SetW4(G4double aVal){W4 = aVal;}
	 inline void SetW5(G4double aVal){W5 = aVal;}
	 inline void SetW6(G4double aVal){W6 = aVal;}
#endif 	
	 void ResetIntegrationParameters();
	 void SelectBulirshStoerMethod();
	 void SelectG4IntegrationMethods();
	 void SetEpsilon(G4double aVal);
	 void SetDeltaChord(G4double aVal);
	 void SetBSMaxStep(G4double aVal);
	 void SetDeltaIntersection(G4double aVal); 
	 
	 void PrintStormParameter();
	 void ReadTSYParameter(G4String nameFile);
	 void ReadIgrfTable(G4String name_file);
	 void ComputeIgrfCoefAccordingToTime();
	 void PrintBfield(G4ThreeVector geo_pos) const; 
	 
	 //Get methods
	 
	 inline MAGCOSEquationOfMotion* GetEquationOfMotion()
	                                {return fEquationOfMotion;}
	 
	 inline G4ThreeVector GetDipoleSchift() const
	                                {return DipoleSchift;}
	 				
	 inline G4double GetDipoleB0() const
	                                {return DipoleB0;}
	 inline G4double GetDipolePhi() const
	                                {return DipolePhi;}
	 inline G4double GetDipoleTheta() const
	                                {return DipoleTheta;}
	 inline G4double GetDipolePS() const 
					{return DipolePS;}
	 
	 inline DateAndTime GetReferenceDate() const
	 				{return ReferenceDate;}
	 inline G4double GetPdyn() const{ return Pdyn ;}
	 inline G4double GetDst() const{ return Dst ;}
	 inline G4double GetImfx() const{ return Imfy ;}
	 inline G4double GetImfy() const{ return Imfy ;}
	 inline G4double GetImfz() const{ return Imfz ;}
	 inline G4double GetStdOff() const{ return StdOff ;}
	 inline G4double GetG1() const{ return G1 ;}
	 inline G4double GetG2() const{ return G2 ;}
	 inline G4double GetNsw() const{ return Nsw ;}
	 inline G4double GetVsw() const{ return Vsw ;}
	 inline G4double Getiopt() const{ return iopt;}
	 inline G4int Getnm_igrf() const{ return nm_igrf;}
#ifdef  USE_TSY04
	 inline G4double GetW1() const{ return W1 ;}
	 inline G4double GetW2() const{ return W2 ;}
	 inline G4double GetW3() const{ return W3 ;}
	 inline G4double GetW4() const{ return W4 ;}
	 inline G4double GetW5() const{ return W5 ;}
	 inline G4double GetW6() const{ return W6 ;}
	 	 
#endif	 											
	 inline void SetConsiderDipoleShift(G4bool abool)
	                               {ConsiderDipoleShift = abool;}
				       
#ifdef  USE_UNILIB
	void SelectUnilibModel(UL_Long kopt, UL_Long kext);
	inline void SetUseUnilibModel(bool aBool){use_unilib_models = aBool;}
	inline bool GetUseUnilibModel() const {return use_unilib_models;}
	inline UL_Long Getul_ifail() const {return ul_ifail;}
	inline double* Getul_param()const {return ul_param;}
	inline double Getul_year() const {return ul_year;}
	inline double Getul_mjd() const {return ul_mjd;}
	inline UL_Long Getul_kopt() const {return ul_kopt;}
	inline UL_Long Getul_kext()const {return ul_kext;}
#endif	  			       
#ifdef  USE_PALEO
	void  SetPaleoYear(double t);
#endif				       				
	 //McIlwain L parameter calculation
	 std::vector< double> ComputeMcIlwainLParameter(double pos[],G4double pitch_angle) ; 
	 std::vector< double> ComputeMcIlwainLParameter(G4ThreeVector pos,G4double pitch_angle); 
	 
protected:

private:

        // messenger
        MAGCOSFieldMessenger* theFieldMessenger;

        //Magnetic field model parameters  
        G4double TimeOfB;
	G4double year_igrf;
	DateAndTime StartDate,ReferenceDate;
	G4bool External,Internal,ConsiderDipoleShift; 
	G4int iopt;
	G4String IntFieldModelName, ExtFieldModelName;
	G4double DipoleB0,DipoleTheta,DipolePhi,DipolePS;
	G4ThreeVector DipoleSchift;
	G4double *times_of_data, *Pdyn_data, *Dst_data, *Imf_gsm_y_data, 
	         *Imf_gsm_z_data,
	         *g1_tsy_param_data,*g2_tsy_param_data;
        G4int n_tsy_data,nm_igrf;
	G4double Pdyn,Dst,Imfx, Imfy,Imfz,G1,G2,Vsw,Nsw,StdOff;
#ifdef  USE_TSY04
	G4double  W1,W2,W3,W4,W5,W6;
#endif

	//Igrf coefficient table
	std::vector<float>  igrf_year;
	std::vector< std::vector<float > > h_nm;
	std::vector< std::vector<float > > g_nm;
	std::vector<float>  dh_nm;
	std::vector<float>  dg_nm;
	
	//coefficient for the recurrence equation of legendre assocated function
	std::vector<float> recurrence_coef;
	G4int nmax;
	
	//Igrf coefficient according to reference time 
	std::vector<float> hh_nm;
	std::vector<float> gg_nm;
	
	std::vector<double> hh_nm1;
	std::vector<double> gg_nm1;
	
	float* hh;
	float* gg;
	float* rec;
	
	
	//equation of motion 
        MAGCOSEquationOfMotion* fEquationOfMotion; 
     
        //attribute for the integration method 
	G4ChordFinder* fChordFinder;
	G4MagIntegratorStepper* fStepper;
	G4double DefaultDeltaChord,DefaultBSMaxStep;
	G4double DeltaChord;
	G4double DefaultDeltaIntersection;
	G4double DefaultEpsilon;
	
	
	SecondInvariant* theSecondInvariantCalculator;

#ifdef  USE_UNILIB
	G4bool use_unilib_models;
	double* ul_param;
	double ul_year;
	double ul_mjd;
	UL_Long ul_kopt;
	UL_Long ul_kext;
	UL_Long ul_ifail;
	UL_Long ul_mfld;
#endif

//private methods       
private:     
        //IGRF model
	G4ThreeVector  GetIGRFFortran(G4ThreeVector pos) const; 
	G4ThreeVector  GetIGRFFortran1(G4ThreeVector pos) const; 
	G4ThreeVector  GetIGRF(G4ThreeVector pos) const; 
	G4ThreeVector  GetIGRF1(G4ThreeVector pos) const;
	G4ThreeVector  GetIGRFDOUBLE(G4ThreeVector pos) const;  
	
	//Geomagnetic dipole in GEO coordinate
	G4ThreeVector  GetGEODipole(G4ThreeVector pos) const;
	
	//Tsyganenko89 model 
	G4ThreeVector  GetTSY89(G4ThreeVector pos) const;
	G4ThreeVector  GetTSYBOB89(G4ThreeVector pos) const;
	
	//Tsyganenko96 model 
	G4ThreeVector  GetTSY96(G4ThreeVector pos) const;
	
	//Tsyganenko2001 model
	G4ThreeVector  GetTSY2001(G4ThreeVector pos) const;
#ifdef  USE_TSY04
	G4ThreeVector  GetTSY2004(G4ThreeVector pos) const;
#endif

#ifdef  USE_PALEO
	G4ThreeVector  GetCALS7K(G4ThreeVector pos) const;
#endif
	//Unilib models
#ifdef  USE_UNILIB	
	G4ThreeVector   GetUniLibBField(G4ThreeVector pos) const;
#endif	
	//pointers on the selected Internal and External field model
	G4ThreeVector (MAGCOSMagneticField::* GetInternalField)(G4ThreeVector) const;
	G4ThreeVector (MAGCOSMagneticField::* GetExternalField)(G4ThreeVector) const;
	
	
	//pointer on magnetopause model
	bool (MAGCOSMagneticField::*SelectedOutsideMagnetosphere)
	                                             (G4ThreeVector pos) const;
        
	
	//IGRF and geomagnetic magnetopause model r>25. re is outside
	//magnetosphere
	bool IGRFOutsideMagnetosphere(G4ThreeVector pos) const; 
	
	
	//tsy89 magnetopause model
	bool TSY89OutsideMagnetosphere(G4ThreeVector pos) const;
	
	//tsy89 magnetopause model
	bool TSY96OutsideMagnetosphere(G4ThreeVector pos) const;
        
	//tsy2001 magnetopause model
	bool TSY2001OutsideMagnetosphere(G4ThreeVector pos) const;
#ifdef  USE_TSY04
	bool TSY2004OutsideMagnetosphere(G4ThreeVector pos) const;
#endif	

	void ComputeTSYParameter(G4double t);
	   
} ;

#endif
