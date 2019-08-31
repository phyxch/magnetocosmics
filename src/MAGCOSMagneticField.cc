// 8/19/2014 Hexc, Olesya
// Checking this code out since the main program hangs on this object.
#include "MAGCOSMagneticField.hh"
#include "MAGCOSFieldMessenger.hh"
#include "globals.hh"
#include "geomdefs.hh"
#include "magneto_fsubroutine_def.hh"
#include "time.h"
#include "G4ios.hh"
#include "fstream"
#include "G4MagIntegratorStepper.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "SpaceCoordinateConvertor.hh"
#include "MAGCOSEquationOfMotion.hh"
#include "MAGCOSUnits.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ProcessTable.hh"
#include "BSIntegrator.hh"
#include "BSEquationInG4.hh"
#include <time.h>
#include "SecondInvariant.hh"

////////////////////////////////////////////////////////////////////////////////
//
MAGCOSMagneticField::MAGCOSMagneticField()
{ 
  //messenger
  theFieldMessenger = new MAGCOSFieldMessenger(this); 
  
  //second invarinat calculator
  
  theSecondInvariantCalculator = new SecondInvariant(this);
  
  //Parameters
  
  Pdyn=2.;
  Dst=0.;
  Imfx=0;
  Imfy=0.;
  Imfz=1.*nanotesla;
  G1=1.;
  G2=0.;
  nm_igrf=13;
  Nsw=5./cm3;
  Vsw=300.*km/s;
  StdOff=10.*re;
  
  times_of_data=0;
  Pdyn_data =0;
  Dst_data=0;
  Imf_gsm_y_data=0;
  Imf_gsm_z_data=0;
  g1_tsy_param_data=0;
  g2_tsy_param_data=0;
  n_tsy_data=0;
  
  
  //coefficients for igrf
  hh=new float[105];
  gg=new float[105];
  rec=new float[105];
 
  G4cout << " Checking the environment variable: " <<  getenv("IGRF_TABLE") << G4endl;
  //Igrf table
  if (getenv("IGRF_TABLE")){
  	ReadIgrfTable(getenv("IGRF_TABLE"));
  }
  else {
   G4cout<<"The environment variable IGRF_TABLE has not been defined"<<std::endl;
  }
    
  G4cout << "Done with reading the Table " << G4endl;
    
  //Set pointers on IGRF table for the SpaceCoordinateConvertor
  SpaceCoordinateConvertor* theConvertor = SpaceCoordinateConvertor::getInstance();
  theConvertor->Setp_hh_nm(&hh_nm);
  theConvertor->Setp_gg_nm(&gg_nm);
  
  // Set the start date and compute IGRF coefficient 
  TimeOfB=0.; 
  StartDate = DateAndTime(2000,1,1);
  ReferenceDate = StartDate;
  ComputeIgrfCoefAccordingToTime();
  theConvertor->SetReferenceDate(StartDate);
 
  
  //dipole parameter
  SetEccentricDipoleParameterFromIGRF();
  GetInternalField=&MAGCOSMagneticField::GetIGRF;
  GetExternalField=&MAGCOSMagneticField::GetTSY89;
  SelectedOutsideMagnetosphere=&MAGCOSMagneticField::IGRFOutsideMagnetosphere;
 
  External=false;
  Internal=true;
  iopt=1;
  
  
  
  //EquationOfMotion
  fEquationOfMotion = new MAGCOSEquationOfMotion(this);
   
  //Integration methods

  fChordFinder=0;
  fStepper=0;
  DefaultDeltaChord=.01*Re;
  DeltaChord=DefaultDeltaChord;
  DefaultDeltaIntersection= .001*re;
  DefaultEpsilon= 1.e-6;
  DefaultBSMaxStep=.5*re;
    
  G4TransportationManager::GetTransportationManager()->GetPropagatorInField()
	                        ->SetLargestAcceptableStep(10.*Re);
  G4TransportationManager::GetTransportationManager()->GetFieldManager()
                              ->SetDetectorField(this);
    
  ResetIntegrationParameters(); 
    
  ConsiderDipoleShift=false;

#ifdef  USE_UNILIB
  use_unilib_models =false;
  SelectUnilibModel(0,0);
  ul_ifail = -999999;
#endif  
  
}

////////////////////////////////////////////////////////////////////////
//
MAGCOSMagneticField::~MAGCOSMagneticField()
{ if (fEquationOfMotion) delete fEquationOfMotion;
  if (fChordFinder) delete fChordFinder;
  if (fStepper) delete fStepper;
}

///////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::GetFieldValue( const G4double yTrack[7],
                                          G4double B[3]         ) const 
{ G4ThreeVector pos=G4ThreeVector(yTrack[0],yTrack[1],yTrack[2]);
  G4ThreeVector Bfield =GetFieldValue(pos);
  B[0]=Bfield.x();
  B[1]=Bfield.y();
  B[2]=Bfield.z();    
  
  return;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetFieldValue
                               (const G4ThreeVector GEOposition) const 
{
 G4ThreeVector Bfield =G4ThreeVector(0.,0.,0.);

#ifdef  USE_UNILIB
  if (use_unilib_models && ul_ifail>= 0){
 	Bfield =  GetUniLibBField(GEOposition);
  }
  else {
#endif  
  	if (Internal) Bfield = (this->*GetInternalField)(GEOposition);
  	if (External) Bfield=Bfield+ (this->*GetExternalField)(GEOposition);
#ifdef  USE_UNILIB
  }
#endif   
  return Bfield;
}			       
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::FindPositionOnMagnetopause(G4double theta, 
                                                              G4double xgsm, 
							      G4double precision) const
{ G4ThreeVector pos1 = G4ThreeVector(xgsm,0.,0.);
 
  SpaceCoordinateConvertor* theCoordinateConvertor
                            = SpaceCoordinateConvertor::getInstance();
  theCoordinateConvertor->SetSystemInAndOut("GSM","GEO");
 
  if  (OutsideMagnetosphere(theCoordinateConvertor->Transform(pos1)))
                                            return G4ThreeVector (0.,0.,0.);
 
  G4ThreeVector pos2 =G4ThreeVector(xgsm,40.*re*std::sin(theta),40.*re*std::cos(theta));
  G4ThreeVector pos3;  
  while ( (pos2-pos1).mag()> precision){
 	pos3= (pos1+pos2)/2.;
	if (OutsideMagnetosphere(theCoordinateConvertor->Transform(pos3))) pos2=pos3;
	else pos1=pos3;
  }
  return (pos1+pos2)/2.;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::FindStandOffPosition(G4double precision ) const
{ G4ThreeVector pos1 = G4ThreeVector(0.,0.,0.)*Re;
  G4ThreeVector pos2 = G4ThreeVector(40.,0.,0.)*Re;
  G4ThreeVector pos3;  
  while ( (pos2-pos1).mag()> precision){
  	pos3= (pos1+pos2)/2.;
	if (OutsideMagnetosphere(pos3)) pos2=pos3;
	else pos1=pos3;
  }
  return pos1;	
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetTimeOfB(G4double val)
{ 
  TimeOfB=val;
  G4double JulDate = StartDate.JulianDate() + ( val / (24. * 3600.) );
  /* G4cout<<StartDate.min<<std::endl;
     G4cout<<StartDate.sec<<std::endl;*/
  ReferenceDate = DateAndTime(JulDate);
  // G4cout<<ReferenceDate.month<<" Reference date"<<std::endl;
  /*G4cout<<ReferenceDate.min<<std::endl;
    G4cout<<ReferenceDate.sec<<std::endl;*/
  ComputeTSYParameter(val);

  // new Igrf Coeff 
  ComputeIgrfCoefAccordingToTime();
  
  // set Reference date for Space coordinate convertor 
  SpaceCoordinateConvertor::getInstance()->SetReferenceDate(ReferenceDate);
  
  // new dipole parameters 
  SetEccentricDipoleParameterFromIGRF();
}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetStartDate(G4int year, G4int month,G4int day,
	                               G4int hour, G4int minute,G4int second )
{
  // Set the start date
  StartDate = DateAndTime(year,month,day,hour,minute,second);
  double val =TimeOfB;
  SetTimeOfB(val);
}		
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetIopt(  G4int val)
{
  iopt=val;
  if (val <1) iopt=1;
  if (val >7) iopt=7;
  return ;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetInternalField(  G4String aString)
{ Internal=true;
  IntFieldModelName =aString;
  if (aString == "DIPOLE") 
  GetInternalField=&MAGCOSMagneticField::GetGEODipole;
 /* else if (aString == "IGRFC++")
  	GetInternalField=&MAGCOSMagneticField::GetIGRF;
  else if (aString == "IGRFC++DOUBLE")
  	GetInternalField=&MAGCOSMagneticField::GetIGRFDOUBLE;	 
  else if (aString == "IGRFFORTRAN")
  	GetInternalField=&MAGCOSMagneticField::GetIGRFFortran;
  else if (aString == "IGRFFORTRAN1")
  	GetInternalField=&MAGCOSMagneticField::GetIGRFFortran1;*/
  else if (aString == "IGRF")
  	GetInternalField=&MAGCOSMagneticField::GetIGRFFortran1;	
#ifdef  USE_PALEO
  else if (aString == "CALS7K")
   	GetInternalField=&MAGCOSMagneticField::GetCALS7K;	
#endif
  else if (aString == "NOFIELD")
  	Internal=false;  
  	   
}
///////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetExternalField(  G4String aString)
{External=false;
 ExtFieldModelName =aString;
 if (aString == "TSY89"){
  	External=true;
  	GetExternalField=&MAGCOSMagneticField::GetTSY89;
  	SetMagnetopauseModel("TSY89");
 }
 else if (aString == "TSYBOB89"){ 
  	External=true;
  	GetExternalField=&MAGCOSMagneticField::GetTSYBOB89;
  	SetMagnetopauseModel("TSY89");
 } 
 else if (aString == "TSY96"){
	External=true;
  	GetExternalField=&MAGCOSMagneticField::GetTSY96;
  	SetMagnetopauseModel("TSY96");
 }
 else if (aString == "TSY2001"){
	External=true;
  	GetExternalField=&MAGCOSMagneticField::GetTSY2001;
  	SetMagnetopauseModel("TSY2001");
 }
#ifdef  USE_TSY04
  else if (aString == "TSY2004"){
	External=true;
  	GetExternalField=&MAGCOSMagneticField::GetTSY2004;
  	SetMagnetopauseModel("TSY2004");
 }
#endif 
 else  SetMagnetopauseModel("NOFIELD");
  
}
////////////////////////////////////////////////////////////////////////////////
//
void  MAGCOSMagneticField::SetMagnetopauseModel( G4String aString) 
{// G4cout<<aString<<std::endl; 
 if (aString == "TSY2001")
 	SelectedOutsideMagnetosphere
                  =&MAGCOSMagneticField::TSY2001OutsideMagnetosphere;
		  
#ifdef  USE_TSY04 
 else if (aString == "TSY2004")
 	SelectedOutsideMagnetosphere
                  =&MAGCOSMagneticField::TSY2004OutsideMagnetosphere;
#endif
 else if (aString == "TSY96")
     	SelectedOutsideMagnetosphere
                  =&MAGCOSMagneticField::TSY96OutsideMagnetosphere; 		  
		  
 else  if (aString == "TSY89")
     	SelectedOutsideMagnetosphere
                  =&MAGCOSMagneticField::TSY89OutsideMagnetosphere; 
	    
 else	
 	SelectedOutsideMagnetosphere
                  =&MAGCOSMagneticField::IGRFOutsideMagnetosphere; 	    

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetDipoleB0(G4double aVal)
{ DipoleB0=aVal;
  bdip_.b0dip=aVal/nanotesla;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetTiltedDipoleParameterFromIGRF()
{
  SpaceCoordinateConvertor* theConvertor = 
                           SpaceCoordinateConvertor::getInstance();
			   
  DipoleB0=theConvertor->GetGeoDipoleB0();
  G4ThreeVector DipoleAxis=theConvertor->GetGeoDipoleAxisInGEO();
  DipoleTheta=DipoleAxis.theta();
  DipolePhi=DipoleAxis.phi();
  bdip_.b0dip=DipoleB0/nanotesla;
  DipolePS=theConvertor->GetTiltAngle();
  DipoleSchift=G4ThreeVector(0.,0.,0.);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetEccentricDipoleParameterFromIGRF()
{ SpaceCoordinateConvertor* theConvertor = 
                           SpaceCoordinateConvertor::getInstance();
  DipoleB0=theConvertor->GetGeoDipoleB0();
  G4ThreeVector DipoleAxis=theConvertor->GetGeoDipoleAxisInGEO();
  DipoleTheta=DipoleAxis.theta();
  DipolePhi=DipoleAxis.phi();
  bdip_.b0dip=DipoleB0/nanotesla;
  DipolePS=theConvertor->GetTiltAngle();
  DipoleSchift=theConvertor->GetGeoDipoleSchiftInGEO();
 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetEpsilon(G4double aVal)
{ //precision for G4Integration methods
  G4TransportationManager::GetTransportationManager()
                                  ->GetPropagatorInField()
                                  ->SetMinimumEpsilonStep(aVal);
  G4TransportationManager::GetTransportationManager()
                                  ->GetPropagatorInField()
                                  ->SetMaximumEpsilonStep(aVal*1.00001);
  //precision for Bulirsh Stoer method
  BSIntegrator::getInstance()->SetEps(aVal);
 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetDeltaChord(G4double aVal)
{ if (fChordFinder) fChordFinder->SetDeltaChord(aVal);
  DeltaChord=aVal;
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetBSMaxStep(G4double aVal)
{ BSIntegrator::getInstance()->SetHmax(aVal);
}

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetDeltaIntersection(G4double aVal)
{ ///////Delta intersection for G4Integration method
  G4TransportationManager::GetTransportationManager()->GetFieldManager()
                                                 ->SetDeltaIntersection(aVal);

 ///////Delta intersection Bulirsh method
 BSEquationInG4* theBSEquation = 
           ((BSEquationInG4*) BSIntegrator::getInstance()
                                ->GetEquationTobeIntegrated());
 if (theBSEquation) theBSEquation->SetCrossingDelta(aVal);

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::ResetIntegrationParameters()
{ SetStepper("CashKarpRKF45");
  SetDeltaChord(DefaultDeltaChord);
  SetDeltaIntersection(DefaultDeltaIntersection);
  SetEpsilon(DefaultEpsilon);
  SetBSMaxStep(DefaultBSMaxStep);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SelectBulirshStoerMethod()
{ G4ProcessTable::GetProcessTable()->SetProcessActivation("MYTransportation",true);
  G4ProcessTable::GetProcessTable()->SetProcessActivation("Transportation",false);
 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SelectG4IntegrationMethods()
{ G4ProcessTable::GetProcessTable()->SetProcessActivation("MYTransportation",false);
  G4ProcessTable::GetProcessTable()->SetProcessActivation("Transportation",true);
 
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::ReadTSYParameter(G4String nameFile)
{
  if (times_of_data) delete times_of_data;
  if (Pdyn_data) delete Pdyn_data;
  if ( Dst_data) delete  Dst_data;
  if ( Imf_gsm_y_data) delete  Imf_gsm_y_data;
  if ( Imf_gsm_z_data) delete  Imf_gsm_z_data;
  if ( g1_tsy_param_data) delete  g1_tsy_param_data;
  if ( g2_tsy_param_data) delete  g2_tsy_param_data;
 
  n_tsy_data=0;
  std::fstream File_Input(nameFile, (std::ios::binary | std::ios::in));
  File_Input>>n_tsy_data;
 
  if (n_tsy_data <= 0){
  	G4cout<<"the file does not exist or is not well written"<<G4endl;
  } 
  else{ 
  	times_of_data= new G4double[n_tsy_data];
   	Pdyn_data= new G4double[n_tsy_data];
   	Dst_data= new G4double[n_tsy_data];
   	Imf_gsm_y_data= new G4double[n_tsy_data];
   	Imf_gsm_z_data= new G4double[n_tsy_data];
   	g1_tsy_param_data= new G4double[n_tsy_data];
   	g2_tsy_param_data= new G4double[n_tsy_data];
 
   	G4int nyear;
   	G4int nmonth;
   	G4int nday;
   	G4int nhour;
   	G4int nminute;
   	G4int nsecond=0;
   	File_Input>>nyear>>nmonth>>nday>>nhour>>nminute;
   	SetStartDate(nyear, nmonth, nday, nhour, nminute, nsecond); 
 
   	for (G4int i=0;i<n_tsy_data;i++){
   		File_Input>>times_of_data[i]
	      			>>Pdyn_data [i]
              			>> Dst_data[i]
              			>> Imf_gsm_y_data[i]
              			>> Imf_gsm_z_data[i]
              			>> g1_tsy_param_data[i]
              			>>g2_tsy_param_data[i];
   	}
  } 
  File_Input.close(); 
  double val =TimeOfB;
  SetTimeOfB(val);
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::ReadIgrfTable(G4String nameFile)
{ 
  n_tsy_data=0;
  std::fstream File_Input(nameFile, std::ios::in);
  
  //clear vectors
  igrf_year.clear();
  h_nm.clear();
  g_nm.clear();
  dh_nm.clear();
  dg_nm.clear();
  hh_nm.clear();
  gg_nm.clear();
  hh_nm1.clear();
  gg_nm1.clear();
  
  recurrence_coef.clear();
  
  //Process the first three lines 
  char Line[300];
  File_Input.getline(Line,300);
  G4cout << Line << G4endl;
  File_Input.getline(Line,300);
  G4cout << Line << G4endl;
  File_Input.getline(Line,300);
  G4cout << Line << G4endl;
  
  char ch;
  std::stringstream astream;
  File_Input.get(ch);
  while (ch != '\n') {
    astream<<ch;
    File_Input.get(ch);
  }

  G4String str1;
  
  astream>>str1>>str1>>str1;
  G4cout << " str1 = " << str1 << G4endl;

  if (str1!="m"){
    G4cout<<"The environment variable IGRF_TABLE is not defined correctly" <<std::endl;
    G4cout<<"The program will be interrupted"<<std::endl;
    exit(0);   
  }

  G4double year=0.;
  while (!astream.eof() && year >= 0.){
    year=-999999.;
    astream>>year;
    //G4cout << " year = " << year << G4endl;
    if (year!=-999999.){ 
      igrf_year.push_back(year);
    }
  }

  //number of column
  
  G4int nyear=igrf_year.size();   // "-2" is for skipping the "SV" column
  //G4cout << " I am here:  " << nyear << G4endl;

  nmax=0;
  G4int n,m;
  
  // initialisation h_nm,g_nm;  
  
  for (int i=0; i <nyear ; i++){
    h_nm.push_back(std::vector<float> ());
    	g_nm.push_back(std::vector<float> ());
  }
  G4cout << "h_nm and g_nm are initialized. " << G4endl;

  G4int nStop = 0;

  while (!File_Input.eof()){
    //if (nStop++ > 200) break;  // temp stuff
    str1="bibab";
    n=-999;
    m=-999;
    File_Input>>str1>>n>>m;
    G4cout<<str1<<'\t'<<n<<'\t'<<m<<'\t' <<G4endl;
    if (str1 != "bibab" && n != -999 && m != -999) {
      //G4cout<<str1<<'\t'<<n<<'\t'<<m<<'\t'<<G4endl;
      if (n>nmax){
	nmax=n;
	for (int i=0; i <nyear ; i++){
	  h_nm[i].insert(h_nm[i].end(),nmax+1,0.);
	  g_nm[i].insert(g_nm[i].end(),nmax+1,0.);
	}
	dg_nm.insert(dg_nm.end(),nmax+1,0.);
	dh_nm.insert(dh_nm.end(),nmax+1,0.);
	hh_nm.insert(hh_nm.end(),nmax+1,0.);
	gg_nm.insert(gg_nm.end(),nmax+1,0.);
	hh_nm1.insert(hh_nm1.end(),nmax+1,0.);
	gg_nm1.insert(gg_nm1.end(),nmax+1,0.);	    
      }
      
      G4int index = ((n+2)  * (n-1)) /2 + m;
      //G4cout<<"index "<<index<<std::endl;
      G4double value;
      int i_start=0;

      // // The following section was for reading old IGRF data table
      /*
      if (n >10) {
	i_start = 20;
	for (int i=0;i<i_start;i++){	
	  h_nm[i][index]  = 0.;
	  g_nm[i][index]  = 0.;
	}
      }		
      
      for (int i= i_start; i <nyear ; i++){
	File_Input>>value;
	//G4cout<<value<<'\t';
	if (str1 == "h") h_nm[i][index] = value;
	else if (str1 == "g") g_nm[i][index] = value;
      }
      if (n<9) {//cautiona the dh and dg parameter are pro year 
	File_Input>>value;
	//G4cout<<value<<'\t';
	if (str1 == "h") dh_nm[index] = value;
	else if (str1 == "g") dg_nm[index] = value;
	
      }
      else {
	dh_nm[index] = 0.;
	dg_nm[index] = 0.;
      }
      //G4cout<<std::endl;	
      */

      for (int i= i_start; i <nyear ; i++){
	File_Input>>value;
	//G4cout<<value<<'\t';
	if ( (str1 == "h") && (i < nyear-1) ) h_nm[i][index] = value;
	if ( (str1 == "g") && (i < nyear-1) ) g_nm[i][index] = value;
	if ( (str1 == "h") && (i == nyear-1) ) dh_nm[index] = value;
	if ( (str1 == "g") && (i == nyear-1) ) dg_nm[index] = value;     	
      }

    }
  }
  
  //  G4cout << " I am here " << year << G4endl;
  
 //Here we calculate the coefficient for the recurrence equation used to compute 
 //the associated Legendre function in the Gauss normalsation. 
 // In this normalisation the equation of recurrence is
 // P_n_m(x) = x * p_n-1_m -(n+m)*(n-m)* p_n-2_m / ((2n-3)*(2n-1)) (1)
 // The coefficient that we calculate here is (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1))
 // The qeaution (1) is a transformation of the classical  equation of recurrence of lengendre assocaited function
 // (n-m) * P_n_m = x  * p_n-1_m *(2n-1) -(n + m -1) * p_m_n-2  (2)
 // by using the Gauss normalisation (2) transform in (1)
 // (1) is better for computing purpose because you have just to compute 
 //the coefficient  (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) at initialisation an multiplied it by
 //   p_m_n-2  when computinhg p_n_m
 // In  case of (2) you have to compute three vector of cooefficient at initilaisation you have to
 // make three multiplications instead of one 


  for (n=1;n<nmax + 1;n++)
    for (int m= 0; m < n+1; m++) 
      recurrence_coef.push_back(
				float ((n-1+m)*(n-1-m)) / float ((2*n-1)*(2*n-3)) );
  
  for (n=1;n<=nmax + 1;n++)
    for (int m= 1; m < n+1; m++){
      int mn=n*(n-1)/2+m;
      igrfcc_.rec[mn-1] = float ((n-m)*(n+m-2)) / float ((2*n-1)*(2*n-3));
      
    };
   
 
  File_Input.close();
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::ComputeIgrfCoefAccordingToTime()
{
  G4int nyear = igrf_year.size();
  G4int n=(nmax+3)*nmax/2;
  
  // G4cout<<"number of year is "<<nyear<<endl;
  
  if (ReferenceDate.year < igrf_year[0]  || 
      ReferenceDate.year> igrf_year[nyear-1] + 5 ){
    G4cout<<"The year should be in the period "<<int(igrf_year[0])
	  <<"-"<<int(igrf_year[nyear-1]) + 5<<std::endl;
    
    G4cout<<"The coefficient will be set according to year "<<
      int(igrf_year[nyear-1])<<std::endl;
    for (int i=0;i<n;i++){
      gg_nm[i]=g_nm[nyear-1][i];
      hh_nm[i]=h_nm[nyear-1][i];
      gg_nm1[i]=g_nm[nyear-1][i];
      hh_nm1[i]=h_nm[nyear-1][i];}
  }
  else {
    G4int index=0;
    
    while (ReferenceDate.year >= (igrf_year)[index] && index < nyear) index++;
    
    // G4cout <<"Selected index"<<index<<std::endl;
    if (index<nyear){
      DateAndTime Date1=DateAndTime(int (igrf_year[index-1]), 1,1);
      DateAndTime Date2=DateAndTime(int (igrf_year[index]), 1,1);
      G4double delta_day12=Date2.DifferenceInDays(Date1);
      G4double delta_day=ReferenceDate.DifferenceInDays(Date1);
      G4double f2= delta_day / delta_day12;
      G4double f1=1.-f2;
      for (int i=0;i<n;i++){
	gg_nm[i]=g_nm[index-1][i]*f1 + g_nm[index][i]*f2;
	hh_nm[i]=h_nm[index-1][i]*f1 + h_nm[index][i]*f2;
	hh_nm1[i]=hh_nm[i];
	gg_nm1[i]=gg_nm[i];
      }
    }
    else {
      DateAndTime Date1=DateAndTime(int (igrf_year[index-1]), 1,1);
      G4double delta_day=ReferenceDate.DifferenceInDays(Date1);
      G4double dt=delta_day / 365.25;
      for (int i=0;i<n;i++){
	gg_nm[i]=g_nm[nyear-1][i] + dt*dg_nm[i];
	hh_nm[i]=h_nm[nyear-1][i] + dt*dh_nm[i];
	hh_nm1[i]=hh_nm[i];
	gg_nm1[i]=gg_nm[i];
      }
    }
  }
  
  //The coefficient are  Schmidt quasi normalised
  // For using the recuurence relation in Gauss normalisation
  // we have to convert the coefficient in Gauss normalisation
  // the g_nm and h_nm coefficient become
  float s=1.;
  for (n=1;n<nmax+1;n++){
    G4int mn = ((n + 1) * n) / 2  -1;
    s*=float(2*n-1)/float(n);
    hh_nm[mn] *= s;
    gg_nm[mn] *= s;
    G4double p=s;
    for (int m=1;m<n+1;m++){
      float aa=1.;
      if (m == 1) aa=2.;
      p *= std::sqrt(aa * float(n - m + 1) / float (n + m) );
      int nmn =  mn + m ;
      hh_nm[nmn] *= p;
      gg_nm[nmn] *= p;
      hh_nm1[nmn]=hh_nm[nmn];
      gg_nm1[nmn]=gg_nm[nmn];
      
    }
  }
  igrfcc_.hh[0]=0.;
  igrfcc_.gg[0]=0.;
  G4cout << "Need to check on the hh and gg array dimension !!! " << hh_nm.size() << G4endl;
  for (unsigned int i=0; i<hh_nm.size(); i++){
    hh[i]=hh_nm[i];
    gg[i]=gg_nm[i];
    igrfcc_.hh[i+1]=hh_nm[i];
    igrfcc_.gg[i+1]=gg_nm[i];
  } 				   			     
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::PrintBfield(G4ThreeVector geo_pos) const
{ 
  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  
  G4ThreeVector Bfield =G4ThreeVector(0.,0.,0.);
  if (Internal) Bfield = (this->*GetInternalField)(geo_pos);
  if (External) Bfield=Bfield+ (this->*GetExternalField)(geo_pos);
  Bfield /= nanotesla;
  G4ThreeVector Bgsm = theConvertor->Transform(Bfield, "GEO", "GSM");
  G4ThreeVector Bmag = theConvertor->Transform(Bfield, "GEO", "MAG");
  G4ThreeVector Bsm = theConvertor->Transform(Bfield, "GEO", "SM");
  
  G4cout<<"Bfield at primary position" <<std::endl; 		   
  G4cout<<"Bfield GEO: "<<Bfield<<" "<<Bfield.mag()<<std::endl;
  G4double altitude, longitude, latitude;
  theConvertor->ComputeGEOIDCoordinatesFromGEOPosition(geo_pos,
						       altitude,longitude,latitude);
  
  G4cout<<"longitude "<<longitude/degree<<std::endl;
  //Northward, eastward, dowmward component
  //---------------------------------------
  
  G4ThreeVector vertical = G4ThreeVector(1.,0.,0.);
  vertical.setRThetaPhi(1.,90.*degree-latitude,longitude);
  G4ThreeVector east_dir =G4ThreeVector(-std::sin(longitude),std::cos(longitude),0. );
  G4ThreeVector north_dir = -east_dir.cross(vertical);
  G4double Bdownward = -Bfield.dot(vertical);
  G4double Bnorthward = Bfield.dot(north_dir);
  G4double Beastward = Bfield.dot(east_dir);
  G4double h=std::sqrt(Bnorthward*Bnorthward+Beastward*Beastward);
  G4double inclination=
    std::acos(h/Bfield.mag());
  if (	Bdownward<0) 	inclination=-inclination;		
  G4double declination=
    std::acos(Bnorthward/h);
  if (Beastward <0) declination=-declination;			
  G4cout<<"Bfield north, east, down: "<<Bnorthward
	<<" "<<Beastward<<" "<<Bdownward<<std::endl;
  G4cout<<"Inclination, declination [degree] : "<<inclination/degree
	<<" "<<declination/degree<<std::endl;		       
  G4cout<<"Bfield GSM: "<<Bgsm<<" "<<Bgsm.mag()<<std::endl;
  G4cout<<"Bfield SM: "<<Bsm<<" "<<Bsm.mag()<<std::endl;
  G4cout<<"Bfield MAG: "<<Bmag<<" "<<Bmag.mag()<<std::endl; 
}
////////////////////////////////////////////////////////////////////////////////
//
#ifdef  USE_UNILIB
void MAGCOSMagneticField::SelectUnilibModel(UL_Long kopt, UL_Long kext)
{
  DateAndTime date1= DateAndTime(ReferenceDate.year, 1, 1);
  DateAndTime date2= DateAndTime(ReferenceDate.year+1, 1, 1); 
  ul_year = ReferenceDate.year + 
    (ReferenceDate.JulianDate()-date1.JulianDate())/
    (date2.JulianDate()-date1.JulianDate());
  UL_Long year=ReferenceDate.year;
  UL_Long month=ReferenceDate.month;
  UL_Long day=ReferenceDate.day;
  UL_Long hour=ReferenceDate.hour;
  UL_Long min=ReferenceDate.min;
  double sec=double (ReferenceDate.sec);
  
  SelectedOutsideMagnetosphere
    =&MAGCOSMagneticField::IGRFOutsideMagnetosphere; 
  
  
  
  ul_ifail= ulcalend( &ul_mjd, &year, 
		      &month, 
		      &day, 
		      &hour, 
		      &min, &sec,UL_COMPUTE_EPOCH); 
  ul_param = new double[15];
  ul_param[0]=iopt;
  ul_param[1]=Dst/nanotesla;
  ul_param[2]=Pdyn;
  ul_param[3]=Nsw*cm3;
  ul_param[4]=Vsw*s/km;
  ul_param[5]=Imfx/nanotesla;
  ul_param[6]=Imfy/nanotesla;
  ul_param[7]=Imfz/nanotesla;
  ul_param[8]=StdOff/re;
  ul_param[9]=DipolePS/degree;
  ul_param[10]=0.;
  ul_param[11]=0.;
  ul_param[12]=DipoleTheta/degree;
  ul_param[13]=DipolePhi/degree;
  ul_param[14]=DipoleB0/gauss;	
  ul_kopt =kopt;
  ul_kext = kext;
  ul_ifail= ulsetmfld( &ul_mfld, ul_kext, ul_year, ul_param, ul_kopt);
  if (ul_ifail<0){
    G4cout<<"The selection of a Unilib magnetic field model failed!"<<std::endl;
    G4cout<<"Therefore no magnetic field model form Unilib will be used."<<std::endl; 
    G4cout<<"The program returned the tag error ifail = "<<ul_ifail<<std::endl;
    G4cout<<"You should look into the Unilib documentation for the meaning of this error."<<std::endl;
  }
  
}
#endif

////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::PrintStormParameter()
{ 
  G4cout<<"iopt "<<iopt<<G4endl;
  G4cout<<"Pdyn "<<Pdyn<<G4endl;
  G4cout<<"Dst "<<Dst/nanotesla<<G4endl;
  G4cout<<"Imfy "<<Imfy/nanotesla<<G4endl;
  G4cout<<"Imfz "<<Imfz/nanotesla<<G4endl;
  G4cout<<"G1 "<<G1<<G4endl;
  G4cout<<"G2 "<<G2<<G4endl;
  G4cout<<"PS "<<DipolePS/degree<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetIGRFFortran(G4ThreeVector pos) const 
{ 
  G4double r=pos.mag()/re;
  G4double t=pos.theta();
  G4double p=pos.phi();
  float Br=1.;
  float Bp=1.;
  float Bt=1.;
  G4int nm=nm_igrf;
  //year_igrf=1980;
  
  /*G4double ndays_per_y  = DateAndTime(ReferenceDate.year+1,1,1)
                                  .DifferenceInDays(DateAndTime(ReferenceDate.year,1,1));
				  
  G4double fraction_of_year =
                     ReferenceDate.DifferenceInDays(
		     DateAndTime(ReferenceDate.year,1,1))/ndays_per_y;*/
  // float yea=float(ReferenceDate.year + fraction_of_year)  ;
  float r1 =float(r);
  float t1 =float(t);
  float p1 =float(p);
  
  if (nm > 10) nm=10;
  
  // igrf1_(&yea,&nm,&r1,&t1,&p1,&Br,&Bt,&Bp);
  
  igrf_to_cc__(&r1,&t1,&p1,&Br,&Bt,&Bp);
  
  G4double Bunit=nanotesla;
  G4double st=sin(t);
  G4double ct=cos(t);
  G4double sp=sin(p);
  G4double cp=cos(p);
  G4double Be=double(Br)*st+double(Bt)*ct;
  G4double Bx=Bunit*(Be*cp-double(Bp)*sp);
  G4double By=Bunit*(Be*sp+double(Bp)*cp);
  G4double Bz=Bunit*(double(Br)*ct-double(Bt)*st);
  if (r < 0.01){   
    Bx=.00000001*tesla;
    By=0.;
    Bz=0.;
  }
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetIGRFFortran1(G4ThreeVector pos) const 
{ 
  G4double r=pos.mag()/re;
  G4double t=pos.theta();
  G4double p=pos.phi();
  float Br=1.;
  float Bp=1.;
  float Bt=1.;
  G4int old_nm_igrf =nm_igrf;
  
  if (ReferenceDate.year <2000. && nm_igrf >10) old_nm_igrf =10; 
  igrfcc_.nm=old_nm_igrf;
  
  float r1 =float(r);
  float t1 =float(t);
  float p1 =float(p);
  
  igrf_to_cc__(&r1,&t1,&p1,&Br,&Bt,&Bp);
  
  G4double Bunit=nanotesla;
  G4double st=sin(t);
  G4double ct=cos(t);
  G4double sp=sin(p);
  G4double cp=cos(p);
  G4double Be=double(Br)*st+double(Bt)*ct;
  G4double Bx=Bunit*(Be*cp-double(Bp)*sp);
  G4double By=Bunit*(Be*sp+double(Bp)*cp);
  G4double Bz=Bunit*(double(Br)*ct-double(Bt)*st);
  if (r < 0.01){   
    Bx=.00000001*tesla;
    By=0.;
    Bz=0.;
  }
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetIGRF(G4ThreeVector pos) const 
{ 
  // This method computes the igrf magnetic field
  // For this purpose associated legendre function are calculated by the recurrence equation
  // P_n_m(x) = x * p_n-1_m(x) - recurrence_coef(n,m) * p_n-2_m(x)  (1)
  // where recurrence_coef = (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) and x= cos_theta
  // By deriving  equation (1) by theta we obtain
  // dp_n_m/dtheta = -sin_theta * p_n-1_m + cos_theta * dp_n-1_m/dtheta  -  
  //                                                        recurrence_coef * dp_n-2_m/dtheta
  // and dp_n_m/dtheta can  also be computed by recurrence. 
  // we should also now p_m_m = (sin (theta))^m
  
  float r=float(pos.mag()/re);
  float theta=float(pos.theta());
  if (theta < 0.01*degree) theta =0.01*degree;
  else if (theta > 179.99*degree) theta =179.99*degree;
  float phi=float(pos.phi());
  
  float pp = float(1.) / r;
  float p =  pp * pp;
    
  float cos_phi = std::cos(phi);
  float two_cos_phi = std::cos(phi) * 2.;
  float sin_phi = std::sin(phi);
  float cos_theta = std::cos(theta);
  float sin_theta = std::sin(theta);
   
   // a_n = (1 /r) ^ (n+1)
  // da_n_dth = -(a_n1) * n+1 
 
  std::vector<float> a;
  std::vector<float> da_dr;
  G4int n;
  G4int old_nm_igrf =nm_igrf;
  if (ReferenceDate.year <2000. && nm_igrf >10) old_nm_igrf =10; 	
    
  for ( n=1;n<old_nm_igrf+1;n++){
    	p *= pp;
     	a.push_back(p);
     	da_dr.push_back( -p * float (n+1) );
  }
    
  // case where m=0, p_1_0 =1.
  //---------------------------
    
  // n=1 
  float p_nm=cos_theta;
  float dp_nm_dth=-sin_theta;
    
  float Br=-da_dr[0]*gg_nm[0]*p_nm;
  float Bt=-a[0]*gg_nm[0]*dp_nm_dth;
  float Bp=0.;
    
  float p_nm1=p_nm; //corresponds to P for n-1 and m
  float p_nm2=1.;   //corresponds to P for n-2 and m
  
  float dp_nm1_dth=dp_nm_dth;
  float dp_nm2_dth=0.;
    
  // n>2
    
 for (n = 2 ; n < old_nm_igrf +1 ; n++ ){
   G4int ii= ((n+2) * (n-1) )/2; 
   p_nm = cos_theta *p_nm1 - recurrence_coef[ii]*p_nm2;
   dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
     - recurrence_coef[ii]*dp_nm2_dth;
   Br-=da_dr[n-1]*gg_nm[ii]*p_nm;
   Bt-=a[n-1]*gg_nm[ii]*dp_nm_dth; 
   p_nm2=p_nm1;
   p_nm1=p_nm;
   dp_nm2_dth=dp_nm1_dth;
   dp_nm1_dth=dp_nm_dth;
 } 
 
 // continue for m>0
 //---------  
 
 float cos_mphi2=cos_phi;
 float cos_mphi1=1.;
 float cos_mphi=cos_phi;
 float sin_mphi2=-sin_phi;
 float sin_mphi1=0.;
 float sin_mphi=sin_phi;
 
 
 float p_mm=1.;
 float dp_mm_dth=1.; 
 float Bpm=0.;
 float Bp_hh,Bp_gg,Br_hh,Br_gg,Bt_hh,Bt_gg;
 for (int m = 1 ; m < old_nm_igrf +1 ; m++ ){
   //Initialisation
   //------------------
   
   Bp_hh=0.;
   Bp_gg=0.;
   Br_hh=0.;
   Br_gg=0.;
   Bt_hh=0.;
   Bt_gg=0.;
   
   //compute cos_mphi and sin_mphi with
   // recurrence formula from numerical recipies in C++ p 184 
   cos_mphi= two_cos_phi*cos_mphi1 -cos_mphi2;
   sin_mphi= two_cos_phi*sin_mphi1 -sin_mphi2;
   
   
   cos_mphi2=cos_mphi1;
   sin_mphi2=sin_mphi1;
   cos_mphi1=cos_mphi;
   sin_mphi1=sin_mphi;
   
   dp_mm_dth= float(m) * p_mm *cos_theta;
   p_mm *= sin_theta;
   Bpm=0;
   // case where n=m
   n=m;
   G4int ii= ((n+3) * (n) )/2 -1; 
   p_nm = p_mm;
   dp_nm_dth = dp_mm_dth;
   
   float gg=gg_nm[ii];
   float hh=hh_nm[ii];
   float aa=a[n-1];
   float da=da_dr[n-1]; 
   
   float ggg = gg*p_nm;
   float hhh = hh*p_nm;
   
   Br_gg-=da*ggg;
   Br_hh-=da*hhh;
   
   float w=aa*dp_nm_dth;
   Bt_gg-=gg*w;
   Bt_hh-=hh*w;
   
   //Bt-=a[n-1]*gh*dp_nm_dth;
   
   //Bpm-=a[n-1]*dgh_dphi*p_nm; 
   
   Bp_gg-=aa*ggg;
   Bp_hh-=aa*hhh; 
   
   
   p_nm1=p_nm; //corresponds to P for n-1 and m
   p_nm2=0.;  //corresponds to P for n-2 and m
   dp_nm1_dth=dp_nm_dth;
   dp_nm2_dth=0.;
   
   for ( n = m + 1 ; n < old_nm_igrf +1 ; n++ ){
     G4int ii= ((n+2) * (n-1) )/2 + m;
     float rc= recurrence_coef[ii];
     p_nm = cos_theta *p_nm1 - rc*p_nm2;
     dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
       - rc*dp_nm2_dth;
     
     
     gg=gg_nm[ii];
     hh=hh_nm[ii];
     
     float ggg = gg*p_nm;
     float hhh = hh*p_nm;
     da=da_dr[n-1];
     aa=a[n-1];
     
     Br_gg-=da*ggg;
     Br_hh-=da*hhh;
     
     w=aa*dp_nm_dth;
     Bt_gg-=w*gg;
     Bt_hh-=w*hh;
     
     Bp_gg-=aa*ggg;
     Bp_hh-=aa*hhh; 
     
     p_nm2=p_nm1;
     p_nm1=p_nm;
     dp_nm2_dth=dp_nm1_dth;
     dp_nm1_dth=dp_nm_dth;
   }
   Br+=Br_gg*cos_mphi + Br_hh*sin_mphi;
   Bt+=Bt_gg*cos_mphi + Bt_hh*sin_mphi;
   Bp+=float(m)*(Bp_hh*cos_mphi - Bp_gg*sin_mphi);
   
 } 
 
 if (sin_theta >= 1.e-8) Bp/=sin_theta; 
 else if (cos_theta<0.) Bp=-Bp;
 
 G4double Be= Br*sin_theta+ Bt*cos_theta;
 G4double Bx= Be*cos_phi- Bp*sin_phi;
 G4double By= Be*sin_phi+ Bp*cos_phi;
 G4double Bz= Br*cos_theta - Bt*sin_theta;
 if (r < 0.01){   
   Bx=.00000001*tesla;
   By=0.;
   Bz=0.;
 }
 
 G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz)*nanotesla;
 return Bfield;
}

////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetIGRFDOUBLE(G4ThreeVector pos) const 
{ 
  // This method computes the igrf magnetic field
  // For this purpose associated legendre function are calculated by the recurrence equation
  // P_n_m(x) = x * p_n-1_m(x) - recurrence_coef(n,m) * p_n-2_m(x)  (1)
  // where recurrence_coef = (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) and x= cos_theta
  // By deriving  equation (1) by theta we obtain
  // dp_n_m/dtheta = -sin_theta * p_n-1_m + cos_theta * dp_n-1_m/dtheta  -  
  //                                                        recurrence_coef * dp_n-2_m/dtheta
  // and dp_n_m/dtheta can  also be computed by recurrence. 
  // we should also now p_m_m = (sin (theta))^m

  double r=double(pos.mag()/re);
  double theta=double(pos.theta());
  if (theta < 0.01*degree) theta =0.01*degree;
  else if (theta > 179.99*degree) theta =179.99*degree;
  double phi=double(pos.phi());
  
  double pp = double(1.) / r;
  double p =  pp * pp;
   
  double cos_phi = std::cos(phi);
  double two_cos_phi = std::cos(phi) * 2.;
  double sin_phi = std::sin(phi);
  double cos_theta = std::cos(theta);
  double sin_theta = std::sin(theta);
   
  // a_n = (1 /r) ^ (n+1)
  // da_n_dth = -(a_n1) * n+1 
 
  std::vector<double> a;
  std::vector<double> da_dr;
  G4int n;
  G4int old_nm_igrf =nm_igrf;
  if (ReferenceDate.year <2000. && nm_igrf >10) old_nm_igrf =10; 	
   
   
  for ( n=1;n<old_nm_igrf+1;n++){
    	p *= pp;
     	a.push_back(p);
     	da_dr.push_back( -p * double (n+1) );
  }
    
  // case where m=0, p_1_0 =1.
  //---------------------------
    
  // n=1 
  double p_nm=cos_theta;
  double dp_nm_dth=-sin_theta;
    
  double Br=-da_dr[0]*gg_nm1[0]*p_nm;
  double Bt=-a[0]*gg_nm1[0]*dp_nm_dth;
  double Bp=0.;
    
  double p_nm1=p_nm; //corresponds to P for n-1 and m
  double p_nm2=1.;   //corresponds to P for n-2 and m
  
  double dp_nm1_dth=dp_nm_dth;
  double dp_nm2_dth=0.;
    
  // n>2
    
 for (n = 2 ; n < old_nm_igrf +1 ; n++ ){
   G4int ii= ((n+2) * (n-1) )/2; 
   p_nm = cos_theta *p_nm1 - recurrence_coef[ii]*p_nm2;
   dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
     - recurrence_coef[ii]*dp_nm2_dth;
   Br-=da_dr[n-1]*gg_nm1[ii]*p_nm;
   Bt-=a[n-1]*gg_nm1[ii]*dp_nm_dth; 
   p_nm2=p_nm1;
   p_nm1=p_nm;
   dp_nm2_dth=dp_nm1_dth;
   dp_nm1_dth=dp_nm_dth;
 } 
 
 // continue for m>0
 //---------  
 
 double cos_mphi2=cos_phi;
 double cos_mphi1=1.;
 double cos_mphi=cos_phi;
 double sin_mphi2=-sin_phi;
 double sin_mphi1=0.;
 double sin_mphi=sin_phi;
 
 double p_mm=1.;
 double dp_mm_dth=1.; 
 double Bpm=0.;
 double Bp_hh,Bp_gg,Br_hh,Br_gg,Bt_hh,Bt_gg;
 for (int m = 1 ; m < old_nm_igrf +1 ; m++ ){
   //Initialisation
   //------------------
   
   Bp_hh=0.;
   Bp_gg=0.;
   Br_hh=0.;
   Br_gg=0.;
   Bt_hh=0.;
   Bt_gg=0.;
   
   //compute cos_mphi and sin_mphi with
   // recurrence formula from numerical recipies in C++ p 184 
   cos_mphi= two_cos_phi*cos_mphi1 -cos_mphi2;
   sin_mphi= two_cos_phi*sin_mphi1 -sin_mphi2;
   
   
   cos_mphi2=cos_mphi1;
   sin_mphi2=sin_mphi1;
   cos_mphi1=cos_mphi;
   sin_mphi1=sin_mphi;
   
   dp_mm_dth= double(m) * p_mm *cos_theta;
   p_mm *= sin_theta;
   Bpm=0;
   // case where n=m
   n=m;
   G4int ii= ((n+3) * (n) )/2 -1; 
   p_nm = p_mm;
   dp_nm_dth = dp_mm_dth;
   
   double gg=gg_nm1[ii];
   double hh=hh_nm1[ii];
   double aa=a[n-1];
   double da=da_dr[n-1]; 
   
   double ggg = gg*p_nm;
   double hhh = hh*p_nm;
   
   Br_gg-=da*ggg;
   Br_hh-=da*hhh;
   
   double w=aa*dp_nm_dth;
   Bt_gg-=gg*w;
   Bt_hh-=hh*w;
   
   //Bt-=a[n-1]*gh*dp_nm_dth;
   
   //Bpm-=a[n-1]*dgh_dphi*p_nm; 
   
   Bp_gg-=aa*ggg;
   Bp_hh-=aa*hhh; 
   
   
   p_nm1=p_nm; //corresponds to P for n-1 and m
   p_nm2=0.;  //corresponds to P for n-2 and m
   dp_nm1_dth=dp_nm_dth;
   dp_nm2_dth=0.;
   
   
   for ( n = m + 1 ; n < old_nm_igrf +1 ; n++ ){
     G4int ii= ((n+2) * (n-1) )/2 + m;
     double rc= recurrence_coef[ii];
     p_nm = cos_theta *p_nm1 - rc*p_nm2;
     dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
       - rc*dp_nm2_dth;
     
     
     gg=gg_nm1[ii];
     hh=hh_nm1[ii];
     
     
     double ggg = gg*p_nm;
     double hhh = hh*p_nm;
     da=da_dr[n-1];
     aa=a[n-1];
     
     Br_gg-=da*ggg;
     Br_hh-=da*hhh;
                
     w=aa*dp_nm_dth;
     Bt_gg-=w*gg;
     Bt_hh-=w*hh;
     
     Bp_gg-=aa*ggg;
     Bp_hh-=aa*hhh; 
     
     p_nm2=p_nm1;
     p_nm1=p_nm;
     dp_nm2_dth=dp_nm1_dth;
     dp_nm1_dth=dp_nm_dth;
   }
   Br+=Br_gg*cos_mphi + Br_hh*sin_mphi;
   Bt+=Bt_gg*cos_mphi + Bt_hh*sin_mphi;
   Bp+=double(m)*(Bp_hh*cos_mphi - Bp_gg*sin_mphi);
   
 } 
 
 if (sin_theta >= 1.e-8) Bp/=sin_theta; 
 else if (cos_theta<0.) Bp=-Bp;
 
 G4double Be= Br*sin_theta+ Bt*cos_theta;
 G4double Bx= Be*cos_phi- Bp*sin_phi;
 G4double By= Be*sin_phi+ Bp*cos_phi;
 G4double Bz= Br*cos_theta - Bt*sin_theta;
 if (r < 0.01){   
   Bx=.00000001*tesla;
   By=0.;
   Bz=0.;
 }
 
 G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz)*nanotesla;
 return Bfield;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetGEODipole(G4ThreeVector pos) const
{ G4double xgeo,ygeo,zgeo; 
  xgeo=pos.x()/re;
  ygeo=pos.y()/re;
  zgeo=pos.z()/re;
  
  G4ThreeVector SchiftedPosition = 
    G4ThreeVector(xgeo,ygeo,zgeo)-(DipoleSchift/re); 
  
  G4double sth=sin(DipoleTheta);
  G4double cth=cos(DipoleTheta);
  G4double sphi=sin(DipolePhi);
  G4double cphi=cos(DipolePhi);
  
  //coordinates in DipoleSystem
  G4double xd=(SchiftedPosition.x()*cphi+SchiftedPosition.y()*sphi)*cth
    -sth*SchiftedPosition.z();
  G4double yd=SchiftedPosition.y()*cphi-SchiftedPosition.x()*sphi;
  G4double zd=(SchiftedPosition.x()*cphi+SchiftedPosition.y()*sphi)*sth
    +cth*SchiftedPosition.z();
  
  
  //component in dipole system		
  G4double Bxd,Byd,Bzd;
  G4double r2=xd*xd +yd *yd +zd*zd; 
  G4double r5=std::pow(r2,2.5);
  G4double factor=DipoleB0/r5;
  Bxd=-3*factor*xd*zd;
  Byd=-3*factor*yd*zd;
  Bzd=factor*(r2-3*zd*zd);
  
  
  //back to Geo coordinate
  G4double Bx,By,Bz;
  Bx=(cth*cphi*Bxd - sphi*Byd  + sth*cphi*Bzd);
  By=(cth*sphi*Bxd + cphi*Byd  + sth*sphi*Bzd);
  Bz=(-sth*Bxd  + cth*Bzd);
  
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}
/////////////////////////////////////////////////////////////////////////////////
//
#ifdef  USE_PALEO
void  MAGCOSMagneticField::SetPaleoYear(double t)
{double year_paleo = t;
  if (getenv("CALS7K_TABLE")){
    calculate_cals7k_coef__(getenv("CALS7K_TABLE"),&year_paleo);
  }
  else {
    G4cout<<"The environment variable CALS7K_TABLE has not been defined"<<std::endl;
  } 
}
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector  MAGCOSMagneticField::GetCALS7K(G4ThreeVector pos) const
{  
  //G4cout<<pos<<std::endl;
  double rgeo= pos.mag()/km;
  //G4cout<<rgeo<<std::endl;
  double theta_geo= pos.theta();
  if (theta_geo < 0.01*degree) theta_geo =0.01*degree;
  else if (theta_geo > 179.99*degree) theta_geo =179.99*degree;
  double phi_geo= pos.phi();
  double Br,Bt,Bp;
  compute_cals7k_field_for_magcos__(&rgeo,&theta_geo,&phi_geo,&Br,&Bt,&Bp);
  
  G4double st=sin(theta_geo);
  G4double ct=cos(theta_geo);
  G4double sp=sin(phi_geo);
  G4double cp=cos(phi_geo);
  G4double Be= Br*st+ Bt *ct;
  G4double Bx= Be*cp- Bp*sp;
  G4double By= Be*sp+ Bp*cp;
  G4double Bz= Br*ct- Bt*st;
  //G4cout<<theta_geo<<'\t'<<phi_geo<<std::endl;
  
  G4ThreeVector Bgeo=
    G4ThreeVector(Bx,By,Bz);
  
  Bgeo *=nanotesla;
  return Bgeo;	     
  
}
#endif
////////////////////////////////////////////////////////////////////////////////
/////////////// External field models

G4ThreeVector MAGCOSMagneticField::GetTSY89( G4ThreeVector pos) const
{  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  
  G4ThreeVector pos_geo=pos;
  if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  
  float xgsm=float(pos_gsm.x()/re);
  float ygsm=float(pos_gsm.y()/re);
  float zgsm=float(pos_gsm.z()/re);
   
  float parmod[10]; 
  float Bxgsm,Bygsm,Bzgsm;
  float ps=float(DipolePS);
  G4int ipar=iopt;
  t89c__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
	&Bxgsm, &Bygsm, &Bzgsm );
  
  G4ThreeVector Bgsm=
    G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
  Bgsm *= nanotesla;
  
  //GSMtoGEO
  return theConvertor->TransformGSMinGEO(Bgsm);
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetTSYBOB89( G4ThreeVector pos) const
{  
  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  
  G4ThreeVector pos_geo=pos;
  if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  
  float xgsm=float(pos_gsm.x()/re);
  float ygsm=float(pos_gsm.y()/re);
  float zgsm=float(pos_gsm.z()/re);
   
  float parmod[10]; 
  parmod[1]=Dst;
  float Bxgsm,Bygsm,Bzgsm;
  float ps=float(DipolePS);
  G4int ipar=iopt;
  t89c_boberg__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
		&Bxgsm, &Bygsm, &Bzgsm );
  
  G4ThreeVector Bgsm=
    G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
  Bgsm *= nanotesla;
  
  //GSMtoGEO
  return theConvertor->TransformGSMinGEO(Bgsm);
}
////////////////////////////////////////////////////////////////////////////
G4ThreeVector MAGCOSMagneticField::GetTSY96(G4ThreeVector pos) const
{ 
  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  
  G4ThreeVector pos_geo=pos;
  if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  
  float xgsm=float(pos_gsm.x()/re);
  float ygsm=float(pos_gsm.y()/re);
  float zgsm=float(pos_gsm.z()/re);
  
  float parmod[10]; 
  parmod[0]=Pdyn;
  parmod[1]=Dst/nanotesla; 
  parmod[2]=Imfy/nanotesla;
  parmod[3]=Imfz/nanotesla;
  
  float Bxgsm,Bygsm,Bzgsm;
  float ps=float(DipolePS);
  G4int ipar=iopt;
  mcos_t96_01__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
		&Bxgsm, &Bygsm, &Bzgsm );
  /* t96_01__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
     &Bxgsm, &Bygsm, &Bzgsm );*/
  G4ThreeVector Bgsm=
    G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
  Bgsm *= nanotesla;
  
  //GSMtoGEO
  return theConvertor->TransformGSMinGEO(Bgsm);
  
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetTSY2001( G4ThreeVector pos) const
{ 
  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  
  
  G4ThreeVector pos_geo=pos;
  
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  
  float xgsm=float(pos_gsm.x()/re);
  float ygsm=float(pos_gsm.y()/re);
  float zgsm=float(pos_gsm.z()/re);
  
  // call tsy2001 fortran subroutine   
  
  float parmod[10]; 
  parmod[0]=Pdyn;
  parmod[1]=Dst/nanotesla; 
  parmod[2]=Imfy/nanotesla;
  parmod[3]=Imfz/nanotesla;
  parmod[4]=G1;
  parmod[5]=G2;
  
  float Bxgsm,Bygsm,Bzgsm;
  float ps=float(DipolePS);
  
  G4int ipar=iopt;
  t01_01__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
	   &Bxgsm, &Bygsm, &Bzgsm );
  
  G4ThreeVector Bgsm=
    G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
  Bgsm *= nanotesla;
  
  //GSMtoGEO
  return theConvertor->TransformGSMinGEO(Bgsm);
}
#ifdef  USE_TSY04
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetTSY2004( G4ThreeVector pos) const
{ 
  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  
  
  G4ThreeVector pos_geo=pos;
  
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  
  float xgsm=float(pos_gsm.x()/re);
  float ygsm=float(pos_gsm.y()/re);
  float zgsm=float(pos_gsm.z()/re);
  
  // call tsy2004 fortran subroutine   
  
  float parmod[10]; 
  parmod[0]=Pdyn;
  parmod[1]=Dst/nanotesla; 
  parmod[2]=Imfy/nanotesla;
  parmod[3]=Imfz/nanotesla;
  parmod[4]=W1;
  parmod[5]=W2;
  parmod[6]=W3;
  parmod[7]=W4;
  parmod[8]=W5;
  parmod[9]=W6;
  
  float Bxgsm,Bygsm,Bzgsm;
  float ps=float(DipolePS);
  
  G4int ipar=iopt;
  t04_s__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
	  &Bxgsm, &Bygsm, &Bzgsm );
  
  G4ThreeVector Bgsm=
    G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
  Bgsm *= nanotesla;
  
  //GSMtoGEO
  return theConvertor->TransformGSMinGEO(Bgsm);
}
#endif 
#ifdef USE_UNILIB
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector MAGCOSMagneticField::GetUniLibBField(G4ThreeVector pos) const
{
  G4ThreeVector Bfield(0.,0.,0);
  double ul_geo_pos[4];
  ul_geo_pos[0]=pos.x()/re;
  ul_geo_pos[1]=pos.y()/re;
  ul_geo_pos[2]=pos.z()/re;
  ul_geo_pos[3]=ul_mjd;
  double ul_geo_bfield[4];
  UL_Long ifail;
  ifail= ulgeomfl( ul_geo_pos, ul_mfld, ul_geo_bfield);
  if (ifail <=0){
    G4cout<<"The magnetic field computation with the unilib package failed!"<<std::endl;
    return Bfield;
  }
  Bfield=G4ThreeVector(ul_geo_bfield[0],
		       ul_geo_bfield[1],
		       ul_geo_bfield[2])*gauss;
  	
  return Bfield;
}
#endif
////////////////////////////////////////////////////////////////////////////////
//-------------------------------Outside magnetosphere
bool MAGCOSMagneticField::OutsideMagnetosphere( G4ThreeVector pos) const
{ return  (this->*SelectedOutsideMagnetosphere)(pos);
}
////////////////////////////////////////////////////////////////////////////////
//
bool MAGCOSMagneticField::IGRFOutsideMagnetosphere( G4ThreeVector pos)  const
{ 
  G4double r=pos.mag()/re;
  bool res=false;
  if (r > 25.) res =true;
  return res;
}
////////////////////////////////////////////////////////////////////////////////
//
bool MAGCOSMagneticField::TSY89OutsideMagnetosphere( G4ThreeVector pos) const
{
  //model from Kobel PhD "Zu die magnetosphärischen Effekten der kosmischen
  //  Strahlung"
  G4double a=-0.0545;
  G4double b_iopt[]={11.7,11.1,10.8,10.4,10.4,10.2,10.2};
  G4double f_iopt[]={20.,15.,6.67,10.,5.,6.,6.};
  
  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos)/re;
  G4double rho2 = pos_gsm.y()*pos_gsm.y() + pos_gsm.z()*pos_gsm.z();
  G4bool res = (rho2>900 || pos_gsm.x() <-60);
  if (!res){ 
    double  k_par = DipolePS/f_iopt[iopt-1];
    G4double cos_k=std::cos(k_par);
    G4double sin_k=std::sin(k_par);
    G4RotationMatrix mat = G4RotationMatrix(G4ThreeVector(cos_k,0,sin_k),
					    G4ThreeVector(0,1,0), 
					    G4ThreeVector(-sin_k,0,cos_k));
    G4ThreeVector pos_rot = mat*pos_gsm;
    G4double rho_rot= pos_rot.y()*pos_rot.y() + pos_rot.z()*pos_rot.z();
    res = (pos_rot.x() > a*rho_rot +b_iopt[iopt-1]); 
  }
  
  return res;
}
///////////////////////////////////////////////////////////////////////////////
//
bool MAGCOSMagneticField::TSY96OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko96 code
  SpaceCoordinateConvertor* theConvertor =
    SpaceCoordinateConvertor::getInstance(); 
  
  G4ThreeVector pos_geo=pos;
  
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  
  
  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  
  // magnetopause model
  G4double am0=70.;
  G4double s0=1.08;
  G4double x00=5.48;
  G4double dsig=.005;
  G4double xappa=std::pow(Pdyn/2.,0.14);
  G4double x0=x00/xappa;
  G4double am=am0/xappa;
  G4double rho2=ygsm*ygsm+zgsm*zgsm;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));
  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
}
///////////////////////////////////////////////////////////////////////////////
//
bool MAGCOSMagneticField::TSY2001OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko2001 code
  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 

  G4ThreeVector pos_geo=pos;
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  

  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  G4double r= pos_gsm.mag()/re;
  
  G4double xss=xgsm;
  G4double zss=zgsm;
  
  G4double rh0=8.94335;
  G4double dd;
  G4double rh2= -5.2;
  G4double sps = std::sin(DipolePS);
  //change of coordinate due to the tail warping
  do{ G4double xsold=xss;
      G4double zsold=zss;
      G4double rh=rh0+rh2*(zss/r)*(zss/r);
      G4double sinpsas=sps/std::pow(1.+std::pow(r/rh,3.),1./3.);
      G4double cospsas=std::sqrt(1. - sinpsas*sinpsas);
      zss=xgsm*sinpsas+zgsm*cospsas;
      xss=xgsm*cospsas-zgsm*sinpsas;
      dd=std::abs(xss-xsold)+std::abs(zss-zsold);
  } while (dd > 1.e-6);
  
  // magnetopause model
  G4double a0_a=34.586;
  G4double s0=1.196;
  G4double a0_x0=3.4397;
  G4double dsig=.003;
  G4double xappa=std::pow(Pdyn/2.,0.15790);
  G4double x0=a0_x0/xappa;
  G4double am=a0_a/xappa;
  G4double rho2=ygsm*ygsm+zss*zss;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));
  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
}
#ifdef  USE_TSY04
///////////////////////////////////////////////////////////////////////////////
//
bool MAGCOSMagneticField::TSY2004OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko2004 code
  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 

  G4ThreeVector pos_geo=pos;
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  

  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  G4double r= pos_gsm.mag()/re;
  
  G4double xss=xgsm;
  G4double zss=zgsm;
  
  G4double rh0=7.5;
  G4double dd;
  G4double rh2= -5.2;
  G4double sps = std::sin(DipolePS);
  //change of coordinate due to the tail warping
  do{ G4double xsold=xss;
      G4double zsold=zss;
      G4double rh=rh0+rh2*(zss/r)*(zss/r);
      G4double sinpsas=sps/std::pow(1.+std::pow(r/rh,3.),1./3.);
      G4double cospsas=std::sqrt(1. - sinpsas*sinpsas);
      zss=xgsm*sinpsas+zgsm*cospsas;
      xss=xgsm*cospsas-zgsm*sinpsas;
      dd=std::abs(xss-xsold)+std::abs(zss-zsold);
  } while (dd > 1.e-6);
  
  // magnetopause model
  G4double a0_a=34.586;
  G4double s0=1.196;
  G4double a0_x0=3.4397;
  G4double dsig=.005;
  G4double xappa=std::pow(Pdyn/2.,0.154269);
  G4double x0=a0_x0/xappa;
  G4double am=a0_a/xappa;
  G4double rho2=ygsm*ygsm+zss*zss;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));
  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
}
#endif 
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::ComputeTSYParameter(G4double t)
{ G4int index =-1;
  if (n_tsy_data>0){
  	for (G4int i=0;i<n_tsy_data-1;i++){
   		if ( (t >= times_of_data[i]) && (t <= times_of_data[i+1])){
       			index=i;
        		i=n_tsy_data-1;
		}
	}
  	if (index > -1){
    		G4double dx=(t-times_of_data[index])
                	/(times_of_data[index+1]-times_of_data[index]);
     		Pdyn=Pdyn_data[index] + dx * (Pdyn_data[index+1]-Pdyn_data[index]);
     		Dst=Dst_data[index] + dx * (Dst_data[index+1]-Dst_data[index]);
     		Dst*=nanotesla;
		Imfy=Imf_gsm_y_data[index] + dx * 
                      (Imf_gsm_y_data[index+1]-Imf_gsm_y_data[index]);
		Imfy*=nanotesla;      
     		Imfz=Imf_gsm_z_data[index] + dx * 
                      (Imf_gsm_z_data[index+1]-Imf_gsm_z_data[index]);
     		Imfz*=nanotesla;
		G1=g1_tsy_param_data[index] + dx * 
                       (g1_tsy_param_data[index+1]-g1_tsy_param_data[index]);
     		G2=g2_tsy_param_data[index] + dx * 
                      (g2_tsy_param_data[index+1]-g2_tsy_param_data[index]);     
  	}
  }  		
}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::SetStepper(G4String aString)
{ if (fStepper) delete fStepper;
   
  if (aString == "ExplicitEuler"){
    	fStepper = new G4ExplicitEuler( fEquationOfMotion );
     	G4cout<<"G4ExplicitEuler is called"<<G4endl;
  }
  else if  (aString == "ImplicitEuler"){
    	fStepper = new G4ImplicitEuler( fEquationOfMotion );
     	G4cout<<"G4ImplicitEuler is called"<<G4endl;
  }
  else if  (aString == "SimpleRunge"){
    	fStepper = new G4SimpleRunge( fEquationOfMotion );
     	G4cout<<"G4SimpleRunge is called"<<G4endl;
  } 
  else if  (aString == "ClassicalRK4"){
    	fStepper = new G4ClassicalRK4( fEquationOfMotion );
     	G4cout<<"G4ClassicalRK4 is called"<<G4endl;
  }  
  else if  (aString == "RKG3_Stepper"){
    	fStepper = new G4RKG3_Stepper( fEquationOfMotion );
     	G4cout<<"G4RKG3_Stepper is called"<<G4endl;
  }
  else if  (aString == "CashKarpRKF45"){
    	fStepper = new G4CashKarpRKF45( fEquationOfMotion );
     	G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  }   
  else {
    	fStepper = new G4CashKarpRKF45( fEquationOfMotion );
     	G4cout<<"The selected stepper is not available"<<G4endl;
     	G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  }
  
  if (fChordFinder) delete fChordFinder;
  G4double min_step= .000001*re;
  fChordFinder = new G4ChordFinder(this,min_step,fStepper);
  fChordFinder->SetDeltaChord(DeltaChord); 
  
  G4FieldManager* fieldMgr= 
  G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetChordFinder(fChordFinder); 

}
////////////////////////////////////////////////////////////////////////////////
//
void MAGCOSMagneticField::ComputeBfieldAtDifferentPosition(G4String cosys_pos, 
				G4double alt0, G4double dAlt, G4int nAlt,
				G4double lat0, G4double dlat, G4int nlat,
				G4double long0, G4double dlong, G4int nlong,
				G4String file_output) const
{ clock_t clock1,clock2;
  
  clock1=clock();
  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 
  std::ofstream OutputFile;
  OutputFile.open(file_output, (std::ios::binary | std::ios::out));
  OutputFile<<"System of coordinate for position: "<<cosys_pos<<std::endl;
  OutputFile<<"altitude [km]"<<'\t'<<"latitude "<<'\t'
				  <<"longitude "<<'\t'
				  <<"Bx"<<'\t'<<"By"<<'\t'
				  <<"Bz"<<'\t'<<"B"<<'\t'
				  <<"Br"<<'\t'<<"Bth"<<'\t'<<"Bphi"<<'\t'
				  <<"Bvert"<<'\t'<<"Bnorth"<<std::endl;
  OutputFile.precision(2);
  OutputFile.setf(std::ios::fixed);	
  if (cosys_pos != "GEOID") theConvertor->SetSystemInAndOut(cosys_pos,"GEO");
  for (double  altitude=alt0;altitude<alt0 + double(nAlt)*dAlt; altitude+=dAlt){
  	for (double  latitude=lat0; latitude< lat0 + double(nlat)*dlat; 
							latitude+=dlat){
		for (double  longitude=long0; longitude< long0 + double(nlong)*dlong;
							 longitude+=dlong){
			G4ThreeVector geo_pos, pos ;
			pos =G4ThreeVector(1.,0.,0.);
			pos.setRThetaPhi(Re+altitude,
						  90*degree-latitude,
						  longitude);
			geo_pos=pos;			  
			/*if (cosys_pos != "GEOID"){
				G4ThreeVector pos=G4ThreeVector(0.,0.,1.);
				pos.setRThetaPhi(Re+altitude,
						  90*degree-latitude,
						  longitude);
				geo_pos = pos;
				if (cosys_pos != "GEO")		  
					geo_pos = theConvertor->Transform(pos); 
	           	}
			else {
				theConvertor->ComputeGEOPositionFromGEOID
	            			(altitude, latitude, longitude, geo_pos);
			}*/
			G4ThreeVector Bfield =G4ThreeVector(0.,0.,0.);
  			if (Internal) Bfield = (this->*GetInternalField)(geo_pos);
  			if (External) Bfield=Bfield+(this->*GetExternalField)(geo_pos);
  			Bfield /= nanotesla;
		/*	G4double geoid_alt,geoid_long,geoid_lat;
			theConvertor->ComputeGEOIDCoordinatesFromGEOPosition
								(geo_pos,
 						     		 geoid_alt,
						     		 geoid_long,
								 geoid_lat);
			G4ThreeVector vertical = G4ThreeVector(1.,0.,0.);
  			vertical.setRThetaPhi(1.,
					      90.*degree-geoid_lat,
					      geoid_long);
  			G4ThreeVector east_dir =G4ThreeVector(-std::sin(longitude),std::cos(longitude),0. );
  			G4ThreeVector north_dir = -east_dir.cross(vertical);
  			G4double Bdown = -Bfield.dot(vertical);
  			G4double Bnorth = Bfield.dot(north_dir);
  			G4double Beast = Bfield.dot(east_dir);
			
			vertical = geo_pos/geo_pos.mag();
			G4double Br= Bfield.dot(vertical);
			G4ThreeVector theta_dir = east_dir.cross(vertical);
			G4double Bth = Bfield.dot(theta_dir);*/
			
			
		/*	OutputFile<<altitude/km<<'\t'<<'\t'
				  <<latitude/degree<<'\t'<<'\t'
				  <<longitude/degree<<'\t'<<'\t'
				  <<Bfield.x()<<'\t'<<Bfield.y()<<'\t'
				  <<Bfield.z()<<'\t'<<Bfield.mag()<<'\t'
				  <<Br<<'\t'<<Bth<<'\t'<<Beast<<'\t'
				  <<Bdown<<'\t'<<Bnorth<<std::endl; */           	
		}
	}
  
  }
  OutputFile.close();
  
  clock2=clock();
  double tclock=double(clock2-clock1)/(double(CLOCKS_PER_SEC));
  G4cout<<"time used for the cutoff computation: "<<tclock<<" s"<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
//
std::vector< double> MAGCOSMagneticField::ComputeMcIlwainLParameter(double pos[],G4double pitch_angle)
{
 bool DidIntegrationSucceed; 
 G4String OldExtFieldModelName = ExtFieldModelName;
 
 //only internal field
 
 SetExternalField(G4String("NOFIELD"));
 G4double I_invar= theSecondInvariantCalculator->ComputeI(DidIntegrationSucceed, pos,pitch_angle);
 
 //Compute L from I and Bm by using the Hilton approximation JGR 76, 28,6952-6954,1971 
 G4double Bm= theSecondInvariantCalculator->GetBm();
 G4double Bm_B0=Bm/DipoleB0;
 G4double X1_3=I_invar*std::pow(Bm_B0,1./3.)/Re;
 G4double X2_3=X1_3*X1_3;
 G4double X=X1_3*X2_3;
 G4double F= 1.+1.35047*X1_3+0.465376*X2_3+0.0475455*X;
 G4double L=std::pow(F/Bm_B0,1./3.);
 //G4cout<<"L "<<L<<std::endl;
 
 std::vector< double> Lvec;
 Lvec.clear(); 
 Lvec.push_back(L);
 Lvec.push_back(L);
 
 if (OldExtFieldModelName == "NOFIELD") return Lvec;
 
 // compute L with external field
 
 SetExternalField(OldExtFieldModelName);
 I_invar= theSecondInvariantCalculator->ComputeI(DidIntegrationSucceed, pos,pitch_angle);
 //Compute L from I and Bm by using the Hilton approximation JGR 76, 28,6952-6954,1971 
 Bm= theSecondInvariantCalculator->GetBm();
 Bm_B0=Bm/DipoleB0;
 X1_3=I_invar*std::pow(Bm_B0,1./3.)/Re;
 X2_3=X1_3*X1_3;
 X=X1_3*X2_3;
 F= 1.+1.35047*X1_3+0.465376*X2_3+0.0475455*X;
 L=std::pow(F/Bm_B0,1./3.);
 Lvec[1]=L;
 return Lvec;
 
}
////////////////////////////////////////////////////////////////////////////////
//
std::vector< double> MAGCOSMagneticField::ComputeMcIlwainLParameter(G4ThreeVector pos,G4double pitch_angle) 
{double pos1[3];
 pos1[0]= pos.x();
 pos1[1]= pos.y();
 pos1[2]= pos.z(); 
 return ComputeMcIlwainLParameter(pos1,pitch_angle);
 
}
