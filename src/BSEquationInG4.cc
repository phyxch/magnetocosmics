#include "BSEquationInG4.hh"
#include "MAGCOSUnits.hh"
////////////////////////////////////////////////////////////////////////////////
//
BSEquationInG4::BSEquationInG4(G4EquationOfMotion* anEquation, G4Navigator* aNavigator)
{ 
  fEquation=anEquation;
  fNavigator=aNavigator;
  nvar=7; 
  xmax=re;
  hinit=10000000000000.;
  crossing_delta= .0001*re;
  BSIntegrator::getInstance()->SetEquationTobeIntegrated(this);
 
  SafetyAtNewCheckPosition=0.;
  SafetyAtLastPosition=0.;
  SafetyAtLastSubPosition=0.;
}
////////////////////////////////////////////////////////////////////////////////
//
BSEquationInG4::BSEquationInG4()
{ 
  fEquation=0;
  fNavigator=0;
  nvar=7;
  xmax=re;
  hinit=100000000000000.;
  crossing_delta= .0001*re;
  SafetyAtNewCheckPosition=0.;
  SafetyAtLastPosition=0.;
  SafetyAtLastSubPosition=0.;  
  BSIntegrator::getInstance()->SetEquationTobeIntegrated(this);   
}
////////////////////////////////////////////////////////////////////////////////
//
BSEquationInG4::~BSEquationInG4()
{ 
  if (fEquation) delete fEquation;
  if (fNavigator) delete fNavigator;
}
////////////////////////////////////////////////////////////////////////////////
//
void BSEquationInG4::derivs(const G4double x, const G4double  y[], G4double dydx[])
{ 
  // lorentz equation where y[0],y[1],y[2] = position
  //                        y[3],y[4],y[5] = momentum given in MeV unit
  //                         y[6]    = for crossing check 
 
  G4ThreeVector Position=G4ThreeVector(y[0],y[1],y[2]);
  double y_vec[6],dy_vec[6];
  
  //lorentz equation
  for (G4int i=0;i<6;i++) y_vec[i]=y[i];
  fEquation->RightHandSide(  y_vec, dy_vec);
  for (G4int i=0;i<6;i++) dydx[i]=dy_vec[i];
 
  //conditions for crossing boundary
  // dydx[6] =1 till the boundary of ther physical volume has been encountered 
  
  G4ThreeVector vector1;
  vector1=Position-NewCheckPosition;
  if (x == 0){ 
    dydx[6]=1.; //at x=0 the point is in its own volume
    DidCrossBoundaryDuringMmid=false;
  } 
  else if (DidCrossBoundaryDuringMmid){ 
    dydx[6]=0.;
  }    
  
  else if (vector1.mag()<SafetyAtNewCheckPosition ) dydx[6]=1.; //in this case it
  // is sure that the particle did not encounter the boundary
  
  else{ 
    if (!WasSafetyAtLastPositionAlreadyComputed){
      
      ///G4cout<<"test compute new saftey"<<std::endl;
      fNavigator->LocateGlobalPointWithinVolume(LastPosition);
      SafetyAtLastPosition=fNavigator->ComputeSafety(LastPosition);
      SafetyAtNewCheckPosition=SafetyAtLastPosition;
      NewCheckPosition=LastPosition;
      /*G4cout<<LastPosition.mag()/(6371200.*m)<<endl;
	G4cout<<" Safety "<<SafetyAtNewCheckPosition/(6371200.*m)<<endl;
	G4cout<<SafetyAtLastPosition/(6371200.*m)<<endl;*/
      WasSafetyAtLastPositionAlreadyComputed=true;
    }
    
    vector1= Position-LastPosition;
    
    if (vector1.mag()<SafetyAtLastPosition) dydx[6]=1.;
    else {
      fNavigator->LocateGlobalPointWithinVolume(LastSubPosition);
      //SafetyAtLastSubPosition=fNavigator->ComputeSafety(LastSubPosition);
      vector1= Position-LastSubPosition;
      double vec1_mag=vector1.mag();
      if  (vec1_mag == 0.){ 
	dydx[6]=1.;
      }
      else {	
	vector1=vector1/vec1_mag; 
	double linear_length= 
	  fNavigator->ComputeStep(LastSubPosition, 
				  vector1,vec1_mag,
				  SafetyAtLastSubPosition);
	
	
	if (linear_length >= vec1_mag )  dydx[6]=1.;
	else {
	  dydx[6]=0.;
	  //G4cout<<" here "<<std::endl;
	  //G4cout<<linear_length/re<<" "<<vec1_mag/re<<std::endl;
	  DidCrossBoundaryDuringMmid=true;
	  DidCrossBoundaryDuringBsstep=true;
	}
	
      } 
      
      
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//
int BSEquationInG4::StopCondition(const G4double x, G4double y[])
{ 
  G4ThreeVector Position=G4ThreeVector(y[0],y[1],y[2]);
  /*G4cout<<"test"<<endl;
  G4cout<<Position.mag()/(6371.2*km)<<endl;
  G4cout<<NewCheckPosition.mag()/(6371.2*km)<<endl;
  G4cout<<LastPosition.mag()/(6371.2*km)<<endl;
  G4cout<<LastSubPosition.mag()/(6371.2*km)<<endl;
  G4cout.precision(20);
  G4cout<<"LastX "<<LastX/(6371.2*km)<<endl;
  G4cout<<x/(6371.2*km)<<endl;*/
 
  /*G4cout<<"Stop test"<<endl;

  G4cout<<Position.mag()/(6371.2*km)<<endl;
  G4cout<<y[6]/(6371.2*km)<<endl;
  G4cout<<(x-y[6])/(6371.2*km)<<endl;
 
  G4cout<<" DidCrossBoundaryDuringBsstep "
       <<DidCrossBoundaryDuringBsstep<<std::endl;
  G4cout<<" DidCrossBoundaryDuringMmid "
       <<DidCrossBoundaryDuringMmid<<std::endl;  */    
  //Check if cross boundary
  //----------------------------------
  G4ThreeVector vector1;
  vector1=Position-NewCheckPosition;
  if (vector1.mag() >= SafetyAtNewCheckPosition) {
    if (!WasSafetyAtLastPositionAlreadyComputed) {
      fNavigator->LocateGlobalPointWithinVolume(LastPosition);
      SafetyAtLastPosition=fNavigator->ComputeSafety(LastPosition);
      SafetyAtNewCheckPosition=SafetyAtLastPosition;
      NewCheckPosition=LastPosition;
      WasSafetyAtLastPositionAlreadyComputed=true;
    }
    
    vector1= Position-LastPosition;
    fNavigator->LocateGlobalPointWithinVolume(LastPosition);
    if (vector1.mag()>=SafetyAtLastPosition) {
      double vec1_mag=vector1.mag();
      if  (vec1_mag == 0.){ 
	G4cout<<vec1_mag<<" vec1_mag "<<G4endl;
	return 1;
      }	
      vector1=vector1/vec1_mag; 
      double linear_length= fNavigator->ComputeStep
	(LastPosition, vector1,vec1_mag,SafetyAtLastPosition);
      
      if (linear_length <= vec1_mag ){ 
	if (std::abs(linear_length-vec1_mag) >= crossing_delta){
	  G4double dx=x-LastX;
	  if (DidCrossBoundaryDuringBsstep)
	    h_before_step = y[6] - LastX;
	  else h_before_step = dx*linear_length/vec1_mag;
	  if (x>xmax) {
	    h_before_step=std::min(h_before_step,xmax-LastX);
	  }
	  return 0;
	}
	else if (x <= xmax) {
	  return 1;
	}
      }   
    }	 
  }    
  //integration will  continue
  //---------------------------     
  if ( !DidCrossBoundaryDuringBsstep) {
    if (x == xmax) return  1;
    else if (x >xmax) {
      h_before_step=xmax-LastX;
      return 0;
    }
    y[6]=x;
    return 2;
  }
  n_crossed_the_boundary+=1;
  h_before_step=std::min(y[6],xmax)-LastX;
  return 0;
} 

////////////////////////////////////////////////////////////////////////////////
//
bool BSEquationInG4::InterruptCondition(const G4double , G4double *)
{ 
  return false;  
}
////////////////////////////////////////////////////////////////////////////////
//
double BSEquationInG4::StepBeforeStop(const G4double , const G4double *,
				      const G4double x, G4double *)
{ 
  //G4cout<<h_before_step/(6371.2*km)<<endl;
  if (h_before_step<= 0 )
    h_before_step=.5*(std::min(x,xmax)-LastX); 
  return h_before_step;
}
////////////////////////////////////////////////////////////////////////////////
//
double BSEquationInG4::ComputeTrajectory(G4FieldTrack& pFieldTrack,
					 double smax, double new_safety)
{
  BSIntegrator* theIntegrator = BSIntegrator::getInstance();
  if ( theIntegrator->GetKmax() < 2)  theIntegrator->SetKmax(2);
  theIntegrator->SetEquationTobeIntegrated(this);
  theIntegrator->SetIsTheEquationOfMotionForG4Tracking(true);
  theIntegrator->SetXmax(smax);
  
  xmax=smax;
  
  G4double y0[7];
  
  G4ThreeVector StartPosition, StartMomentum;
  StartPosition = pFieldTrack.GetPosition(); 
  // <<StartPosition.mag()/(6371.2*km)<<endl;
  
  StartMomentum= pFieldTrack.GetMomentum();
  
  SafetyAtNewCheckPosition=new_safety;
  NewCheckPosition=StartPosition;
  
  SetLastPosition(NewCheckPosition);
  SafetyAtLastPosition=new_safety;
  WasSafetyAtLastPositionAlreadyComputed=true;
  
  y0[0]=StartPosition.x();
  y0[1]=StartPosition.y();
  y0[2]=StartPosition.z();
  y0[3]=StartMomentum.x();
  y0[4]=StartMomentum.y();
  y0[5]=StartMomentum.z();
  y0[6]=0.;

  integration_precision = theIntegrator->GetEps();
 
  double sstart= pFieldTrack.GetCurveLength();
  theIntegrator->SetXmax(xmax+sstart);
  xmax+=sstart;
  int nok,nbad;
  y0[6]=sstart;
  

  // Lorentz Motion
  hinit=xmax;
  n_crossed_the_boundary=0;
  // G4cout<<"TestI1"<<std::endl;
  theIntegrator->do_integration(y0,sstart,hinit,nok,nbad);
  // G4cout<<"TestI2"<<std::endl;
  
  
  G4double slast= theIntegrator->GetXlast(); 
  slast= y0[6];
  pFieldTrack.SetPosition(G4ThreeVector(y0[0],y0[1],y0[2]));
  pFieldTrack.SetMomentum(G4ThreeVector(y0[3],y0[4],y0[5]));                                      
  pFieldTrack.SetCurveLength(slast);
  // G4cout<<slast<<" "<<sstart<<endl;
  return slast-sstart;
}
