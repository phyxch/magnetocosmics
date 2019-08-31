// 9/17/2014: Hexc & Olesya:  Checking on the code
//
#include "SecondInvariant.hh"
#include "MAGCOSUnits.hh"
#include "MAGCOSMagneticField.hh"
////////////////////////////////////////////////////////////////////////////////
//
SecondInvariant::SecondInvariant(MAGCOSMagneticField* aField)
{ 
  theMagneticField=aField;
  nvar=5; 
  Bm_relative_precision=.0000000001;
  top_of_atmosphere=0.;
}
////////////////////////////////////////////////////////////////////////////////
//
SecondInvariant::~SecondInvariant()
{ 
}
////////////////////////////////////////////////////////////////////////////////
//
void SecondInvariant::derivs(const G4double , const G4double  y[], G4double dydx[])
{ 
  double pos[3];
  pos[0]=y[0];
  pos[1]=y[1];
  pos[2]=y[2];
  
  
  G4double B[3]; 
  theMagneticField->GetFieldValue(pos,B);
  double BB=std::sqrt( B[0] * B[0] + B[1] * B[1] + B[2] * B[2]); 

  dydx[0]=direction*B[0]/BB;
  dydx[1]=direction*B[1]/BB;
  dydx[2]=direction*B[2]/BB;
  dydx[3]=0;
  
  if ( NeedtoComputeI && BB <= Bm ) dydx[3]=direction*std::sqrt(1-BB/Bm); 
  
  dydx[4] =1.;
  if (BB > Bm) dydx[4]=0.;

  
}
////////////////////////////////////////////////////////////////////////////////
//
int SecondInvariant::StopCondition(const G4double x, G4double y[])
{ return StopAtMirrorPoint(x,y); 
}
////////////////////////////////////////////////////////////////////////////////
//
int SecondInvariant::StopAtMirrorPoint(const G4double , G4double y[])
{ double pos[3];
  pos[0]=y[0];
  pos[1]=y[1];
  pos[2]=y[2];
  /*G4cout<<"pos "<<pos[0]/re<<'\t'<<pos[1]/re<<'\t'<<pos[2]/re<<std::endl;
  G4cout<<"Bm "<<Bm/nanotesla<<std::endl;*/
   
  
  G4double B[3]; 
  theMagneticField->GetFieldValue(pos,B);
  Bnew=std::sqrt( B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
  //G4cout<<"Bnew "<<Bnew/nanotesla<<std::endl;
  G4double* ylast= BSIntegrator1::getInstance()->GetYlast();
  double h=y[4]-ylast[4];
  G4double deltaB_rel=(Bm-Bnew)/Bm;
 /* G4cout<<"deltaB_rel "<<deltaB_rel<<std::endl;
  G4cout<<"y[3] "<<y[3]<<std::endl;
  G4cout<<"h "<<h/re<<std::endl;*/
  if (y[3] == 0 &&  deltaB_rel <0) return 1;
  if (  std::abs(deltaB_rel) <= Bm_relative_precision && (y[3] != 0.)
                                    && (std::abs(y[4])>.00001 *re)) return 1;
  if ( deltaB_rel <0 && (y[3] != 0.) && (std::abs(h/Re) <1e-9)) return 1;				    
  if ( deltaB_rel <0 && (y[3] != 0.)) return 0;
  Blast=Bnew;
  return 2;				    
} 
////////////////////////////////////////////////////////////////////////////////
//
bool SecondInvariant::InterruptCondition(const G4double, G4double y[])
{ 
  G4double r=sqrt(y[0]*y[0] + y[1]*y[1] + y[2] * y[2]);
  G4double re= 6371200.*m; 
  G4double rmin=re/20.;
  //G4cout<<r/rmin<<std::endl;
  if( r<rmin){
    LastInterruptCondition=true;
    return true;
  }
  else { 
    LastInterruptCondition=false;
    return false;
  }  
}

////////////////////////////////////////////////////////////////////////////////
//
double SecondInvariant::StepBeforeStop(const G4double , const G4double ylast[],
				       const G4double, G4double y[])
{ 
  double h=y[4]-ylast[4];
  // G4cout<<"h "<<h<<endl;
  if (h==0.) h=1.e-9*km;
  return h;
}

////////////////////////////////////////////////////////////////////////////////
//
G4double SecondInvariant::ComputeI(double pos[], bool& DidIntegrationSucceed)
{ 
  // Test is B< or > Bm
   //--------------------
   
   G4double B[3]; 
   theMagneticField->GetFieldValue(pos,B);
   double BB=std::sqrt( B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
   if (BB>Bm) {
   	G4cout<<"The start position is at a position with B>Bm."<<std::endl;
	G4cout<<"No second invariant integration will be computed"<<std::endl;
   	DidIntegrationSucceed = false;
	return -1;
   }
   NeedtoComputeI = true;
   
   BSIntegrator1* theIntegrator = BSIntegrator1::getInstance();
   G4int  max_step = theIntegrator->GetMaxstep();
   G4double  eps =  theIntegrator->GetEps();
   G4double  hmax =  theIntegrator->GetHmax();
   
   theIntegrator->SetMaxstep(20000);
   theIntegrator->SetHmax(5.*Re);
   theIntegrator->SetEps(1.e-6);
   
   if ( theIntegrator->GetKmax() < 2)  theIntegrator->SetKmax(2);
   theIntegrator->SetEquationTobeIntegrated(this);
   theIntegrator->SetIsTheEquationOfMotionForG4Tracking(false);
   double ystart[5];
   ystart[0]=pos[0];
   ystart[1]=pos[1];
   ystart[2]=pos[2];
   ystart[3]=0.;
   ystart[4]=0.;
   
   double sstart=0.;
   int nok,nbad;
   direction = 1.;
   DidIntegrationSucceed = theIntegrator->do_integration(ystart,sstart,10000000000.,nok,nbad);
   double* ylast = theIntegrator->GetYlast();
   G4double I= ylast[3];
   if (DidIntegrationSucceed) DidIntegrationSucceed = !LastInterruptCondition;
   //south  part of the line 
   ystart[0]=pos[0];
   ystart[1]=pos[1];
   ystart[2]=pos[2];
   ystart[3]=0.;
   ystart[4]=0.;
   sstart=0.;
   direction = -1.;
   G4bool DidIntegrationSucceed1 = theIntegrator->do_integration(ystart,sstart,10000000000.,nok,nbad);
   I+=-ylast[3];
   if (DidIntegrationSucceed1) DidIntegrationSucceed1=!LastInterruptCondition; 
   DidIntegrationSucceed =  DidIntegrationSucceed && DidIntegrationSucceed1; 
   
   
   theIntegrator->SetMaxstep(max_step);
   theIntegrator->SetHmax(hmax);
   theIntegrator->SetEps(eps);
   
   
   return std::abs(I); 
   
       
}
////////////////////////////////////////////////////////////////////////////////
//
G4double SecondInvariant::ComputeI(bool& DidIntegrationSucceed, double pos[],double alpha)
{// compute Bm 
   G4double B[3];
   theMagneticField->GetFieldValue(pos,B);
   double BB=std::sqrt( B[0] * B[0] + B[1] * B[1] + B[2] * B[2]); 
   double sina=std::sin(alpha);
   Bm=BB/(sina*sina);
   return   ComputeI(pos, DidIntegrationSucceed);
   

}
