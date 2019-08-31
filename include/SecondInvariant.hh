// 9/17/2014: Hexc, Oleya - clean up the code
#ifndef SecondInvariant_h
#define SecondInvariant_h 1

#include "BSEquation.hh"
#include "BSIntegrator1.hh"
#include "globals.hh"
#include "G4Navigator.hh"
#include "G4FieldTrack.hh"

class MAGCOSMagneticField;
class SecondInvariant: public BSEquation
{
public:
  
  SecondInvariant(MAGCOSMagneticField* aField);
  ~SecondInvariant();
  
  void derivs(const G4double x, const G4double  y[], G4double dydx[]);
  int StopCondition(const G4double x, G4double y[]);
  int StopAtMirrorPoint(const G4double x, G4double y[]);
  bool InterruptCondition(const G4double x, G4double y[]);
  double StepBeforeStop(const G4double xlast, const G4double ylast[],
			const G4double x, G4double y[]);
 
  //Set and Get method
  //------------------
  inline void SetXmax(double val) {xmax=val;}
  inline double GetXmax() {return xmax;}
  inline void Sethinit(double val) {hinit=val;} 
  inline void SetXstart( double val){xstart= val;}
  inline void SetBm_relative_precision(double val) {Bm_relative_precision=val;}
  inline void SetBm(double aVal){Bm=aVal;} 
  inline double GetBm(){return Bm;} 			  
  
  // methods to compute  the mirror point and I, and to find Bmin Position		      
  G4double ComputeI(double pos[], bool& DidIntegrationSucceed);
  G4double ComputeI(bool& DidIntegrationSucceed, double pos[], double alpha);
	      
 private:
  double xstart,B0,Bm,alpha0,top_of_atmosphere,Blast,Bnew; 
  double Bm_relative_precision;
  MAGCOSMagneticField* theMagneticField;
  double NeedtoComputeI;
  bool LastInterruptCondition;
  G4bool InterruptCalledFromEquation;
  double xmax,hinit; 
  double direction;
  
  double h_before_step;
  
  // crossing precision
  //-------------
  G4double crossing_delta;
  
  //precision
  //-------------
  G4double integration_precision;
  
};
#endif
