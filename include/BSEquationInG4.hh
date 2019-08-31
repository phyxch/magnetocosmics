// 9/18/2014: Hexc, Olesya - Check on the code and clean up the code.
#ifndef BSEquationInG4_h
#define BSEquationInG4_h 1

#include "BSEquation.hh"
#include "BSIntegrator.hh"
#include "G4EquationOfMotion.hh"
#include "globals.hh"
#include "G4Navigator.hh"
#include "G4FieldTrack.hh"

class BSEquationInG4: public BSEquation
{
public:
  
  BSEquationInG4(G4EquationOfMotion* anEquation, G4Navigator* aNavigator);
  BSEquationInG4();
  ~BSEquationInG4();
  
  
  void derivs(const G4double x, const G4double  y[], G4double dydx[]);
  int StopCondition(const G4double x, G4double y[]);
  bool InterruptCondition(const G4double x, G4double y[]);
  double StepBeforeStop(const G4double xlast, const G4double ylast[],
                                                   const G4double x, G4double y[]);
 
  //Set and Get method
  //------------------
  inline void SetXmax(double val) {xmax=val;}
  inline double GetXmax() {return xmax;}
  inline void  SetCrossingDelta(double val) {crossing_delta=val;}
  
  inline void Sethinit(double val) {hinit=val;}
  inline void SetDirectionOfIntegration(double val) 
  {direction_of_integration=val;}
  
  inline void SetEquationOfMotion(G4EquationOfMotion* anEquation )
  {fEquation = anEquation;}
  inline void SetNavigator(G4Navigator* aNavigator )
  {fNavigator = aNavigator;}
  
  inline void SetDidCrossBoundaryDuringMmid(G4bool aBool)
  {DidCrossBoundaryDuringMmid =aBool;}
  
  inline void SetDidCrossBoundaryDuringBsstep(G4bool aBool)
  {DidCrossBoundaryDuringBsstep =aBool;}			  			    			    
  inline G4EquationOfMotion* GetEquationOfMotion( )
  {return fEquation;}	
  
  inline G4bool GetDidCrossBoundaryDuringMmid( )
  {return DidCrossBoundaryDuringMmid;}
  
  inline G4bool GetDidCrossBoundaryDuringBsstep( )
  {return DidCrossBoundaryDuringBsstep;}			  
  
  // method to compute the trajectory		      
  
  double ComputeTrajectory(G4FieldTrack& pFieldTrack,double xmax,double
			   new_safety);
  
private:
  double xmax,hinit; 
  double direction_of_integration;
  
  G4EquationOfMotion* fEquation; 
  //Navigator used to detect boundary
  G4Navigator   *fNavigator; 
  
  // Safety Variable
  //-----------------
  double SafetyAtNewCheckPosition;
  double SafetyAtLastPosition;
  double SafetyAtLastSubPosition;
   
  //Entering exiting 
  //-------------
  G4bool IsEnteringDaughterOrExitingMother;
   
  double h_before_step;
  
  //Close to the boundary 
  //-------------
  G4int n_crossed_the_boundary;
  
  // crossing precision
  //-------------
  G4double crossing_delta;
  
  //precision
  //-------------
  G4double integration_precision;
  
  //check of Boundary crossing
  //---------------------------------
  G4bool DidCrossBoundaryDuringMmid;
  G4bool DidCrossBoundaryDuringBsstep;
};
#endif
