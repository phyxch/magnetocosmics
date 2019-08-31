// 9/18/2014: Hexc, Olesya - Check on the code and clean up the code.
#ifndef BSEquation_h
#define BSEquation_h 1

#include "G4ThreeVector.hh"

class BSEquation
{
public:
  BSEquation();
  BSEquation(int nval);
  virtual ~BSEquation()=0;
  
  //method
public:
  
  virtual void derivs(const G4double x, const G4double  y[], G4double dydx[]) =0;
  virtual int StopCondition(const G4double x, G4double y[]) =0;
  virtual bool InterruptCondition(const G4double x, G4double y[]) =0;
  virtual double StepBeforeStop(const G4double xlast, const G4double ylast[],
				const G4double x, G4double y[]) =0;
  
  //inline methods
  //--------------
  inline int GetNvar() {return nvar;}
  
  inline  void SetLastPosition(G4ThreeVector aVec){
    LastPosition=aVec;
    LastSubPosition=aVec;
    WasSafetyAtLastPositionAlreadyComputed=false;
  }
  
  inline  void SetLastSubPosition(G4ThreeVector aVec) {
    LastSubPosition=aVec;
  }
  inline void SetLastSubX(G4double aVal){
    LastSubX=aVal;
  }
  inline void SetLastX(G4double aVal){LastX=aVal;
    LastSubX=aVal;
  }
  
protected:
  int nvar;
  
  //Position
   //-----------------------
  G4ThreeVector LastPosition , NewCheckPosition;
  G4ThreeVector LastSubPosition ;
  G4double LastSubX,LastX;
   
  G4bool WasSafetyAtLastPositionAlreadyComputed;  
  
};
#endif
