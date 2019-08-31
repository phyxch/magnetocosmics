#ifndef MAGCOSEQUATIONOFMOTION
#define MAGCOSEQUATIONOFMOTION 

#include "G4Mag_EqRhs.hh"
#include "G4MagneticField.hh"
#include "G4String.hh"

// DESCRIPTION
// -----------
//
// This class defines the equation of motion used in the simulation. 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// MAGCOSEquationOfMotion(G4MagneticField* MagField )
//    Constructor.
//
// ~MAGCOSEquationOfMotion()
//    Destructor.
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

class MAGCOSEquationOfMotion : public G4Mag_EqRhs
{
   public:  // with description

     // Constructor and destructor.
     MAGCOSEquationOfMotion( G4MagneticField* MagField );
    ~MAGCOSEquationOfMotion() {;}
    
     // Given the value of the magnetic field B, this function 
     // calculates the value of the derivative dydx.

     void EvaluateRhsGivenB( const G4double y[],
			     const G4double B[3],
				   G4double dydx[] ) const;
     
     // this is a pointer on one of the equation of motion methods
     
     void (MAGCOSEquationOfMotion::*Equation)( const G4double y[],
			                       const G4double B[3],
				               G4double dydx[] ) const;	
				   
     //Lorentz equation of motion				   
     void LorentzMotion( const G4double y[],
			  const G4double B[3],
		          G4double dydx[] ) const;

     //Equation for tracing magnetic field line			  
     void MotionAlongBfieldLine( const G4double y[],
			         const G4double B[3],
				 G4double dydx[]) const;			   			   		   
     
     //Set methods
     
     void SetReverseTimeMode(G4bool abool);  
     void SetEquationType(G4String aString);
     
    
     
    private:
     G4bool reverse_time;
     G4double direction;
     
     
};

#endif /* G4MAG_USUAL_EQRHS */
